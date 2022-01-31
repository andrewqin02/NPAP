make_syn<-function(index,
                   foo,
                   predictors,
                   predictors.op,
                   dependent,
                   unit.variable,
                   time.variable,
                   special.predictors,
                   treated.units,
                   control.units,
                   time.predictors.prior,
                   time.optimize.ssr,
                   unit.names.variable,
                   time.plot,
                   Sigf.ipop){
  
  treated.units <- treated.units[index]
  
  dataprep.out<-Synth::dataprep(
    foo = foo,
    predictors = predictors,
    predictors.op = predictors.op,
    dependent = dependent,
    unit.variable = unit.variable,
    time.variable = time.variable,
    special.predictors = special.predictors,
    treatment.identifier = treated.units,
    controls.identifier = control.units,
    time.predictors.prior = time.predictors.prior,
    time.optimize.ssr = time.optimize.ssr,
    unit.names.variable = unit.names.variable,
    time.plot = time.plot)
  
  synth.out<-Synth::synth(data.prep.obj = dataprep.out, optimxmethod = "All")
  
  a<-data.frame(dataprep.out$Y0plot %*% synth.out$solution.w)
  
  colnames(a)<-paste('unit',index,sep='')
  
  out<-list(a = a, synth.out=synth.out, dataprep.out = dataprep.out)
  
  return(out)
}

# Function used in generate,placebos - edit to expand optimization method
# 
syn_plac <- function(i, 
                     dataprep.out, 
                     Sigf.ipop, 
                     tr, 
                     names.and.numbers) {
  cases <- as.numeric(dataprep.out$tag$controls.identifier)
  X0 <- dataprep.out$X0[, -i]
  X1 <- matrix(dataprep.out$X0[, i, drop = F])
  Z0 <- dataprep.out$Z0[, -i]
  Z1 <- matrix(dataprep.out$Z0[, i, drop = F])
  Y0plot <- dataprep.out$Y0plot[, -i]
  Y1plot <- dataprep.out$Y0plot[, i, drop = F]
  foo <- character(length = length(dataprep.out$tag$foo))
  id <- as.numeric(dataprep.out$tag$unit.variable)
  d <-
    which(regexpr(tr, stringr::str_split(dataprep.out$tag$foo[id],
                                         ', ')[[1]]) == T)
  for (j in 1:length(dataprep.out$tag$foo)) {
    foo[j] <-
      paste(stringr::str_split(dataprep.out$tag$foo[j],
                               ', ')[[1]][-c(d[1]:d[length(d)])], 
            collapse =', ')
  }
  predictors <- dataprep.out$tag$predictors
  predictors.op <- dataprep.out$tag$predictors.op
  special.predictors <- dataprep.out$tag$special.predictors
  dependent <- dataprep.out$tag$dependent
  unit.variable <- dataprep.out$tag$unit.variable
  time.variable <- dataprep.out$tag$time.variable
  treatment.identifier <- i
  controls.identifier <- dataprep.out$tag$controls.identifier[-i]
  time.predictors.prior <- dataprep.out$tag$time.predictors.prior
  time.optimize.ssr <- dataprep.out$tag$time.optimize.ssr
  time.plot <- dataprep.out$tag$time.plot
  unit.names.variable <- dataprep.out$tag$unit.names.variable
  tag <- list(
    foo = as.character(foo),
    predictors = predictors,
    predictors.op = predictors.op,
    special.predictors = special.predictors,
    dependent = dependent,
    unit.variable = unit.variable,
    time.variable = time.variable,
    treatment.identifier = treatment.identifier,
    controls.identifier = controls.identifier,
    time.predictors.prior = time.predictors.prior,
    time.optimize.ssr = time.optimize.ssr,
    time.plot = time.plot,
    unit.names.variable = unit.names.variable
  )
  dp <- list(
    X0 = X0,
    X1 = X1,
    Z0 = Z0,
    Z1 = Z1,
    Y0plot = Y0plot,
    Y1plot = Y1plot,
    names.and.numbers = names.and.numbers,
    tag = tag
  )
  s.out <- Synth::synth(data.prep.obj = dp, 
                        Sigf.ipop = Sigf.ipop, optimxmethod = "All")
  a <- data.frame(dp$Y0plot %*% s.out$solution.w)
  s.mspe <- s.out$loss.v
  res <- list(a = a, s.mspe = s.mspe)
  return(res)
}

synth_object_names <- c("solution.v",
                        "solution.w",
                        "loss.v",
                        "loss.w",
                        "custom.v",
                        "rgV.optim")

dataprep_object_names <- c( "X0", 
                            "X1", 
                            "Z0",
                            "Z1",
                            "Y0plot",
                            "Y1plot",
                            "names.and.numbers",
                            "tag" )
generate.placebos <- function(dataprep.out,
                              synth.out,
                              Sigf.ipop = 5,
                              strategy = "sequential") {
  
  # Inputs Checkss
  if(!all(names(dataprep.out) %in% dataprep_object_names)){
    stop("Invalid dataprep object. Have you run `dataprep` on your data?")
  }
  
  if(!all(names(synth.out) %in% synth_object_names)){
    stop("Invalid synth output. Have you run the `synth` function?")
  }
  
  if(is.integer(Sigf.ipop) | Sigf.ipop <= 0){
    stop("You have not passed a valid argument for Signf.ipop. Please pass a positive integer")
  }
  
  strategy_match <- match.arg(strategy, c("sequential", "multiprocess"))
  
  unit.numbers <- NULL
  
  tr <- as.numeric(dataprep.out$tag$treatment.identifier)
  
  names.and.numbers <-
    subset(dataprep.out$names.and.numbers, unit.numbers != tr)
  
  n <- length(dataprep.out$tag$controls.identifier)
  
  b <-
    data.frame(matrix(
      0,
      ncol = n,
      nrow = length(dataprep.out$tag$time.plot)
    ))
  
  #mspe.placs <- data.frame(matrix(0, ncol = 1, nrow = n))
  
  # for (i in 1:n) {
  #   temp <- syn_plac(i, dataprep.out, Sigf.ipop, tr, names.and.numbers)
  #   b[, i] <- temp$a
  #   colnames(b)[i] <-
  #     paste('synthetic', as.character(names.and.numbers[i, 2]), sep = '.')
  #   mspe.placs[i, ] <- temp$s.mspe
  # }
  
  strategy <- paste0("future::", strategy_match)
  
  plan(strategy)
  oplan <- plan()
  
  mspe2 <- furrr::future_map(1:n, ~syn_plac(.x, 
                                            dataprep.out, 
                                            Sigf.ipop, 
                                            tr, names.and.numbers))
  
  b <- dplyr::bind_cols(purrr::map(mspe2,"a"))
  b <- setNames(b, paste('synthetic', as.character(names.and.numbers[ ,2]), sep = '.'))
  
  mspe.placs <- purrr::map(mspe2, "s.mspe")
  mspe.placs <- as.data.frame(unlist(mspe.placs))
  
  
  on.exit(plan(oplan), add = TRUE)
  
  df <-
    cbind(
      b,
      dataprep.out$Y0,
      dataprep.out$Y1,
      dataprep.out$Y0plot %*% synth.out$solution.w,
      dataprep.out$tag$time.plot
    )
  colnames(df)[(ncol(df) - 2):ncol(df)] <-c('Y1', 'synthetic.Y1', 'year')
  
  t0 <- as.numeric(dataprep.out$tag$time.optimize.ssr[1])
  t1 <-
    as.numeric(dataprep.out$tag$time.optimize.ssr[length(dataprep.out$tag$time.optimize.ssr)]) +
    1
  treated.name <-
    as.character(dataprep.out$names.and.numbers$unit.names[dataprep.out$names.and.numbers[, 2] %in% dataprep.out$tag$treatment.identifier])
  
  loss.v <- synth.out$loss.v
  
  res2 <-
    list(
      df = df,
      mspe.placs = mspe.placs,
      t0 = t0,
      t1 = t1,
      tr = tr,
      names.and.numbers = names.and.numbers,
      n = n,
      treated.name = treated.name,
      loss.v = loss.v
    )
  
  class(res2) <- append(class(res2),"tdf")
  
  return(res2)
}