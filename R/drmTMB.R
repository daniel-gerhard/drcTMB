#' Fit a drmTMB model
#' 
#' @param form formula with name of response and dose scale
#' @param fform list with fixed-effect right-hand side formulas for each model parameter
#' @param rform random effect formula of the form b1 + b2 + b3 + b4 + b5 ~ 1 | cluster
#' @param data a data.frame object
#' @param start a vector with starting values for each model parameter. If NA a self-starting function is called.
#' @param fix a binary vector with same length as start. If TRUE, a model coefficient is fixed to the starting value.
#' @param lower vector of lower bounds, passed to nlminb.
#' @param upper vector of upper bounds, passed to nlminb.
#' @param family Specify the distributional assumptions. One of "gaussian", "binomial", and "poisson".
#' @param model Specify the nonlinear model. One of "logistic", "loglogistic", "weibull1", "weibull2", "lognormal".
#'
#' @return A drmTMB object
#' @useDynLib drmTMB

drmTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, lower=-Inf, upper=Inf, family="gaussian", model="L5"){
  # response
  mf <- model.frame(form, data=data)
  # fixed effects
  if (is.null(fform)){
    fform <- vector("list", 5)
    names(fform) <- paste("X", 1:5, sep="")
    for (i in 1:5) fform[[i]] <- ~ 1
  } else {
    names(fform) <- paste("X", 1:5, sep="")
  }
  Xs <- lapply(fform, function(x) model.matrix(x, data=data))
  npars <- sapply(Xs, ncol)
  npi <- rep(1:length(npars), npars)
  if (is.null(start)) start <- rep(NA, sum(npars))
  if (is.null(fix)) fix <- rep(FALSE, sum(npars))
  
  # model
  mod <- which(model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal"))
  
  # random effects
  if (!is.null(rform)){
    rformrh <- rformlh <- rform
    rformrh[2] <- NULL
    rf <- model.frame(subbars(rformrh), data=data)
    rt <- mkReTrms(findbars(rform), rf)
    Z <- t(rt$Zt)
    
    # ... allow user-specific term names
    rformlh[3] <- NULL
    rcoefs <- all.vars(rformlh)
    ind <- as.numeric(paste("b", 1:5, sep="") %in% rcoefs)
  }
  # response
  y <- mr <- model.response(mf, "any")
  if (family == "binomial" & ncol(as.matrix(mr)) == 2){
    bn <- mr[,1] + mr[,2]
    y <- mr[,1]
  } else {
    bn <- rep(0, length(y))  
  }
  
  # data list for TMB
  if (is.null(rform)){
    dlist <- list(y=y, x=mf[,2], bn=bn, mod=mod)
  } else {
    dlist <- list(y=y, x=mf[,2], bn=bn, mod=mod, ind=ind, z0=rep(0, nrow(data)), Z=Z)
  }
  datalist <- c(dlist, Xs)
  
  # Starting values
  sb1 <- start[which(npi == 1)]
  if (is.na(sb1[1]) & !fix[which(npi == 1)[1]]) sb1[1] <- mean(dlist$y[which(dlist$x == max(dlist$x))])
  sb1[is.na(sb1)] <- 0
  sb2 <- start[which(npi == 2)]
  if (is.na(sb2[1]) & !fix[which(npi == 2)[1]]) sb2[1] <- mean(dlist$y[which(dlist$x == min(dlist$x))])
  sb2[is.na(sb2)] <- 0
  sb3 <- start[which(npi == 3)]
  if (is.na(sb3[1]) & !fix[which(npi == 3)[1]]) sb3[1] <- -coefficients(lm(dlist$y ~ dlist$x))[2]
  sb3[is.na(sb3)] <- 0
  sb4 <- start[which(npi == 4)]
  if (is.na(sb4[1]) & !fix[which(npi == 4)[1]]) sb4[1] <- min(dlist$x[rank(dlist$y)/length(dlist$y) > 0.5])
  sb4[is.na(sb4)] <- 0
  sb5 <- start[which(npi == 5)]
  if (is.na(sb5[1]) & !fix[which(npi == 5)[1]]) sb5[1] <- 1
  sb5[is.na(sb5)] <- 0
  
  if (is.null(rform)){
    parameters <- list(b1=sb1, b2=sb2, b3=sb3, b4=sb4, b5=sb5)
    if (family == "gaussian") parameters <- c(parameters, list(log_sigma=1))
  } else {
    nvc <- sum(ind > 0)
    m <- ncol(Z)
    vtheta <- rep(0.001, nvc)
    ctheta <- rep(0, nvc*(nvc-1)/2)
    stheta <- c(vtheta, ctheta)
    up <- matrix(rnorm(nvc*m, 0, rep(vtheta, each=m)), ncol=nvc)
    parameters <- list(b1=sb1, b2=sb2, b3=sb3, b4=sb4, b5=sb5, u=as.vector(up), theta=stheta)
    if (family == "gaussian") parameters <- c(parameters, list(log_sigma=0))
  }
  
  # fix parameters to starting values
  if (!any(fix)){
    if (is.null(rform)){
      if (family == "gaussian") obj <- MakeADFun(c(model = "drcnormfix", datalist), parameters, DLL="drcTMB_TMBExports")
      if (family == "poisson") obj <- MakeADFun(c(model = "drcpoisfix", datalist), parameters, DLL="drcTMB_TMBExports")
      if (family == "binomial") obj <- MakeADFun(c(model = "drcbinomfix", datalist), parameters, DLL="drcTMB_TMBExports")
    } else {
      if (family == "gaussian") obj <- MakeADFun(c(model = "drcnorm", datalist), parameters, random="u", DLL="drcTMB_TMBExports")
      if (family == "poisson") obj <- MakeADFun(c(model = "drcpois", datalist), parameters, random="u", DLL="drcTMB_TMBExports")
      if (family == "binomial") obj <- MakeADFun(c(model = "drcbinom", datalist), parameters, random="u", DLL="drcTMB_TMBExports")
    }
  } else {
    lmap <- vector("list", 5)
    names(lmap) <- paste("b", 1:5, sep="")
    idp <- 1:sum(npars) 
    idp[fix] <- NA
    for (i in 1:5){
      lmap[[i]] <- idp[npi == i]  
    }
    lmap <- lapply(lmap, function(x) as.factor(x))
    if (is.null(rform)){
      if (family == "gaussian") obj <- MakeADFun(c(model = "drcnormfix", datalist), parameters, DLL="drcTMB_TMBExports", map=lmap)
      if (family == "poisson") obj <- MakeADFun(c(model = "drcpoisfix", datalist), parameters, DLL="drcTMB_TMBExports", map=lmap)
      if (family == "binomial") obj <- MakeADFun(c(model = "drcbinomfix", datalist), parameters, DLL="drcTMB_TMBExports", map=lmap)
    } else {
      if (family == "gaussian") obj <- MakeADFun(c(model = "drcnorm", datalist), parameters, random="u", DLL="drcTMB_TMBExports", map=lmap)
      if (family == "poisson") obj <- MakeADFun(c(model = "drcpois", datalist), parameters, random="u", DLL="drcTMB_TMBExports", map=lmap)
      if (family == "binomial") obj <- MakeADFun(c(model = "drcbinom", datalist), parameters, random="u", DLL="drcTMB_TMBExports", map=lmap)
    }
  }
  obj$hessian <- TRUE
  #opt <- do.call("optim", obj)
  if (is.null(rform)){
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper)
  } else {
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper)
  }
  res <- list()
  res$family <- family
  res$model <- model
  res$start <- start
  res$fix <- fix
  res$Xs <- Xs
  res$rform <- rform
  res$estimates <- opt$par
  res$opt <- opt
  res$obj <- obj
  res$sdrl <- sdreport(obj)
  res$ssdrl <- summary(res$sdrl)
  class(res) <- "drmTMB"
  return(res)
}