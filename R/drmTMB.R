#' Fit a drmTMB model
#'
#' @param form formula with name of response and dose scale
#' @param fform list with fixed-effect right-hand side formulas for each model parameter
#' @param rform random effect formula of the form b1 + b2 + b3 + b4 + b5 ~ 1 | cluster
#' @param data a data.frame object
#' @param start a vector with starting values for each model parameter. If NA a self-starting function is called.
#' @param fix a binary vector with same length as start. If TRUE, a model coefficient is fixed to the starting value.
#' @param weights a vector of weights for each observation.
#' @param lower vector of lower bounds, passed to nlminb.
#' @param upper vector of upper bounds, passed to nlminb.
#' @param family Specify the distributional assumptions. One of "gaussian", "binomial", "beta", and "poisson".
#' @param model Specify the nonlinear model. One of "logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz".
#' @param link Specify the link function for estimating asymptote parameters. One of "identity", "log", "logit". 
#' @param gsl a vector of lower bounds for a grid search. If NULL, no grid search is performed.
#' @param gsu a vector of upper bounds for a grid search. If NULL, no grid search is performed.
#' @param gsn a vector of number of grid points for a grid search. If NULL, no grid search is performed.
#' @param control a list of control parameters. See \code{\link{drmTMBcontrol}} for details.
#'
#' @return A drmTMB object
#' @useDynLib drmTMB

drmTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, weights=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", link="identity", gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  # control list
  control <- do.call(drmTMBcontrol, control)
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
  if (is.null(weights)) weights <- rep(1, nrow(data))
  
  # model
  mod <- which(model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
  lnk <- which(link == c("identity", "log", "logit"))
  fam <- which(family == c("gaussian", "poisson", "binomial", "beta"))
  
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
  
  # no random effects: empty design matrix
  if (is.null(rform)){
    ind <- rep(0, 5)
    Z <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(length(y), 0))
  }

  # data list for TMB
  dlist <- list(model="drcTMB", y=y, x=mf[,2], w=weights, bn=bn, fam=fam, mod=mod, lnk=lnk, ind=ind, Z=Z)
  datalist <- c(dlist, Xs)
  
  # Starting values
  #drcfunction <- logistic
  #drcstart <- start[c(which(npi == 3), which(npi == 1), which(npi == 2), which(npi == 4), which(npi == 5))]
  #drcstart[which(npi == 5)] <- exp(drcstart[which(npi == 5)])
  #drcstart[!fix] <- NA
  #fformdrc <- fform[c(which(npi == 3), which(npi == 1), which(npi == 2), which(npi == 4), which(npi == 5))]
  #drcoef <- coef(drm(form, fct=drcfunction(fixed=drcstart), data=data, pmodels=fformdrc))

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
  if (is.na(sb5[1]) & !fix[which(npi == 5)[1]]) sb5[1] <- 0
  sb5[is.na(sb5)] <- 0
  
  if (is.null(rform)){
    up <- numeric(0)
    stheta <- numeric(0)
  } else {
    nvc <- sum(ind > 0)
    m <- ncol(Z)
    vtheta <- rep(0.001, nvc)
    ctheta <- rep(0, nvc*(nvc-1)/2)
    stheta <- c(vtheta, ctheta)
    up <- matrix(rnorm(nvc*m, 0, rep(vtheta, each=m)), ncol=nvc)
  }
  # log_sigma: residual sd (gaussian) or precision (beta); fixed for poisson and binomial
  if (family == "gaussian") slsigma <- if (is.null(rform)) 1 else 0
  if (family == "beta") slsigma <- 10
  if (family %in% c("poisson", "binomial")) slsigma <- 0
  parameters <- list(b1=sb1, b2=sb2, b3=sb3, b4=sb4, b5=sb5, u=as.vector(up), theta=stheta, log_sigma=slsigma)

  # fix parameters to starting values
  lmap <- list()
  if (any(fix)){
    idp <- 1:sum(npars)
    idp[fix] <- NA
    for (i in 1:5){
      lmap[[paste("b", i, sep="")]] <- as.factor(idp[npi == i])
    }
  }
  if (family %in% c("poisson", "binomial")) lmap$log_sigma <- factor(NA)

  obj <- MakeADFun(datalist, parameters, map=lmap, random=if (!is.null(rform)) "u", DLL="drcTMB_TMBExports", silent=TRUE)

  ## grid search
  if (!(is.null(gsl) & is.null(gsu) & is.null(gsn))){
    if (length(gsn) == 1) gsn <- rep(gsn, length(gsl))
    if (length(gsl) != length(gsu)) stop("gsl and gsu must have the same length.")
    grid <- lapply(1:length(gsl), function(i) seq(gsl[i], gsu[i], length.out=gsn[i]))
    gridpar <- expand.grid(grid)
    objvals <- apply(gridpar, 1, function(x){
      par <- obj$par
      par[1:length(x)] <- x
      obj$fn(par)
    })
    best <- which.min(objvals)
    obj$par[1:length(gridpar[best,])] <- as.numeric(gridpar[best,])
  }

  if (is.null(rform)){
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper)
  } else {
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper)
  }
  res <- list()
  res$data <- data
  res$formula <- form
  res$family <- family
  res$model <- model
  res$link <- link
  res$start <- start
  res$fix <- fix
  res$Xs <- Xs
  res$fform <- fform
  res$rform <- rform
  res$estimates <- opt$par
  res$opt <- opt
  res$obj <- obj
  res$sdrl <- sdreport(obj, par.fixed=opt$par)
  res$ssdrl <- summary(res$sdrl)
  class(res) <- "drmTMB"
  return(res)
}
