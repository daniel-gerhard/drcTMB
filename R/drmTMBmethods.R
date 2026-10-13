#' drmTMB Methods
#'
#' @param x,object a drmTMB object
#' @param newdata a data frame for prediction. If NULL, the original data is used.
#' @param random logical, if TRUE, the random effects are included in the predictions. For simulate, if TRUE, new random effects are drawn for each simulation, otherwise the fitted random effects are kept.
#' @param nsim number of response vectors to simulate
#' @param seed an optional seed for the random number generator
#' @param ... anything else
#'
#' @return Estimated model coefficients, variance-covariance matrix, or predictions.
#' simulate returns a data frame with nsim columns of simulated responses at the estimated
#' parameters; for a binomial response given as cbind(successes, failures), the number of successes.
#' @export
#'
#' @name drmTMB-methods
#' @aliases print.drmTMB predict.drmTMB coef.drmTMB vcov.drmTMB simulate.drmTMB
#' 
#' @rdname drmTMB-methods
print.drmTMB <- function(x, ...){
  ests <- x$ssdrl
  rn <- rownames(ests)
  coefs <- ests[rn %in% paste("b", 1:5, sep=""),,drop=FALSE]
  #lsigma <- ests[rn %in% "log_sigma",,drop=FALSE]
  cat("Coefficients:\n")
  print(coefs)
  for (vc in recov(x)){
    pn <- paste("b", vc$pars, sep="")
    names(vc$sd) <- rownames(vc$corr) <- colnames(vc$corr) <- pn
    cat("\nRandom effects:", vc$group, "\nStd.Dev.:\n")
    print(vc$sd)
    if (length(pn) > 1){
      cat("Correlation:\n")
      print(vc$corr)
    }
  }
}

#' @rdname drmTMB-methods
predict.drmTMB <- function(object, newdata=NULL, random=TRUE){
  if (is.null(newdata)) newdata <- object$data
  form <- object$formula
  form[2] <- NULL
  mf <- model.frame(form, data=newdata)
  x <- mf[,1]
  Xs <- lapply(object$fform, function(x) model.matrix(x, data=newdata))
  npars <- sapply(Xs, ncol)
  npi <- rep(1:length(npars), npars)
  ests <- object$ssdrl
  rn <- rownames(ests)
  coefs <- ests[rn %in% paste("b", 1:5, sep=""),,drop=FALSE]
  cest <- object$start
  cest[!object$fix] <- coefs[,1]
  pxl <- lapply(1:5, function(i) Xs[[i]] %*% cest[npi == i])
  
  if (random && !is.null(object$rform)){
    # random effects of the fitted levels; levels not in the fitted data get zero
    re <- ests[rn %in% "u", 1]
    uoff <- 0
    for (tm in object$reterms){
      nl <- length(tm$levels)
      np <- length(tm$pars)
      U <- matrix(re[uoff + seq_len(nl*np)], nrow=nl)
      uoff <- uoff + nl*np
      rt <- mkReTrms(findbars(tm$formula), model.frame(subbars(tm$formula), data=newdata))
      Zt <- rt$Ztlist[[tm$name]]
      lid <- match(rownames(Zt), tm$levels)
      Znew <- t(Zt[!is.na(lid), , drop=FALSE])
      for (k in seq_len(np)){
        pxl[[tm$pars[k]]] <- pxl[[tm$pars[k]]] + as.vector(Znew %*% U[lid[!is.na(lid)], k])
      }
    }
  }

  if (object$link == "log"){
    pxl[[1]] <- exp(pxl[[1]])
    pxl[[2]] <- exp(pxl[[2]])
  }
  if (object$link == "logit"){
    pxl[[1]] <- 1/(1 + exp(-1*pxl[[1]]))
    pxl[[2]] <- 1/(1 + exp(-1*pxl[[2]]))
  }
  
  if (object$model == "logistic"){
    flx <- 1 + exp(pxl[[3]]*(x - pxl[[4]]))
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / (flx^exp(pxl[[5]]))
  }
  if (object$model == "loglogistic"){
    flx <- 1 + exp(pxl[[3]]*(log(x) - log(pxl[[4]])))
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / (flx^exp(pxl[[5]]))
  }  
  if (object$model == "weibull1"){
    flx <- exp(-1*exp(pxl[[3]]*(log(x) - log(pxl[[4]]))))
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) * flx
  } 
  if (object$model == "weibull2"){
    flx <- 1 - exp(-1*exp(pxl[[3]]*(log(x) - log(pxl[[4]]))))
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) * flx
  }   
  if (object$model == "lognormal"){
    flx <- pxl[[3]]*(log(x) - log(pxl[[4]]))
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) * pnorm(flx)
  }
  if (object$model == "michaelismenten"){
    flx <- 1 + pxl[[4]] / x
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / flx
  }
  if (object$model == "asymptoticregression"){
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / (1 - exp(-x / pxl[[4]]))
  }
  if (object$model == "exponentialdecay"){
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / exp(-x / pxl[[4]])
  }
  if (object$model == "gompertz"){
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) * exp(-exp(pxl[[3]] * (x - pxl[[4]])))
  }
  return(as.vector(fx))
}

#' @rdname drmTMB-methods
coef.drmTMB <- function(object, ...){
    return(object$estimates)
}

#' @rdname drmTMB-methods
vcov.drmTMB <- function(object, ...){
    return(object$sdrl$cov.fixed)
}

#' @rdname drmTMB-methods
simulate.drmTMB <- function(object, nsim=1, seed=NULL, random=TRUE, ...){
  if (!is.null(seed)) set.seed(seed)
  obj <- object$obj
  par <- obj$env$last.par.best
  if (!is.null(object$rform) && !random){
    # keep the fitted random effects: rebuild the objective with simre = 0
    dat <- obj$env$data
    dat$simre <- 0
    obj <- MakeADFun(dat, obj$env$parList(par=par), map=obj$env$map, random="u", DLL="drcTMB_TMBExports", silent=TRUE)
  }
  sims <- lapply(seq_len(nsim), function(i) obj$simulate(par)$y)
  nbad <- sum(sapply(sims, function(y) any(!is.finite(y))))
  if (nbad > 0) warning(nbad, " of ", nsim, " simulations contain non-finite values; simulated random effects may fall outside the valid parameter range (e.g., b4 <= 0).")
  sims <- as.data.frame(sims, col.names=paste("sim", seq_len(nsim), sep="_"))
  attr(sims, "seed") <- seed
  sims
}


# standard deviations, correlation and covariance matrix of each random effect term
recov <- function(object){
  ests <- object$ssdrl
  rn <- rownames(ests)
  thetas <- ests[rn %in% "theta", 1]
  corrs <- ests[rn %in% "corr", 1]
  toff <- 0
  coff <- 0
  lapply(object$reterms, function(tm){
    np <- length(tm$pars)
    sdm <- exp(thetas[toff + seq_len(np)])
    cr <- matrix(corrs[coff + seq_len(np*np)], nrow=np)
    toff <<- toff + np + np*(np-1)/2
    coff <<- coff + np*np
    list(pars=tm$pars, group=tm$group, sd=sdm, corr=cr, cov=diag(sdm, nrow=np) %*% cr %*% diag(sdm, nrow=np))
  })
}
