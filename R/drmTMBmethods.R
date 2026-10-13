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
  if (!is.null(x$rform)){
    thetas <- ests[rn %in% "theta",,drop=FALSE]
    corrs <- ests[rn %in% "corr",,drop=FALSE]
    cr <- matrix(corrs[,1], nrow=sqrt(nrow(corrs)))
    sdm <- exp(thetas[1:sqrt(nrow(corrs))])
    names(sdm) <- rownames(cr) <- colnames(cr) <- all.vars(x$rform[[2]])
    cat("\nStd.Dev.:\n")
    print(sdm)
    cat("\nCorrelation:\n")
    print(cr)
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
    rformrh <- rformlh <- object$rform
    rformrh[2] <- NULL
    rf <- model.frame(subbars(rformrh), data=newdata)
    rt <- mkReTrms(findbars(object$rform), rf)
    Z <- t(rt$Zt)
    re <- ests[rn %in% "u",,drop=FALSE]
    rformlh[3] <- NULL
    rcoefs <- all.vars(rformlh)
    rind <- which(paste("b", 1:5, sep="") %in% rcoefs)
    umat <- matrix(re[,1], ncol=length(rind))
    pxl <- lapply(1:5, function(i){
      if (i %in% rind){
        pxl[[i]] <- pxl[[i]] + Z %*% umat[,which(i == rind), drop=FALSE]
      } else {
        return(pxl[[i]]) 
      }
    })
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
