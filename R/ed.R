#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level in percent, default is 50. The effective dose ED p is the dose at which the response has changed by p percent of the difference between the asymptote at low doses (the control) and the asymptote at high doses.
#' @param linfct a matrix of linear functions to indentify the parameters of the dose-response curve (each row vector results in parameters b1-b5), default is the identity matrix
#' @param random logical, if TRUE, the variance-covariance matrix of the random effects is included in the standard error estimation
#' @param residual logical, if TRUE, the residual variance is included in the standard error estimation
#'
#' @return a edTMB object
#' @export
#'
ed <- function(x, respLev=50, linfct=diag(5), random=TRUE, residual=FALSE){ 
  if (!is.list(linfct)) linfct <- list(linfct)
  if (length(respLev) != length(linfct)){
    if (length(respLev) == 1){
      respLev <- rep(respLev, length(linfct))
    }
    if (length(linfct) == 1){
      linfct <- rep(linfct, length(respLev))
    } 
    if (length(respLev) != length(linfct)){
      stop("respLev and linfct must have the same length or one of them must be of length 1.")
    }
  }
  if (is.null(names(linfct))) names(linfct) <- paste("lf", 1:length(linfct), sep="")
  
  est <- x$sdrl$par.fixed
  covest <- x$sdrl$cov.fixed
  rn <- names(est)
  id <- rn %in% paste("b", 1:5, sep="")
  est <- est[id,drop=FALSE]
  covest <- covest[id,id,drop=FALSE]
  
  pv <- numeric(length=length(x$start))
  pv[!x$fix] <- est
  pv[x$fix] <- x$start[x$fix]
  pvar <- matrix(0, nrow=length(x$start), ncol=length(x$start))
  pvar[!x$fix, !x$fix] <- covest

  ests <- x$ssdrl

  # random effects: sum of the covariance matrices of all terms
  rvc <- matrix(0, nrow=5, ncol=5)
  if (!is.null(x$rform) & random == TRUE){
    for (vc in recov(x)){
      rvc[vc$pars, vc$pars] <- rvc[vc$pars, vc$pars] + vc$cov
    }
  }

  if (residual == TRUE){
      rvar <- exp(x$sdrl$par.fixed["log_sigma"])^2
  } else {
      rvar <- 0
  }

  mod <- which(x$model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
  lnk <- which(x$link == c("identity", "log", "logit"))

  # gradient of each ED w.r.t. (parameters, u1-u5, r)
  edl <- numeric(length=length(linfct))
  ql <- numeric(length=length(linfct))
  grad <- matrix(0, nrow=length(linfct), ncol=length(pv) + 6)

  for (i in seq_along(linfct)){
    lpv <- as.vector(linfct[[i]] %*% pv)

    # level of the increasing function F of the model at the ED: the control is b1 for an
    # increasing curve (b3 > 0) and for the models without b3, and b2 for a decreasing curve
    q <- respLev[i]/100
    if (!(mod %in% c(6, 7, 8)) && lpv[3] < 0) q <- 1 - q
    ql[i] <- q

    # ED and its gradient w.r.t. (b1-b5, u1-u5, r) from the drcED TMB template
    eobj <- MakeADFun(list(model="drcED", mod=mod, lnk=lnk, q=q),
                      list(b=lpv, u=rep(0, 5), r=0),
                      ADreport=TRUE, DLL="drcTMB_TMBExports", silent=TRUE)
    edl[i] <- eobj$fn()
    egrad <- eobj$gr()
    grad[i, ] <- c(egrad[1:5] %*% linfct[[i]], egrad[6:11])
  }

  cvc <- as.matrix(bdiag(pvar, rvc, rvar))
  edvar <- grad %*% cvc %*% t(grad)
  names(edl) <- paste(respLev, names(linfct), sep=":")
  colnames(edvar) <- rownames(edvar) <- paste(respLev, names(linfct), sep=":")
  # the fit, linear functions and levels of F are kept for profile
  names(ql) <- names(edl)
  res <- list(estimates=edl, variance=edvar, fit=x, linfct=linfct, q=ql)
  class(res) <- "edTMB"
  return(res)
}
