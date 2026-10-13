#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level, default is 50
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

  if (!is.null(x$rform) & random == TRUE){
    tn <- rownames(ests)
    thetas <- ests[tn %in% "theta",,drop=FALSE]
    corrs <- ests[tn %in% "corr",,drop=FALSE]
    cr <- matrix(corrs[,1], nrow=sqrt(nrow(corrs)))
    sdm <- exp(thetas[1:sqrt(nrow(corrs))])
    rname <- all.vars(x$rform[[2]])
    if (length(sdm) == 1) sdm <- as.matrix(sdm)
    vc <- diag(sdm) %*% cr %*% diag(sdm)
    idr <- paste("b", 1:5, sep="") %in% rname
    rvc <- matrix(0, nrow=5, ncol=5)
    rvc[idr, idr] <- vc
  } else {
    rvc <- matrix(0, nrow=5, ncol=5)
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
  grad <- matrix(0, nrow=length(linfct), ncol=length(pv) + 6)

  for (i in seq_along(linfct)){
    lpv <- as.vector(linfct[[i]] %*% pv)

    if ((lpv[2] > lpv[1] & lpv[3] < 0) | (lpv[2] < lpv[1] & lpv[3] > 0)){
      p <- respLev[i]/100
    } else {
      p <- 1 - respLev[i]/100
    }

    # ED and its gradient w.r.t. (b1-b5, u1-u5, r) from the drcED TMB template
    eobj <- MakeADFun(list(model="drcED", mod=mod, lnk=lnk, p=p),
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
  res <- list(estimates=edl, variance=edvar)
  class(res) <- "edTMB"
  return(res)
}
