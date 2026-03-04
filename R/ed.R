#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level, default is 50
#' @param linfct a matrix of linear functions to indentify the parameters of the dose-response curve (each row vector results in parameters b1-b5), default is the identity matrix
#' @param random logical, if TRUE, the variance-covariance matrix of the random effects is included in the standard error estimation
#' @param residual logical, if TRUE, the residual variance is included in the standard error estimation
#'
#' @return a vector with estimates
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

  edl <- numeric(length=length(linfct))
  if (x$model %in% c("logistic", "loglogistic")){
    grad <- matrix(0, nrow=length(linfct), ncol=11)
  } else {
    grad <- matrix(0, nrow=length(linfct), ncol=9)
  }
  

  for (i in seq_along(linfct)){   
    p <- respLev[i]/100
    lpv <- as.vector(linfct[[i]] %*% pv)
    lpvar <- linfct[[i]] %*% pvar %*% t(linfct[[i]])

    cvc <- bdiag(lpvar, rvc, rvar)
    cvc4 <- cvc[-c(5, 10), -c(5, 10), drop=FALSE]

    b1 <- lpv[1]
    b2 <- lpv[2]
    b3 <- lpv[3]
    b4 <- lpv[4]
    b5 <- lpv[5]
    u1 <- 0
    u2 <- 0
    u3 <- 0
    u4 <- 0
    u5 <- 0
    r <- 0


    # logistic
    if (x$model == "logistic"){
      edls <- eval(deriv(expression((b4+u4) + 1/(b3+u3) * log((((b1+u1)-(b2+u2))/(((b1+u1)-(b2+u2))*p - r))^(1/exp(b5 + u5)) - 1)), c("b1", "b2", "b3", "b4", "b5", "u1", "u2", "u3", "u4", "u5", "r")))
      edl[i] <- edls
      grad[i, ] <- as.vector(attr(edls, "gradient"))
    }
   
    # loglogistic   
    if (x$model == "loglogistic"){
      edls <- eval(deriv(expression((b4+u4) * ((((b1 + u1)-(b2 + u2))/(((b1+u1)-(b2+u2))*p - r))^(1/exp(b5 + u5)) - 1)^(1/(b3 + u3))), c("b1", "b2", "b3", "b4", "b5", "u1", "u2", "u3", "u4", "u5", "r")))
      edl[i] <- edls
      grad[i, ] <- as.vector(attr(edls, "gradient"))
    }
  
    # Weibull
    if (x$model == "weibull1"){    
      edls <- eval(deriv(expression((b4+u4) * (-1*log(p - (r / ((b1+u1) - (b2+u2)))))^(1/(b3+u3))), c("b1", "b2", "b3", "b4", "u1", "u2", "u3", "u4", "r")))
      edl[i] <- edls
      grad[i, ] <- as.vector(attr(edls, "gradient"))
    }
  
    if (x$model == "weibull2"){
      edls <- eval(deriv(expression((b4+u4) * (-1*log(1 - p + (r / ((b1+u1) - (b2+u2)))))^(-1/(b3+u3))), c("b1", "b2", "b3", "b4", "u1", "u2", "u3", "u4", "r")))
      edl[i] <- edls
      grad[i, ] <- as.vector(attr(edls, "gradient"))
    }

    # lognormal
    if (x$model == "lognormal"){
      edl[i] <- b4 * exp((-1/b3)*qnorm(p - (r / (b1 - b2))))
      grad[i, ] <- eval(Deriv(expression((b4+u4) * exp((-1/(b3+u3))*qnorm(p - (r / ((b1+u1) - (b2+u2)))))), c("b1", "b2", "b3", "b4", "u1", "u2", "u3", "u4", "r")))
    }
  }

  if (x$model %in% c("logistic", "loglogistic")){
    edvar <- grad %*% cvc %*% t(grad)
  } else {
    edvar <- grad %*% cvc4 %*% t(grad)
  }
  names(edl) <- paste(respLev, names(linfct), sep=":")
  colnames(edvar) <- rownames(edvar) <- paste(respLev, names(linfct), sep=":")
  res <- list(estimates=edl, variance=edvar)
  return(res)
}
