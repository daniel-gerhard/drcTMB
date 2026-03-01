#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level, default is 50
#' @param random logical, if TRUE, the variance-covariance matrix of the random effects is included in the standard error estimation
#' @param residual logical, if TRUE, the residual variance is included in the standard error estimation
#'
#' @return a vector with estimates
#' @export
#'
ed <- function(x, respLev=50, random=TRUE, residual=FALSE){ 
  p <- respLev/100
  est <- x$sdrl$par.fixed
  covest <- x$sdrl$cov.fixed
  rn <- names(est)
  id <- rn %in% paste("b", 1:5, sep="")
  est <- est[id,drop=FALSE]
  covest <- covest[id,id,drop=FALSE]
  
  pv <- numeric(length=5)
  pv[!x$fix] <- est
  pv[x$fix] <- x$start[x$fix]
  pvar <- matrix(0, nrow=5, ncol=5)
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

  cvc <- bdiag(pvar, rvc, rvar)

  b1 <- pv[1]
  b2 <- pv[2]
  b3 <- pv[3]
  b4 <- pv[4]
  b5 <- pv[5]
  u1 <- 0
  u2 <- 0
  u3 <- 0
  u4 <- 0
  u5 <- 0
  r <- 0


  # logistic
  if (x$model == "logistic"){
    edl <- eval(deriv(expression((b4+u4) + 1/(b3+u3) * log((((b1+u1)-(b2+u2))/(((b1+u1)-(b2+u2))*p - r))^(1/exp(b5 + u5)) - 1)), c("b1", "b2", "b3", "b4", "b5", "u1", "u2", "u3", "u4", "u5", "r")))
    grad <- attr(edl, "gradient")
    std <- sqrt(grad %*% cvc %*% t(grad))
  }
   
  # loglogistic   
  if (x$model == "loglogistic"){
    edl <- eval(deriv(expression((b4+u4) * ((((b1 + u1)-(b2 + u2))/(((b1+u1)-(b2+u2))*p - r))^(1/exp(b5 + u5)) - 1)^(1/(b3 + u3))), c("b1", "b2", "b3", "b4", "b5", "u1", "u2", "u3", "u4", "u5", "r")))
    grad <- attr(edl, "gradient")
    std <- sqrt(grad %*% cvc %*% t(grad))
  }
  
  # Weibull
  if (x$model == "weibull1"){
    cvc4 <- cvc[-c(5, 10), -c(5, 10), drop=FALSE]
    edl <- eval(deriv(expression((b4+u4) * (-1*log(p - (r / ((b1+u1) - (b2+u2)))))^(1/(b3+u3))), c("b1", "b2", "b3", "b4", "u1", "u2", "u3", "u4", "r")))
    grad <- attr(edl, "gradient")
    std <- sqrt(grad %*% cvc4 %*% t(grad))
  }
  
  if (x$model == "weibull2"){
    cvc4 <- cvc[-c(5, 10), -c(5, 10), drop=FALSE]
    edl <- eval(deriv(expression((b4+u4) * (-1*log(1 - p + (r / ((b1+u1) - (b2+u2)))))^(-1/(b3+u3))), c("b1", "b2", "b3", "b4", "u1", "u2", "u3", "u4", "r")))
    grad <- attr(edl, "gradient")
    std <- sqrt(grad %*% cvc4 %*% t(grad))
  }

  # lognormal
  if (x$model == "lognormal"){
    cvc4 <- cvc[-c(5, 10), -c(5, 10), drop=FALSE]
    edl <- b4 * exp((-1/b3)*qnorm(p - (r / (b1 - b2))))
    grad <- matrix(eval(Deriv(expression((b4+u4) * exp((-1/(b3+u3))*qnorm(p - (r / ((b1+u1) - (b2+u2)))))), c("b1", "b2", "b3", "b4", "u1", "u2", "u3", "u4", "r"))), nrow=1)
    std <- sqrt(grad %*% cvc4 %*% t(grad))
  }

  res <- c(edl, as.vector(std))
  names(res) <- c("Estimate", "Std. Error")
  return(res)
}
