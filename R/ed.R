#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level
#'
#' @return a vector with estimates
#' @export
#'
ed <- function(x, respLev){ 
  p <- respLev
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
  
  # logistic
  if (x$model == "logistic"){
    expr2 <- 100/p
    expr4 <- expr2^(1/pv[5])
    expr5 <- expr4 - 1
    expr6 <- log(expr5)
    value <- pv[4] + expr6/-pv[3]
    grad <- array(0, c(length(value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
    grad[, "b"] <- -(expr6/pv[3]^2)
    grad[, "c"] <- 0
    grad[, "d"] <- 0
    grad[, "e"] <- 1
    grad[, "f"] <- -(expr4 * (log(expr2) * (1/pv[5]^2))/expr5/-pv[3])
    std <- sqrt(grad %*% pvar %*% t(grad)) 
  }
   
  # loglogistic   
  if (x$model == "loglogistic"){
    tempVal <- log((100 - p)/100)
    value <- pv[4] * (exp(-tempVal/pv[5]) - 1)^(1/-pv[3])
    grad <- value * rbind(c(-log(exp(-tempVal/pv[5]) - 1)/(-pv[3]^2), 
                     0, 0, 1/pv[4], 
                     exp(-tempVal/pv[5]) * tempVal/(pv[5]^2) * (1/-pv[3]) * ((exp(-tempVal/pv[5]) - 1)^(-1))))
    std <- sqrt(grad %*% pvar %*% t(grad)) 
  }
  
  # Weibull
  if (x$model == "weibull1" || x$model == "weibull2"){
    pvw <- pvar[1:4, 1:4]
    if (x$model == "weibull1" && pv[3] > 0) p <- 100 - p
    tempVal <- log(-log((100 - p)/100))
    value <- exp(tempVal/pv[3] + log(pv[4]))
    grad <- value * rbind(c(-tempVal/(pv[3]^2), 0, 0, 1/pv[4]))
    std <- sqrt(grad %*% pvw %*% t(grad)) 
  }
  
  # lognormal
  if (x$model == "lognormal"){
    pvw <- pvar[1:4, 1:4]
    p <- 100 - p
    pProp <- 1 - (100 - p)/100
    expr2 <- exp(qnorm(pProp)/-pv[3])
    value <- pv[4] * expr2
    grad <- array(0, c(length(value), 4L), list(NULL, c("b", "c", "d", "e")))
    grad[, "b"] <- -(pv[4] * (expr2 * (qnorm(pProp)/(-pv[3]^2))))
    grad[, "c"] <- 0
    grad[, "d"] <- 0
    grad[, "e"] <- expr2
    std <- sqrt(grad %*% pvw %*% t(grad)) 
  }
  return(c(value, std))
}
