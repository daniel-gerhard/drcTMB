#' drmTMB Methods
#'
#' @param x a drmTMB object
#' @param ... anything else
#'
#' @return Estimated model coefficients
#' @export
#'
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
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / (flx^pxl[[5]])
  }
  if (object$model == "loglogistic"){
    flx <- 1 + exp(pxl[[3]]*(log(x) - log(pxl[[4]])))
    fx <- pxl[[2]] + (pxl[[1]] - pxl[[2]]) / (flx^pxl[[5]])
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
  return(as.vector(fx))
}

