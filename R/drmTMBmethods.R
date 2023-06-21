#' Title
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
  lsigma <- ests[rn %in% "log_sigma",,drop=FALSE]
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
