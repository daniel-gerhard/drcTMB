#' edTMB Methods
#'
#' @param x,object an edTMB object
#' @param ... anything else
#'
#' @return Estimated effective dose coefficients
#' @export
#'
#' @name edTMB-methods
#' @aliases print.edTMB coef.edTMB vcov.edTMB

#' @rdname edTMB-methods
print.edTMB <- function(x, ...){
    ests <- x$estimates
    vc <- x$variance
    std <- sqrt(diag(vc))
    rdat <- data.frame(Estimate=ests, Std.Error=std)
    cat("Effective dose estimates:\n\n")
    print(rdat, digits=3)
}

#' @rdname edTMB-methods
coef.edTMB <- function(object, ...){
    return(object$estimates)
}

#' @rdname edTMB-methods
vcov.edTMB <- function(object, ...){
    return(object$variance)
}
