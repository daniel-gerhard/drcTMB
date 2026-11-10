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
    cat(if (identical(x$type, "peak")) "Peak of the hormesis curves, dose and response:\n\n"
        else if (identical(x$type, "bmd")) paste("Benchmark doses, ", c(excess="excess risk", additional="additional risk", relative="relative deviation", sd="standard deviations")[[x$def]],
                                                 if (length(x$fit$reterms) > 0) paste(",", x$scale, "BMD"), ":\n\n", sep="")
        else "Effective dose estimates:\n\n")
    print(rdat, digits=3)
    invisible(x)
}

#' @rdname edTMB-methods
#' @export
coef.edTMB <- function(object, ...){
    return(object$estimates)
}

#' @rdname edTMB-methods
#' @export
vcov.edTMB <- function(object, ...){
    return(object$variance)
}
