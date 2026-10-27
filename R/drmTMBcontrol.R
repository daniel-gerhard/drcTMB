#' Set control parameters for drcTMB
#'
#' @param optimiser character string specifying the optimiser to use
#' @param iter.max maximum number of iterations of the optimiser, see \code{\link{nlminb}}.
#' @param eval.max maximum number of evaluations of the objective function, see \code{\link{nlminb}}.
#' @param smallsample if TRUE (default), gaussian models without random effects use the unbiased
#'   residual variance with n - p residual degrees of freedom and t quantiles, as lm, nls and drc; if
#'   FALSE, the maximum likelihood estimate and normal quantiles. See \code{\link{drmTMB}}.
#'
#' @return a list with control parameters
#' @export
#'

drmTMBcontrol <- function(optimiser="nlminb", iter.max=1000, eval.max=2000, smallsample=TRUE){
    list(optimiser=optimiser, iter.max=iter.max, eval.max=eval.max, smallsample=smallsample)
}
