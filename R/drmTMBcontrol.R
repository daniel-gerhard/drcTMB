#' Set control parameters for drcTMB
#'
#' @param optimiser character string specifying the optimiser to use
#' @param iter.max maximum number of iterations of the optimiser, see \code{\link{nlminb}}.
#' @param eval.max maximum number of evaluations of the objective function, see \code{\link{nlminb}}.
#'
#' @return a list with control parameters
#' @export
#'

drmTMBcontrol <- function(optimiser="nlminb", iter.max=1000, eval.max=2000){
    list(optimiser=optimiser, iter.max=iter.max, eval.max=eval.max)
}
