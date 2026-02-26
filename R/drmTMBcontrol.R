#' Set control parameters for drcTMB
#' 
#' @param optimiser character string specifying the optimiser to use
#'
#' @return a list with control parameters
#' @export
#' 

drmTMBcontrol <- function(optimiser="nlminb"){
    list(optimiser=optimiser)
}
