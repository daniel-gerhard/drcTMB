#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib drcTMB, .registration = TRUE
#' @useDynLib drcTMB_TMBExports
#' @importFrom stats AIC BIC approx coef coefficients confint fitted lm lm.fit logLik
#' @importFrom stats model.frame model.matrix model.response nlminb nobs pbeta pbinom
#' @importFrom stats plogis pnorm ppois predict printCoefmat profile qchisq qlogis qnorm
#' @importFrom stats quantile residuals runif sd setNames simulate splinefun vcov
#' @importFrom stats complete.cases terms .getXlevels na.pass
#' @importFrom graphics plot
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix t bdiag
#' @importFrom reformulas subbars mkReTrms findbars
#' @importFrom TMB MakeADFun sdreport
## usethis namespace: end
NULL
