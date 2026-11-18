#' @details The main functions are:
#' \itemize{
#'   \item \code{\link{drmTMB}}: fit a dose-response model; methods in \code{\link{drmTMB-methods}},
#'     diagnostics with \code{\link{plot.drmTMB}};
#'   \item \code{\link{ed}}, \code{\link{peak}}: effective doses and the peak of hormesis curves;
#'   \item \code{\link{bmd}}, \code{\link{bmdma}}: benchmark doses and their model averaging;
#'   \item \code{\link{bmdbayes}}, \code{\link{bmdefsa}}: Bayesian model averaging of benchmark doses;
#'   \item \code{\link{confint.drmTMB}}, \code{\link{profile.drmTMB}}, \code{\link{bootstrap}},
#'     \code{\link{vcovrobust}}: intervals and standard errors;
#'   \item \code{\link{anova.drmTMB}}: model comparison.
#' }
#' \code{vignette("drcTMB", package = "drcTMB")} introduces the package with examples.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib drcTMB, .registration = TRUE
#' @useDynLib drcTMB_TMBExports
#' @importFrom stats AIC BIC approx coef coefficients confint cor fitted lm lm.fit lm.wfit logLik
#' @importFrom stats model.frame model.matrix model.response nlminb optim nobs pbeta pbinom
#' @importFrom stats dbeta dlnorm dnorm dt plnorm pt
#' @importFrom stats plogis pnorm ppois predict printCoefmat profile qchisq qlogis qnorm rlnorm
#' @importFrom stats quantile residuals runif sd setNames simulate splinefun vcov
#' @importFrom stats complete.cases terms .getXlevels na.pass as.formula
#' @importFrom stats rnorm rpois rbinom rbeta rnbinom pnbinom anova pchisq pf lowess ppoints qt qf df.residual uniroot optimise pgamma qgamma cov.wt rchisq ave dist median
#' @importFrom graphics plot
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix t bdiag
#' @importFrom reformulas subbars mkReTrms findbars
#' @importFrom TMB MakeADFun sdreport
#' @importFrom methods as
## usethis namespace: end
NULL
