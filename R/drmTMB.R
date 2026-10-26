#' Fit a drmTMB model
#'
#' @param form formula with name of response and dose scale
#' @param fform fixed effects for the model parameters: a formula of the form b2 + b4 ~ Stress, or a list of such formulas, e.g. list(b2 ~ Stress, b4 ~ Stress + Time). Parameters not listed have a single intercept (~ 1). A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~ 1), is also accepted.
#' @param rform random intercepts: a formula of the form b1 + b2 + b4 ~ 1 | cluster, or a list of such formulas, e.g. list(b1 + b4 ~ 1 | exper, b1 ~ 1 | plate). Random intercepts within a formula are correlated across the listed parameters; different formulas, and the levels of nested factors (1 | exper/plate), are independent terms.
#' @param data a data.frame object. Observations with missing values in the variables of the model or in the weights are left out.
#' @param start a vector with starting values for each model coefficient. NA entries are replaced by the starting values of \code{\link{drmTMBstart}}. Alternatively, a vector named by parameter, e.g. c(b1=0, b5=0), gives the same starting value to all coefficients of each named parameter.
#' @param fix a logical vector with one entry for each model coefficient. If TRUE, a model coefficient is fixed to its value in start. Alternatively, a vector named by parameter, e.g. c(b1=TRUE, b5=TRUE), fixes all coefficients of the named parameters; their values must then be given in start, e.g. start=c(b1=0, b5=0). Parameters not used by the model (b5 except for logistic and loglogistic, b3 for michaelismenten, asymptoticregression and exponentialdecay) are always fixed, at 0 unless given in start.
#' @param weights a vector of weights, one for each row of data.
#' @param lower lower bounds for the model coefficients b1-b5, a single value or a vector with the same length as start. Bounds of fixed coefficients are ignored; variance parameters are unbounded.
#' @param upper upper bounds for the model coefficients b1-b5, a single value or a vector with the same length as start.
#' @param family Specify the distributional assumptions. One of "gaussian", "binomial", "beta", and "poisson".
#' @param model Specify the nonlinear model. One of "logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz".
#' @param link Specify the link function for estimating asymptote parameters. One of "identity", "log", "logit". 
#' @param b4link link function for b4, "identity" (default) or "log", see Details. Only for the models with a positive b4, not for logistic and gompertz.
#' @param gsl a vector of lower bounds for a grid search, same length as start; NA for the default range. If gsl, gsu and gsn are all NULL, no grid search is performed. See \code{\link{drmTMBstart}}.
#' @param gsu a vector of upper bounds for a grid search, same length as start; NA for the default range.
#' @param gsn number of grid points for each parameter, a single number or a vector with the same length as start. \code{gsn=10} alone runs a grid search over the default ranges.
#' @param control a list of control parameters. See \code{\link{drmTMBcontrol}} for details.
#'
#' @details All models are of the form f(x) = b1 + (b2 - b1) F(x), with a function F increasing
#' from 0 to 1. For the models with steepness b3 (logistic, loglogistic, weibull1, weibull2,
#' lognormal, gompertz), F depends on b3 (x - b4) or b3 (log(x) - log(b4)); b1 is the lower and b2
#' the upper asymptote, and b3 is positive for an increasing and negative for a decreasing curve.
#' For michaelismenten, asymptoticregression and exponentialdecay (b3 and b5 are not used),
#' b1 is the response at dose zero and b2 the asymptote at high doses. The asymmetry parameter b5
#' of the logistic and loglogistic models is on the log scale. See the vignette for the model functions.
#'
#' With \code{b4link="log"}, the linear predictor of b4, including its random effects, is log(b4):
#' the coefficients of b4 in \code{coef}, \code{start}, \code{fix}, \code{lower}, \code{upper} and
#' \code{linfct} of \code{\link{ed}} are on the log scale. Random effects on b4 then vary b4 by a
#' factor between groups, and the estimate of b4 is the median over the groups. For the models on the
#' log-dose scale, this is the same model as the logistic (or gompertz) model of log(dose) for
#' loglogistic (weibull1), with b4 the log ED50, but it keeps observations at dose zero and reports
#' effective doses on the dose scale. Without random effects on b4, both links give the same fit.
#'
#' @return A drmTMB object
#' @export

drmTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, weights=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", link="identity", b4link="identity", gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  # control list
  control <- do.call(drmTMBcontrol, control)
  dd <- drmTMBdata(form, fform=fform, rform=rform, data=data, weights=weights, family=family, model=model, link=link, b4link=b4link)
  fixnamed <- !is.null(names(fix))
  start <- parvector(start, dd$npi, NA, "start")
  fix <- parvector(fix, dd$npi, FALSE, "fix")
  if (fixnamed && any(fix & is.na(start))){
    stop("give the values of the fixed parameters in start, e.g. start=c(b5=0), fix=c(b5=TRUE).")
  }
  fix[dd$npi %in% unusedpars(dd$mod)] <- TRUE

  # starting values for parameters with NA in start
  method <- if (is.null(gsl) & is.null(gsu) & is.null(gsn)) "lm" else "grid"
  start <- drmTMBstartvals(dd, start, fix, method=method, gsl=gsl, gsu=gsu, gsn=gsn)

  obj <- drmTMBobj(dd, start, fix, rform=rform)

  # bounds apply to the free coefficients only, not to the variance parameters
  bounds <- function(b, inf){
    if (length(b) == 1) b <- rep(b, length(start))
    if (length(b) != length(start)) stop("lower and upper must be a single value or have the same length as start.")
    c(b[!fix], rep(inf, length(obj$par) - sum(!fix)))
  }
  lower <- bounds(lower, -Inf)
  upper <- bounds(upper, Inf)

  ctrl <- list(iter.max=control$iter.max, eval.max=control$eval.max)
  if (is.null(rform)){
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper, control=ctrl)
  } else {
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper, control=ctrl)
  }
  if (!is.finite(opt$objective)){
    warning("the negative log-likelihood is not finite at the estimates; some observations may be impossible under the model, e.g. responses at dose zero with an asymptote fixed at 0.", call.=FALSE)
  } else if (opt$convergence != 0){
    warning("the optimiser did not converge: ", opt$message, call.=FALSE)
  }
  res <- list()
  res$call <- match.call()
  res$data <- dd$data
  res$formula <- form
  res$family <- family
  res$model <- model
  res$link <- link
  res$b4link <- b4link
  res$start <- start
  res$fix <- fix
  res$Xs <- dd$Xs
  res$fform <- dd$fform
  res$fterms <- dd$fterms
  res$coefnames <- coefnames(dd$Xs)
  res$rform <- rform
  res$reterms <- dd$re$terms
  res$estimates <- opt$par
  res$lower <- lower
  res$upper <- upper
  res$opt <- opt
  res$control <- control
  res$obj <- obj
  res$sdrl <- sdreport(obj, par.fixed=opt$par)
  res$ssdrl <- summary(res$sdrl)
  class(res) <- "drmTMB"
  return(res)
}
