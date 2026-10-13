#' Fit a drmTMB model
#'
#' @param form formula with name of response and dose scale
#' @param fform list with fixed-effect right-hand side formulas for each model parameter
#' @param rform random intercepts: a formula of the form b1 + b2 + b4 ~ 1 | cluster, or a list of such formulas, e.g. list(b1 + b4 ~ 1 | exper, b1 ~ 1 | plate). Random intercepts within a formula are correlated across the listed parameters; different formulas, and the levels of nested factors (1 | exper/plate), are independent terms.
#' @param data a data.frame object
#' @param start a vector with starting values for each model parameter. NA entries are replaced by the starting values of \code{\link{drmTMBstart}}.
#' @param fix a binary vector with same length as start. If TRUE, a model coefficient is fixed to the starting value.
#' @param weights a vector of weights for each observation.
#' @param lower lower bounds for the model coefficients b1-b5, a single value or a vector with the same length as start. Bounds of fixed coefficients are ignored; variance parameters are unbounded.
#' @param upper upper bounds for the model coefficients b1-b5, a single value or a vector with the same length as start.
#' @param family Specify the distributional assumptions. One of "gaussian", "binomial", "beta", and "poisson".
#' @param model Specify the nonlinear model. One of "logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz".
#' @param link Specify the link function for estimating asymptote parameters. One of "identity", "log", "logit". 
#' @param gsl a vector of lower bounds for a grid search, same length as start; NA for the default range. If gsl, gsu and gsn are all NULL, no grid search is performed. See \code{\link{drmTMBstart}}.
#' @param gsu a vector of upper bounds for a grid search, same length as start; NA for the default range.
#' @param gsn number of grid points for each parameter, a single number or a vector with the same length as start. \code{gsn=10} alone runs a grid search over the default ranges.
#' @param control a list of control parameters. See \code{\link{drmTMBcontrol}} for details.
#'
#' @return A drmTMB object
#' @useDynLib drmTMB

drmTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, weights=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", link="identity", gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  # control list
  control <- do.call(drmTMBcontrol, control)
  dd <- drmTMBdata(form, fform=fform, rform=rform, data=data, weights=weights, family=family, model=model, link=link)
  if (is.null(start)) start <- rep(NA, sum(dd$npars))
  if (is.null(fix)) fix <- rep(FALSE, sum(dd$npars))

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

  if (is.null(rform)){
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=obj$he, lower=lower, upper=upper)
  } else {
    opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=lower, upper=upper)
  }
  res <- list()
  res$data <- data
  res$formula <- form
  res$family <- family
  res$model <- model
  res$link <- link
  res$start <- start
  res$fix <- fix
  res$Xs <- dd$Xs
  res$fform <- dd$fform
  res$rform <- rform
  res$reterms <- dd$re$terms
  res$estimates <- opt$par
  res$opt <- opt
  res$obj <- obj
  res$sdrl <- sdreport(obj, par.fixed=opt$par)
  res$ssdrl <- summary(res$sdrl)
  class(res) <- "drmTMB"
  return(res)
}
