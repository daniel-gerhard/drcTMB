#' Fit a drmTMB model
#'
#' @param form formula with name of response and dose scale
#' @param fform fixed effects for the model parameters: a formula of the form b2 + b4 ~ Stress, or a list of such formulas, e.g. list(b2 ~ Stress, b4 ~ Stress + Time). Parameters not listed have a single intercept (~ 1). A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~ 1), is also accepted.
#' @param rform random effects: a formula of the form b1 + b2 + b4 ~ 1 | cluster for random intercepts, or b4 ~ x | cluster for random intercepts and slopes in a covariate x, or a list of such formulas, e.g. list(b1 + b4 ~ 1 | exper, b1 ~ 1 | plate). The random effects within a formula are correlated across the listed parameters and their intercepts and slopes; different formulas, the levels of nested factors (1 | exper/plate), and intercepts and slopes with || (b4 ~ x || cluster) are independent terms. The random slopes of a parameter usually go with a fixed effect of the covariate in fform.
#' @param data a data.frame object. Observations with missing values in the variables of the model or in the weights are left out.
#' @param start a vector with starting values for each model coefficient. NA entries are replaced by the starting values of \code{\link{drmTMBstart}}. Alternatively, a vector named by parameter, e.g. c(b1=0, b5=0), gives the same starting value to all coefficients of each named parameter.
#' @param fix a logical vector with one entry for each model coefficient. If TRUE, a model coefficient is fixed to its value in start. Alternatively, a vector named by parameter, e.g. c(b1=TRUE, b5=TRUE), fixes all coefficients of the named parameters; their values must then be given in start, e.g. start=c(b1=0, b5=0). Parameters not used by the model (b5 except for logistic, loglogistic, braincousens, cedergreen, multistage and nec, b3 for michaelismenten, asymptoticregression and exponentialdecay, b2-b5 for constant) are always fixed, at 0 unless given in start.
#' @param weights a vector of weights, one for each row of data.
#' @param lower lower bounds for the model coefficients b1-b5, a single value or a vector with the same length as start. Bounds of fixed coefficients are ignored; variance parameters are unbounded. For the multistage model, the default is 0 for the coefficients b3-b5 without covariates in fform, and for the nec model 0 for the rate b3 and, with the identity link for b4, the threshold b4 without covariates, see Details.
#' @param upper upper bounds for the model coefficients b1-b5, a single value or a vector with the same length as start.
#' @param family Specify the distributional assumptions. One of "gaussian", "binomial", "betabinomial", "beta", "poisson", "nbinom2", and "nbinom1". The negative binomial families are for overdispersed counts, with variance mu (1 + mu/theta) for nbinom2 and mu (1 + phi) for nbinom1, and the beta-binomial family for overdispersed binomial counts, with variance n mu (1 - mu) (1 + (n - 1)/(phi + 1)) for n trials, as in glmmTMB.
#' @param model Specify the nonlinear model. One of "logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz", "gamma", "multistage", "nec", "braincousens", "cedergreen", "constant". The nec model has a threshold, the braincousens and cedergreen models describe hormesis, the constant model has no dose effect, see Details.
#' @param alpha the exponent of the hormesis term of the cedergreen model, a positive number, e.g. 1, 0.5 or 0.25 as the models CRS.4a, CRS.4b and CRS.4c of drc (CRS.5a-c with b1 estimated); only for the cedergreen model. See Details.
#' @param smooth smoothing of the threshold of the nec model: the width of the smoothed threshold relative to the smallest interval between doses, a number of at least 0; 0 is an exact threshold. The default is 0.1 with random effects on b4, and 0 otherwise. Only for the nec model, see Details.
#' @param link links of the parameters b1-b5: NULL for the defaults of the model, a single link for the asymptotes b1 and b2, or links named by parameter, e.g. c(b1="logit", b2="logit", b4="log"), where the parameters not named keep their defaults. b1 and b2 have the link "identity" (default), "log" or "logit"; b4 has "identity" (default) or "log" for the models with a positive b4 (not logistic, gompertz and multistage); b5 has the log link for the asymmetry of logistic and loglogistic and the shape of nec; the other parameters have the identity link. The coefficients, and start, fix, lower, upper and linfct of \code{\link{ed}}, are on the scale of the links, see Details.
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
#' The gamma model, F(x) = pgamma(x/b4, shape=b3), is the gamma model of the benchmark dose software
#' of the US EPA (BMDS) and \code{gammadr} of drc (b3 = e, b4 = 1/b). As for michaelismenten, b1 is the
#' response at dose zero and b2 the asymptote at high doses; for quantal responses, b1 is the
#' background risk and b2 is fixed at 1, e.g. \code{start=c(b2=1), fix=c(b2=TRUE)}. The shape b3 is
#' positive; with b3 = 1, the gamma model is the asymptoticregression model (the quantal-linear model
#' of BMDS for quantal responses), so that \code{\link{anova.drmTMB}} tests the shape. BMDS restricts
#' the shape to b3 >= 1 by default, which can be set by \code{lower}.
#'
#' The multistage model of BMDS, F(x) = 1 - exp(-(b3 x + b4 x^2 + b5 x^3)), has the coefficients
#' beta1-beta3 of BMDS in b3-b5; as for the gamma model, b1 is the response at dose zero (the
#' background risk) and b2 the asymptote at high doses, fixed at 1 for quantal responses. All three
#' coefficients are estimated by default; the model of degree 2 fixes b5 at 0, and of degree 1 also
#' b4, e.g. \code{start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE)}. As in BMDS, b3-b5 are bounded
#' below by 0 by default (if \code{lower} is not given, and for parameters without covariates in
#' fform), so that the curve increases; estimates at the bound are common, and are listed by
#' \code{summary}, as their standard errors and Wald intervals do not account for the bound. Profile
#' intervals of the effective doses respect the bounds. The log link for b4 is not available.
#'
#' The threshold model nec has no effect up to a threshold b4, the no-effect concentration (NEC), and
#' F(x) = 1 - exp(-(b3 (x - b4))^exp(b5)) above it: b1 is the response up to the threshold and b2 the
#' asymptote at high doses, b3 > 0 the rate and b5 the shape on the log scale; b3 and b4 are bounded
#' below by 0 by default (if \code{lower} is not given, for parameters without covariates in fform,
#' and for b4 with the identity link). With b5 = 0, it is the NEC model of drc, NEC.4 (b1 = d,
#' b2 = c, b3 = b, b4 = e); with b2 fixed at 0, NEC.3, e.g. \code{start=c(b2=0, b5=0),
#' fix=c(b2=TRUE, b5=TRUE)}, and with b1 also fixed at 1, NEC.2. With b5 estimated, the curve is
#' sigmoidal above the threshold for b5 > 0, as the necsigm model of bayesnec (with the rate inside
#' the power), and \code{\link{anova.drmTMB}} tests the shape. Effective doses are
#' b4 + (-log(1 - p))^exp(-b5) / b3, above the threshold. The likelihood is smooth in b4 between the
#' doses, has a kink where b4 equals a dose, and often has several local maxima; the starting values
#' therefore come from thresholds within each interval between the doses. A warning is given for a
#' threshold estimated at a dose, where the standard errors are not valid (profile or bootstrap
#' intervals remain available), or at or above the highest dose, where it is not identified. With
#' b5 > 0, the curve starts flat above the threshold, and b4 and b5 are often hard to separate.
#'
#' With random effects on b4, the likelihood of each group with an exact threshold has kinks at the
#' doses, so that the Laplace approximation is unreliable and its inner optimisation has several modes. The
#' threshold is therefore smoothed by default with random effects on b4: the dose above the
#' threshold, max(x - b4, 0), is replaced by s log(1 + exp((x - b4)/s)), with a width s of
#' \code{smooth} (default 0.1) times the smallest interval between the doses. The curve then leaves
#' b1 within a few widths s around b4, and effective doses are b4 + h + s log(1 - exp(-h/s)), with h
#' the term above. In simulations with random thresholds, fits with an exact threshold mostly failed,
#' and smoothed fits converged with little bias. \code{smooth=0} gives an exact threshold also with
#' random effects, and \code{smooth} > 0 a smoothed threshold without them; the width is not
#' estimated. \code{summary} reports the width.
#'
#' The hormesis models braincousens (Brain and Cousens, 1989) and cedergreen (Cedergreen, Ritz and
#' Streibig, 2005) describe a stimulation at low doses before a decreasing curve:
#' f(x) = b1 + (b2 - b1 + b5 g(x)) / (1 + exp(-b3 (log(x) - log(b4)))), with g(x) = x for
#' braincousens and g(x) = exp(-1/x^alpha) for cedergreen. b2 is the response at dose zero (the
#' control), b1 the lower asymptote at high doses, and b5 the size of the hormesis effect on the scale
#' of the response; b5 = 0 gives the loglogistic model, so that \code{\link{anova.drmTMB}} tests for
#' hormesis. They are the models BC.5 and CRS.5a-c of drc, with b3 = -b; with b1 fixed at 0, BC.4 and
#' CRS.4a-c. The curve is decreasing at high doses for b3 < 0, and returns to b1 for the braincousens
#' model only for b3 < -1. The models are meant for decreasing curves; a warning is given otherwise.
#' The exponent alpha of the cedergreen model is not estimated, as it is hardly identifiable together
#' with b5.
#'
#' The constant model, f(x) = b1, has the same mean at all doses; b2-b5 are not used. Compared with
#' a dose-response model by \code{\link{anova.drmTMB}}, it tests for a dose effect, as
#' \code{noEffect} in drc. Under the constant model, the parameters of the curve other than an
#' asymptote are not identifiable, so the chi-square distribution of the likelihood ratio statistic
#' is only an approximation.
#'
#' With the log link for b4, \code{link=c(b4="log")}, the linear predictor of b4, including its random effects, is log(b4):
#' the coefficients of b4 in \code{coef}, \code{start}, \code{fix}, \code{lower}, \code{upper} and
#' \code{linfct} of \code{\link{ed}} are on the log scale. Random effects on b4 then vary b4 by a
#' factor between groups, and the estimate of b4 is the median over the groups. For the models on the
#' log-dose scale, this is the same model as the logistic (or gompertz) model of log(dose) for
#' loglogistic (weibull1), with b4 the log ED50, but it keeps observations at dose zero and reports
#' effective doses on the dose scale. Without random effects on b4, both links give the same fit.
#'
#' Gaussian models without random effects are estimated by maximum likelihood, but by default
#' (\code{smallsample=TRUE} in \code{\link{drmTMBcontrol}}) their inference follows lm, nls and drc:
#' the residual standard deviation is the square root of the unbiased residual variance, the weighted
#' residual sum of squares divided by n - p residual degrees of freedom, for n observations with
#' positive weights and p free coefficients b1-b5; the covariance of the coefficients
#' (\code{vcov}, \code{summary}, \code{\link{ed}}, \code{predict}) is scaled accordingly, the Wald
#' and prediction intervals use t quantiles with n - p degrees of freedom, and the profile intervals of
#' the coefficients and effective doses an F-based cutoff. The log-likelihood, the estimate log_sigma in
#' \code{coef}, AIC, \code{anova}, \code{simulate} and \code{\link{bootstrap}} remain maximum likelihood.
#' The other families and models with random effects use maximum likelihood and normal quantiles.
#'
#' @return A drmTMB object
#' @export

drmTMB <- function(form, fform=NULL, rform=NULL, data, start=NULL, fix=NULL, weights=NULL, lower=-Inf, upper=Inf, family="gaussian", model="logistic", alpha=NULL, smooth=NULL, link=NULL, gsl=NULL, gsu=NULL, gsn=NULL, control=list()){
  # control list
  control <- do.call(drmTMBcontrol, control)
  dd <- drmTMBdata(form, fform=fform, rform=rform, data=data, weights=weights, family=family, model=model, alpha=alpha, smooth=smooth, link=link)
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
  # default lower bounds of 0, e.g. for the coefficients b3-b5 of the polynomial of the multistage
  # model, as in BMDS, and for the rate b3 and threshold b4 of nec; not needed for b4 on the log scale
  m <- drcmodel(dd$mod)
  if (!is.null(m$lower0) && identical(lower, -Inf)){
    lower <- rep(-Inf, length(start))
    bpars <- if (dd$lnk[4] == 2) setdiff(m$lower0, 4) else m$lower0
    lower[dd$npi %in% which(dd$npars == 1) & dd$npi %in% bpars] <- 0
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
  res$alpha <- alpha
  res$smooth <- dd$smooth
  res$link <- dd$link
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
  # without random effects, the Hessian by automatic differentiation: the default finite differences of
  # sdreport use absolute steps, which fail for coefficients on a small scale, e.g. of a multistage model
  res$sdrl <- sdreport(obj, par.fixed=opt$par, hessian.fixed=if (is.null(rform)) obj$he(opt$par))
  # NaN standard errors of a Hessian that is not positive definite are reported by summary.drmTMB
  res$ssdrl <- suppressWarnings(summary(res$sdrl))
  # residual degrees of freedom for small-sample inference of gaussian models without random effects:
  # observations with positive weights minus free coefficients
  if (isTRUE(control$smallsample) && family == "gaussian" && length(dd$re$terms) == 0){
    rdf <- sum(dd$w > 0) - sum(!fix)
    if (rdf > 0) res$rdf <- rdf else warning("no residual degrees of freedom; the standard errors use the maximum likelihood estimate of the residual variance.", call.=FALSE)
  }
  class(res) <- "drmTMB"
  # hormesis models: a decreasing curve, which returns to b1 at high doses for braincousens (b3 < -1)
  if (!is.null(m$hormesis)){
    pv <- start
    pv[!fix] <- opt$par[seq_len(sum(!fix))]
    b3 <- as.vector(dd$Xs[[3]] %*% pv[dd$npi == 3])
    b3max <- m$hormesis$b3max
    if (any(b3 >= b3max)){
      warning("the ", model, " model is meant for decreasing curves, with b3 < ", b3max,
              "; the estimate of b3 is ", if (b3max == 0) "not negative" else paste(b3max, "or larger"), ".", call.=FALSE)
    }
  }
  # threshold models: the likelihood has a kink where the threshold equals a dose, and is flat in a
  # threshold above the highest dose
  if (m$b4 %in% "threshold" && any(!fix[dd$npi == 4])) necwarn(res)
  return(res)
}

# warnings for exact thresholds b4 of a nec fit at a dose, where the Hessian and Wald standard errors
# are not valid, and for thresholds above the highest dose, where b4 is not identified
necwarn <- function(fit){
  pv <- fit$start
  pv[!fit$fix] <- fit$opt$par[seq_len(sum(!fit$fix))]
  npi <- rep(1:5, sapply(fit$Xs, ncol))
  b4 <- unique(linkfun(linkindex(fit$link)[4], as.vector(fit$Xs[[4]] %*% pv[npi == 4]), inverse=TRUE))
  x <- fit$obj$env$data$x
  xd <- unique(x)
  tol <- 1e-6 * diff(range(x))
  # a smoothed threshold has no kink
  atdose <- if (smoothwidth(fit$smooth) == 0) b4[sapply(b4, function(b) any(abs(b - xd[xd > 0]) < tol))]
  if (length(atdose) > 0){
    warning("the threshold b4 is at the dose ", paste(signif(atdose, 4), collapse=", "), ", where the likelihood has a kink; ",
            "the standard errors are not reliable, use profile or bootstrap intervals.", call.=FALSE)
  }
  if (any(b4 >= max(x) - tol)){
    warning("the threshold b4 is at or above the highest dose; it is not identified by the data.", call.=FALSE)
  }
}
