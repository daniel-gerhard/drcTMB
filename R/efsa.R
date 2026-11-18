#' Default priors of the EFSA Bayesian benchmark dose modelling
#'
#' The default (weakly informative) priors of the Bayesian benchmark dose modelling of EFSA (2022), as in
#' its R package BMABMDR (Kremer et al. 2026), for a model in the BMD parameterisation of
#' \code{\link{drmTMB}} (\code{bmdpar}): PERT priors for the response at dose zero (background), the
#' response at high doses (maximum response) and the BMD, lognormal priors for the steepness and for the
#' residual standard deviation. The result gives the arguments of \code{drmTMB}; \code{\link{bmdefsa}}
#' fits a set of models with these priors and averages them.
#'
#' The models are those of EFSA (2022) available in drcTMB: weibull2 (exponential), weibull1 (inverse
#' exponential), loglogistic (Hill), lognormal (log-normal), gamma, and multistage of degree 2 (quadratic
#' exponential, LMS two-stage), and for quantal responses also probit and logit (the probit and logistic
#' models of EFSA, pnorm(b1 + b4 x^b3) and plogis(b1 + b4 x^b3): the power model with
#' \code{curvelink="probit"} or \code{"logit"}). Families: "gaussian" (normal distribution), "lognormal"
#' (fitted with \code{curvelink="log"}: the curve is the log of the median, as the lognormal models of
#' EFSA) and "binomial" (quantal responses, with b2 fixed at 1).
#'
#' Default priors, with the response means at dose zero (y0) and at the highest dose (yN) (geometric means
#' for lognormal), the benchmark response q and the highest dose dN:
#' \itemize{
#'   \item background, increasing curve: pert(0.001, y0, 2 y0); decreasing: pert(0.5 y0, y0, 2 y0);
#'     quantal: pert(max(L/2, 1/(10 n0)), max(p0, 1/(5 n0)), min(3 U/2, 1 - 1/(10 n0))) with the observed
#'     proportion p0 of the n0 controls and its confidence interval (L, U) of prop.test (lower bound 1e-8
#'     without events in the controls); shape 4. For the
#'     probit and logit models the prior is on the background probability pnorm(b1) or plogis(b1).
#'   \item maximum response (continuous): increasing pert(y0 (1.01 + q), yN, 2 yN), or pert(y0 (1.01 + q),
#'     3 yN, 6 yN) if the curve does not level off at the high doses; decreasing pert(0.5 yN, yN,
#'     y0 (1 - q - 0.01)), or pert(0.1 yN, 0.5 yN, y0 (1 - q - 0.01)); shape 4. Whether the curve levels off
#'     follows from a fractional polynomial of degree 2 of the means: its last slope on the log-dose scale is
#'     below half its largest slope. EFSA puts this prior on the fold change, the maximum response divided by
#'     the mode of the background; here it is on the maximum response (the asymptote) itself.
#'   \item BMD: pert(0, dN/2, 3 dN) (\code{extended=TRUE}) or pert(0, dN/2, dN), with shape 0.0001, nearly
#'     uniform.
#'   \item steepness d: |b3| lognormal(1, 1) (EFSA truncates log d at 5, which removes a probability of
#'     3e-5 and is left out) or, with \code{dprior="alternative"}, lognormal(0.4, sqrt(0.5)), for a
#'     sensitivity analysis (the power b3 of the probit and logit models has the same prior); for
#'     decreasing curves of the log-dose models b3 is negative, with the prior
#'     \code{negate(lognormal(...))}. The multistage model of degree 2 in the BMD parameterisation is the
#'     quadratic exponential model of EFSA, with the technical parameter log_d (in place of b3), with the
#'     prior normal(0, 1) of EFSA.
#'   \item residual standard deviation: lognormal(log(1.5 s), 0.5), with the mean s of the standard
#'     deviations within the doses (of log(y) for lognormal), i.e. log(1/sigma^2) normal with mean
#'     -2 log(1.5 s) and standard deviation 1.
#' }
#' For the lognormal family, the priors of the background and the maximum response are priors of the
#' median, exp(b1) or exp(b2). Informative priors are given as \code{c(min, mode, max)} (mode NA for the
#' middle of the range) or as a prior, e.g. from \code{\link{fitpert}}.
#'
#' @param form formula of the response and the dose, as in \code{\link{drmTMB}}
#' @param data a data.frame
#' @param model one of "weibull2", "weibull1", "loglogistic", "lognormal", "gamma", "multistage", and for
#'   the binomial family "probit" and "logit"
#' @param family "gaussian", "lognormal" or "binomial"
#' @param bmr,def the benchmark response and its definition (see \code{\link{bmd}}); by default
#'   "relative" for gaussian and lognormal, "excess" for binomial
#' @param background,maxresponse,bmd informative priors: \code{c(min, mode, max)} or a prior (see
#'   \code{\link{priors}}); NULL for the defaults
#' @param dprior "default" or "alternative" prior of the steepness
#' @param extended logical, the upper bound of the default BMD prior is three times the highest dose (TRUE)
#'   or the highest dose
#' @param x an efsapriors object
#' @param ... not used
#'
#' @references EFSA Scientific Committee (2022). Guidance on the use of the benchmark dose approach in
#'   risk assessment. EFSA Journal 20(10), 7584.
#'
#'   Kremer, C., Saenen, N. D., Bijnens, K., Faes, C., De Boever, P., Smeets, K. and Aerts, M. (2026).
#'   Compilation of NTP and published dose-response data: study of informative priors in the EFSA platform
#'   for Bayesian benchmark dose analysis. EFSA Supporting Publications, EN-9902.
#'
#' @return an efsapriors object, a list with the EFSA model (name), the arguments model, family,
#'   curvelink, start, fix, prior and bmdpar of \code{drmTMB}, the direction of the curve (increasing) and whether it levels off (flat)
#' @export
#'
#' @examples
#' # default priors of EFSA (2022) for the log-logistic model of the ryegrass data
#' efsapriors(rootl ~ conc, data = ryegrass, model = "loglogistic", family = "gaussian")
#' # an informative prior of the BMD, by its minimum, mode and maximum
#' efsapriors(rootl ~ conc, data = ryegrass, model = "loglogistic", family = "gaussian",
#'            bmd = c(0.1, 1, 5))
efsapriors <- function(form, data, model, family=c("gaussian", "lognormal", "binomial"), bmr=NULL, def=NULL,
                       background=NULL, maxresponse=NULL, bmd=NULL, dprior=c("default", "alternative"), extended=TRUE){
  family <- match.arg(family)
  dprior <- match.arg(dprior)
  models <- c("weibull2", "weibull1", "loglogistic", "lognormal", "gamma", "multistage", "probit", "logit")
  if (length(model) != 1 || !(model %in% models)) stop("model must be one of ", paste(models, collapse=", "), ".")
  quantal <- family == "binomial"
  # probit and logit: the power model with a probit or logit curve link, for quantal responses only
  qlink <- model %in% c("probit", "logit")
  if (qlink && !quantal) stop("the ", model, " model is for quantal responses (binomial family).")
  if (is.null(def)) def <- if (quantal) "excess" else "relative"
  if (is.null(bmr)) bmr <- if (quantal) 0.1 else 0.05

  # response means at each dose: means (gaussian), geometric means (lognormal) or proportions (binomial)
  mf <- model.frame(form, data=data)
  yr <- model.response(mf, "any")
  x <- mf[, 2]
  xd <- sort(unique(x))
  if (length(xd) < 3 || xd[1] != 0) stop("the data need a control group at dose zero and at least two other doses.")
  dN <- max(xd)
  if (quantal){
    if (ncol(as.matrix(yr)) != 2) stop("give the quantal response as cbind(affected, unaffected).")
    yy <- tapply(yr[, 1], x, sum)
    nn <- tapply(yr[, 1] + yr[, 2], x, sum)
    md <- as.vector(yy/nn)
  } else {
    y <- as.vector(yr)
    if (any(y <= 0)) stop("the EFSA priors need a positive response.")
    z <- if (family == "lognormal") log(y) else y
    md <- as.vector(tapply(z, x, mean))
    sds <- as.vector(tapply(z, x, sd))
    if (family == "lognormal") md <- exp(md)
    if (all(is.na(sds))) stop("the prior of the residual standard deviation needs replicates at the doses.")
  }
  y0 <- md[1]
  yN <- md[length(md)]
  inc <- yN > y0
  if (quantal && !inc) stop("quantal responses must increase with the dose (adverse events).")
  flat <- if (quantal) NA else fpflat(xd/dN, md, inc)

  # parameters of the background (response at dose zero) and of the maximum response: the asymptotes, for
  # the log-dose models with b3 < 0 for decreasing curves (b2 at dose zero)
  logdose <- model %in% c("weibull2", "weibull1", "loglogistic", "lognormal")
  pbkg <- if (logdose && !inc) "b2" else "b1"
  pmax <- if (pbkg == "b1") "b2" else "b1"

  # background
  if (quantal){
    n0 <- nn[[1]]
    ci <- stats::prop.test(yy[[1]], n0)$conf.int
    bk <- c(if (yy[[1]] > 0) max(ci[1]/2, 1/(10*n0)) else 1e-8, max(y0, 1/(5*n0)), min(3*ci[2]/2, 1 - 1/(10*n0)))
  } else {
    bk <- if (inc) c(0.001, y0, 2*y0) else c(0.5*y0, y0, 2*y0)
  }
  pr <- list()
  pr[[pbkg]] <- efsapert(background, bk, 4, "background")
  if (qlink) pr$b1 <- oninv(pr$b1, if (model == "probit") "pnorm" else "plogis")
  # maximum response
  if (!quantal){
    mx <- if (inc){
      if (flat) c(y0*(1.01 + bmr), yN, 2*yN) else c(y0*(1.01 + bmr), 3*yN, 6*yN)
    } else {
      if (flat) c(0.5*yN, yN, y0*(1 - bmr - 0.01)) else c(0.1*yN, 0.5*yN, y0*(1 - bmr - 0.01))
    }
    pr[[pmax]] <- efsapert(maxresponse, mx, 4, "maximum response")
  } else if (!is.null(maxresponse)){
    stop("quantal models have no maximum response parameter (b2 is fixed at 1).")
  }
  if (family == "lognormal"){
    pr[[pbkg]] <- onexp(pr[[pbkg]])
    pr[[pmax]] <- onexp(pr[[pmax]])
  }
  # BMD
  pr$bmd <- efsapert(bmd, c(0, dN/2, if (extended) 3*dN else dN), 0.0001, "BMD")
  # steepness
  dp <- if (dprior == "default") lognormal(1, 1) else lognormal(0.4, sqrt(0.5))
  # quadratic exponential: b3 is log d of the BMD parameterisation (see drmTMB), with EFSA's normal(0, 1)
  pr$b3 <- if (model == "multistage") normal(0, 1) else if (logdose && !inc) negate(dp) else dp
  # residual standard deviation
  if (!quantal) pr$sigma <- lognormal(log(1.5*mean(sds, na.rm=TRUE)), 0.5)

  start <- c(b5=0)
  fix <- c(b5=TRUE)
  if (quantal){
    start <- c(start, b2=1)
    fix <- c(fix, b2=TRUE)
  }
  res <- list(name=model, model=if (qlink) "power" else model, family=family,
              curvelink=if (qlink) model else if (family == "lognormal") "log" else "identity", start=start, fix=fix,
              prior=pr[c(intersect(c("b1", "b2"), names(pr)), "b3", "bmd", if (!quantal) "sigma")],
              bmdpar=list(bmr=bmr, def=def), increasing=inc, flat=flat)
  class(res) <- "efsapriors"
  res
}

#' @rdname efsapriors
#' @export
print.efsapriors <- function(x, ...){
  cat("EFSA priors for the ", x$name, " model, ", x$family, " family", if (x$curvelink != "identity") paste(" (", x$curvelink, " curve link)", sep=""),
      ", ", if (x$increasing) "increasing" else "decreasing", " curve", if (isTRUE(x$flat)) ", levelling off" else if (isFALSE(x$flat)) ", not levelling off",
      "\nBMD: ", x$bmdpar$def, ", bmr = ", format(x$bmdpar$bmr), "\n\n", sep="")
  for (nm in names(x$prior)) cat("  ", nm, " ~ ", x$prior[[nm]]$name, "\n", sep="")
  invisible(x)
}

# a PERT prior from user input p (c(min, mode, max) with mode NA for the middle, or a prior), or the
# default c(min, mode, max) d with shape; a mode outside the range is moved inside, as in BMABMDR
efsapert <- function(p, d, shape, what){
  if (inherits(p, "drcprior")) return(p)
  if (!is.null(p)){
    if (!is.numeric(p) || length(p) != 3) stop("the prior of the ", what, " must be c(min, mode, max) or a prior.")
    if (is.na(p[2])) p[2] <- (p[1] + p[3])/2
    d <- p
  }
  if (!(d[1] < d[3])) stop("the prior of the ", what, " needs min < max: ", paste(signif(d, 4), collapse=", "), ".")
  if (d[2] >= d[3]) d[2] <- d[3] - 0.05*(d[3] - d[1])
  if (d[2] <= d[1]) d[2] <- d[1] + 0.05*(d[3] - d[1])
  pert(d[1], d[2], d[3], shape=shape)
}

# whether the dose-response curve of the means m at the doses x (divided by the highest dose) levels off,
# as flat() of BMABMDR: a fractional polynomial of degree 2, whose last slope on the log-dose scale is
# below half its largest slope (absolute slopes for decreasing curves)
fpflat <- function(x, m, inc){
  xx <- x + 0.0001
  pw <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  fpterm <- function(p) if (p == 0) log(xx) else xx^p
  best <- NULL
  for (i in seq_along(pw)) for (j in i:length(pw)){
    t1 <- fpterm(pw[i])
    t2 <- if (i == j) t1*log(xx) else fpterm(pw[j])
    f <- tryCatch(lm.fit(cbind(1, t1, t2), m), error=function(e) NULL)
    if (is.null(f) || anyNA(f$coefficients)) next
    rss <- sum(f$residuals^2)
    if (is.null(best) || rss < best$rss) best <- list(rss=rss, fit=m - f$residuals)
  }
  if (is.null(best)) return(FALSE)
  sl <- diff(best$fit)/diff(log(xx))
  if (!inc) sl <- abs(sl)
  sl[length(sl)]/max(sl) < 0.5
}

#' Bayesian model-averaged BMD with the EFSA defaults
#'
#' Fits the models of the Bayesian benchmark dose modelling of EFSA (2022) available in drcTMB, in the BMD
#' parameterisation with the priors of \code{\link{efsapriors}}, for each model and family, and averages
#' them by \code{\link{bmdbayes}}.
#'
#' @param form,data formula and data, as in \code{\link{drmTMB}}
#' @param bmr,def the benchmark response and its definition (see \code{\link{efsapriors}})
#' @param models the models (see \code{\link{efsapriors}}), by default all: six for continuous responses,
#'   and eight, including probit and logit, for quantal responses
#' @param family the families: for continuous responses "gaussian" and "lognormal" (both by default),
#'   for quantal responses "binomial"
#' @param ... further arguments of \code{\link{efsapriors}}: background, maxresponse, bmd, dprior, extended
#' @param prior.weights,method,n,df see \code{\link{bmdbayes}}; prior weights in the order of the models
#'   (for each model the families). The default method is importance sampling ("is"): with the PERT priors,
#'   the posterior modes are often close to a bound of a prior, e.g. a background close to 0, where the
#'   Laplace approximation is poor.
#' @param control a list of control parameters of \code{\link{drmTMB}}
#'
#' @return a bmdbayes object, with the fits (\code{fits}) and the priors (\code{efsa}) of the models; models
#'   whose fit fails are left out with a warning
#' @export
#'
#' @examples
#' \donttest{
#' # quantal: liver tumours after exposure to TCDD
#' set.seed(1)
#' (ba <- bmdefsa(cbind(incidence, total - incidence) ~ conc, data = TCDD, family = "binomial"))
#' }
bmdefsa <- function(form, data, bmr=NULL, def=NULL, models=NULL,
                    family=c("gaussian", "lognormal"), ..., prior.weights=NULL, method=c("is", "laplace"), n=5000, df=4, control=list()){
  method <- match.arg(method)
  if (!all(family %in% c("gaussian", "lognormal", "binomial"))) stop("family must be gaussian and/or lognormal, or binomial.")
  if ("binomial" %in% family && length(family) > 1) stop("quantal (binomial) and continuous families cannot be averaged together.")
  if (is.null(models)) models <- c("weibull2", "weibull1", "loglogistic", "lognormal", "gamma", "multistage",
                                   if (identical(family, "binomial")) c("probit", "logit"))
  grid <- expand.grid(family=family, model=models, stringsAsFactors=FALSE)
  lab <- paste(grid$model, c(gaussian="N", lognormal="LN", binomial="Q")[grid$family], sep="_")
  if (!is.null(prior.weights) && length(prior.weights) != nrow(grid)) stop("prior.weights must have ", nrow(grid), " entries, one for each model and family.")
  efs <- lapply(seq_len(nrow(grid)), function(k) efsapriors(form, data, grid$model[k], grid$family[k], bmr=bmr, def=def, ...))
  fits <- lapply(seq_len(nrow(grid)), function(k){
    e <- efs[[k]]
    tryCatch(drmTMB(form, data=data, model=e$model, family=e$family, curvelink=e$curvelink, start=e$start, fix=e$fix,
                    prior=e$prior, bmdpar=e$bmdpar, control=control),
             error=function(err){
               warning("model ", lab[k], " is left out: ", conditionMessage(err), call.=FALSE)
               NULL
             })
  })
  names(fits) <- names(efs) <- lab
  ok <- !vapply(fits, is.null, logical(1))
  if (!any(ok)) stop("none of the models could be fitted.")
  args <- c(fits[ok], list(prior.weights=if (!is.null(prior.weights)) prior.weights[ok], method=method, n=n, df=df))
  res <- do.call(bmdbayes, args)
  res$models$model <- grid$model[ok]
  res$fits <- fits[ok]
  res$efsa <- efs
  res
}

#' A PERT prior from a posterior distribution
#'
#' Fits a PERT distribution to the posterior distribution of a parameter, e.g. the model-averaged BMD of
#' \code{\link{bmdbayes}}, as an informative prior for the analysis of a new study (Kremer et al. 2026,
#' EFSA EN-9902): the posterior is summarised by 1000 quantiles; the minimum and maximum of the PERT are
#' their smallest and largest value, the mode is the mode of their kernel density, and the shape is
#' estimated by maximum likelihood (between 0 and 20). Several posteriors, e.g. of several historical
#' studies, are mixed with equal weights, the recommended way to combine studies. With discounting, the
#' shape is reduced so that the variance of the PERT doubles.
#'
#' @param ... bmdbayes objects (the posterior of the BMD), or numeric vectors of posterior draws
#' @param weights for a single numeric vector of draws, their weights (e.g. of importance sampling)
#' @param discount logical, discount the prior (double its variance)
#'
#' @return a prior (see \code{\link{priors}}), e.g. for \code{prior=list(bmd=...)} of \code{drmTMB} or the
#'   bmd argument of \code{\link{efsapriors}}
#' @export
#'
#' @examples
#' # a PERT prior from posterior draws of a BMD
#' set.seed(1)
#' draws <- exp(rnorm(2000, log(2), 0.3))
#' fitpert(draws)
#'
#' # from the model-averaged posterior of the BMD of bmdbayes (or bmdefsa), e.g. of a historical
#' # study, as an informative prior of the BMD for a current study
#' pr <- list(bmd = pert(0.1, 1, 5), b2 = normal(8, 1), b3 = normal(-2, 2), sigma = half_normal(1))
#' fit <- function(model)
#'   drmTMB(rootl ~ conc, data = ryegrass, model = model,
#'          start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE),
#'          bmdpar = list(bmr = 0.05, def = "relative"), prior = pr)
#' ba <- bmdbayes(LL = fit("loglogistic"), W2 = fit("weibull2"))
#' fitpert(ba)
fitpert <- function(..., weights=NULL, discount=FALSE){
  post <- list(...)
  if (length(post) == 0) stop("give at least one posterior.")
  probs <- (seq_len(1000) - 0.5)/1000
  qs <- lapply(post, function(p){
    if (inherits(p, "bmdbayes")) return(exp(bmdbayesq(p, probs)))
    if (!is.numeric(p)) stop("the posteriors must be bmdbayes objects or numeric draws.")
    if (!is.null(weights)){
      if (length(post) > 1 || length(weights) != length(p)) stop("weights are for a single vector of draws.")
      return(wquantile(p, weights, probs))
    }
    unname(quantile(p, probs))
  })
  # mixture with equal weights: the quantiles of the pooled quantiles
  y <- unlist(qs)
  if (length(qs) > 1) y <- unname(quantile(y, probs))
  lo <- min(y) - 1e-9
  hi <- max(y) + 1e-9
  dx <- stats::density(y)
  mo <- dx$x[which.max(dx$y)]
  mo <- min(max(mo, lo + 1e-6*(hi - lo)), hi - 1e-6*(hi - lo))
  nll <- function(s){
    a <- 1 + s*(mo - lo)/(hi - lo)
    b <- 1 + s*(hi - mo)/(hi - lo)
    -sum(dbeta((y - lo)/(hi - lo), a, b, log=TRUE))
  }
  s <- stats::optimize(nll, c(1e-9, 20))$minimum
  if (discount){
    mu <- (lo + s*mo + hi)/(s + 2)
    v <- 2*(mu - lo)*(hi - mu)/(s + 3)
    s <- max((mu - lo)*(hi - mu)/v - 3, 0.0001)
  }
  pert(lo, mo, hi, shape=s)
}
