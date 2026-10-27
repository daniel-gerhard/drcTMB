#' Likelihood ratio tests for drmTMB models
#'
#' Compares two or more models fitted to the same data by likelihood ratio tests, e.g. a model with
#' and without covariates in \code{fform}, with and without random effects, or a poisson and a
#' negative binomial model.
#'
#' The models are ordered by their number of parameters, and each model is compared with the
#' previous one: the likelihood ratio statistic is twice the difference of their log-likelihoods
#' (the Laplace approximation of the marginal log-likelihoods for models with random effects), with
#' the difference in the number of parameters as degrees of freedom. The test assumes that the
#' smaller model is nested in the larger one, which is not checked. For parameters on the boundary
#' of their range under the smaller model, e.g. a standard deviation of random effects of zero, or an
#' infinite size of the negative binomial for the poisson model, the chi-square distribution is
#' conservative: the p-value is approximately twice the correct value for a single standard
#' deviation.
#'
#' @param object a drmTMB object
#' @param ... further drmTMB objects
#'
#' @return a data frame of class anova with the number of parameters (npar), AIC, BIC, log-likelihood
#'   and deviance (minus twice the log-likelihood) of each model, and the likelihood ratio statistic
#'   (Chisq), its degrees of freedom (Df) and p-value (Pr(>Chisq)) for the comparison with the
#'   previous model
#' @export
anova.drmTMB <- function(object, ...){
  fits <- c(list(object), list(...))
  if (!all(vapply(fits, inherits, logical(1), "drmTMB"))) stop("all arguments must be drmTMB objects.")
  if (length(fits) < 2) stop("anova needs at least two models fitted to the same data, e.g. anova(fit0, fit1).")
  nm <- vapply(as.list(match.call())[-1], deparse1, character(1))

  # the same observations, trials and weights
  dat <- lapply(fits, function(f) f$obj$env$data[c("y", "bn", "w")])
  if (!all(vapply(dat[-1], identical, logical(1), dat[[1]]))) stop("the models are not all fitted to the same data.")

  ll <- lapply(fits, logLik)
  np <- vapply(ll, function(l) attr(l, "df"), numeric(1))
  ord <- order(np)
  ll <- vapply(ll, as.numeric, numeric(1))[ord]
  np <- np[ord]
  chisq <- c(NA, 2*diff(ll))
  df <- c(NA, diff(np))
  if (any(chisq < -1e-6, na.rm=TRUE)) warning("a model with more parameters has a lower log-likelihood; check the convergence of the fits.", call.=FALSE)
  pval <- ifelse(!is.na(df) & df > 0, pchisq(pmax(chisq, 0), df, lower.tail=FALSE), NA)

  tab <- data.frame(npar=np, AIC=vapply(fits[ord], AIC, numeric(1)), BIC=vapply(fits[ord], BIC, numeric(1)),
                    logLik=ll, deviance=-2*ll, Chisq=chisq, Df=df, "Pr(>Chisq)"=pval,
                    row.names=nm[ord], check.names=FALSE)
  structure(tab, heading="Likelihood ratio tests of drmTMB models\n", class=c("anova", "data.frame"))
}
