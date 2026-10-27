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
#' For gaussian models without random effects, \code{test="F"} gives F-tests, as drc does for
#' continuous responses: the reduction of the (weighted) residual sum of squares per additional
#' parameter, divided by the residual mean square of the largest model, with n - p residual degrees
#' of freedom for p coefficients b1-b5. The F-test is exact for linear models and accounts for the
#' estimation of the residual variance, which the likelihood ratio test does not, so it is usually
#' preferable for small experiments.
#'
#' @param object a drmTMB object
#' @param ... further drmTMB objects
#' @param test "Chisq" (default) for likelihood ratio tests, or "F" for F-tests of gaussian models
#'   without random effects
#'
#' @return a data frame of class anova with the number of parameters (npar), AIC, BIC, log-likelihood
#'   and deviance (minus twice the log-likelihood) of each model, and the likelihood ratio statistic
#'   (Chisq), its degrees of freedom (Df) and p-value (Pr(>Chisq)) for the comparison with the
#'   previous model. For \code{test="F"}, the residual degrees of freedom (Res.Df) and residual sum of
#'   squares (RSS) instead of the deviance, and the F statistic (F) and its p-value (Pr(>F)).
#' @export
anova.drmTMB <- function(object, ..., test=c("Chisq", "F")){
  test <- match.arg(test)
  fits <- c(list(object), list(...))
  if (!all(vapply(fits, inherits, logical(1), "drmTMB"))) stop("all arguments must be drmTMB objects.")
  if (length(fits) < 2) stop("anova needs at least two models fitted to the same data, e.g. anova(fit0, fit1).")
  mc <- as.list(match.call())[-1]
  nm <- vapply(mc[names(mc) != "test"], deparse1, character(1))

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
  ic <- data.frame(npar=np, AIC=vapply(fits[ord], AIC, numeric(1)), BIC=vapply(fits[ord], BIC, numeric(1)), logLik=ll)

  if (test == "F"){
    if (!all(vapply(fits, function(f) f$family == "gaussian" && length(f$reterms) == 0, logical(1)))){
      stop("F-tests are only available for gaussian models without random effects; use test=\"Chisq\".")
    }
    # weighted residual sums of squares, residual degrees of freedom without the residual sd
    rss <- vapply(fits[ord], function(f) sum(f$obj$env$data$w * residuals(f, type="response")^2), numeric(1))
    rdf <- nobs(object) - (np - 1)
    k <- length(fits)
    fval <- c(NA, -diff(rss)/df[-1] / (rss[k]/rdf[k]))
    tab <- data.frame(ic, Res.Df=rdf, RSS=rss, F=fval, Df=df,
                      "Pr(>F)"=ifelse(!is.na(df) & df > 0, pf(pmax(fval, 0), df, rdf[k], lower.tail=FALSE), NA),
                      row.names=nm[ord], check.names=FALSE)
    return(structure(tab, heading="F-tests of drmTMB models\n", class=c("anova", "data.frame")))
  }

  pval <- ifelse(!is.na(df) & df > 0, pchisq(pmax(chisq, 0), df, lower.tail=FALSE), NA)
  tab <- data.frame(ic, deviance=-2*ll, Chisq=chisq, Df=df, "Pr(>Chisq)"=pval,
                    row.names=nm[ord], check.names=FALSE)
  structure(tab, heading="Likelihood ratio tests of drmTMB models\n", class=c("anova", "data.frame"))
}
