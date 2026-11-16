#' Marginal likelihood of a drmTMB model with priors
#'
#' The log marginal likelihood log p(y), the likelihood integrated over the prior of all parameters,
#' for Bayesian model comparison and model averaging: the posterior probability of a model is
#' proportional to its prior probability times its marginal likelihood.
#'
#' The fit needs a proper prior for every estimated parameter (see \code{\link{priors}}): all free
#' coefficients, the standard deviations of the random effects and the residual parameter (sigma, phi
#' or theta); with a flat prior, the marginal likelihood is not defined. Correlations of random effects
#' have no priors, so that random effects of several parameters or with slopes need independent terms
#' (\code{||} or separate formulas). The random effects are integrated out by the Laplace approximation
#' of the fit.
#'
#' \code{method="laplace"} approximates the log posterior of the parameters (on the scale of their
#' optimisation, with the log Jacobians of the priors) by a normal distribution at the posterior mode:
#' log p(y) = -objective + k/2 log(2 pi) - 1/2 log det H, with k parameters and the Hessian H of the
#' objective, as in the Bayesian benchmark dose modelling of EFSA (2022). The random effects are
#' integrated out first, so that this is a nested Laplace approximation. It is poor if the posterior
#' mode is close to a bound, e.g. of a pert prior or a multistage coefficient at 0, which gives a
#' warning.
#'
#' \code{method="is"} estimates the marginal likelihood by importance sampling from a multivariate t
#' distribution with \code{df} degrees of freedom, centred at the posterior mode with the covariance of
#' the fit as its scale matrix. Draws outside the bounds of the parameters have weight 0. The result
#' reports the effective sample size of the weights and the Monte Carlo standard error of the log
#' marginal likelihood; a small effective sample size (e.g. below 100) means the estimate is unreliable.
#' With random effects, each draw needs an inner optimisation, so that importance sampling is slower.
#'
#' @param object a drmTMB object, fitted with priors for all parameters
#' @param method "laplace" (default) or "is" (importance sampling)
#' @param n number of draws for importance sampling
#' @param df degrees of freedom of the multivariate t proposal of importance sampling
#'
#' @return the log marginal likelihood, a number of class "drcmarglik" with the attributes method, the
#'   number of parameters k, and for importance sampling n, the effective sample size ess, the Monte Carlo
#'   standard error se, the number of draws whose objective could not be evaluated (failed, given weight
#'   0), and the Laplace approximation laplace.
#' @export
marglik <- function(object, method=c("laplace", "is"), n=2000, df=4){
  if (!inherits(object, "drmTMB")) stop("object must be a drmTMB fit.")
  method <- match.arg(method)
  obj <- object$obj
  par <- object$opt$par
  pn <- parnames(object)
  k <- length(par)

  # proper priors for all optimised parameters
  pdat <- obj$env$data[c("prwhat", "pridx", "prtrans", "prtype", "prhyp")]
  pos <- if (length(pdat$prwhat) > 0) priorpos(pdat, object$fix, names(par)) else numeric(0)
  miss <- setdiff(seq_len(k), pos)
  if (length(miss) > 0){
    mn <- pn[miss]
    corr <- grepl("^corr_par:", mn)
    stop("the marginal likelihood needs a proper prior for every estimated parameter; ",
         if (any(!corr)) paste("no prior for ", paste(sub("^log_", "", mn[!corr]), collapse=", "), sep=""),
         if (any(!corr) && any(corr)) "; ",
         if (any(corr)) "correlations of random effects have no priors, use independent random effects (||)",
         ".", call.=FALSE)
  }
  if (!is.finite(object$opt$objective)) stop("the objective is not finite at the estimates.")
  if (!isTRUE(object$sdrl$pdHess)) stop("the Hessian at the posterior mode is not positive definite.")
  V <- object$sdrl$cov.fixed
  R <- chol(V)
  logdetV <- 2*sum(log(diag(R)))
  laplace <- -object$opt$objective + k/2*log(2*pi) + logdetV/2

  # posterior modes close to a bound: the normal approximation puts mass outside the parameter space
  lo <- if (is.null(object$lower)) rep(-Inf, k) else object$lower
  up <- if (is.null(object$upper)) rep(Inf, k) else object$upper
  se <- sqrt(diag(V))
  near <- pmin(par - lo, up - par)/se < 2
  if (method == "laplace" && any(near)){
    warning("the posterior mode of ", paste(pn[near], collapse=", "), " is within 2 standard errors of a bound; ",
            "the Laplace approximation may be poor, compare with method=\"is\".", call.=FALSE)
  }
  if (method == "laplace") return(structure(laplace, method="laplace", k=k, class="drcmarglik"))

  # importance sampling: keep the state of the TMB object, as tmbprofile does
  keep <- c("last.par.best", "random.start", "value.best", "last.par", "inner.control", "tracemgc")
  old <- mget(keep, envir=obj$env, ifnotfound=list(NULL))
  on.exit(for (v in keep) assign(v, old[[v]], envir=obj$env))
  obj$env$inner.control$trace <- FALSE
  obj$env$tracemgc <- FALSE

  z <- matrix(rnorm(n*k), n, k)
  s <- sqrt(df/rchisq(n, df))
  th <- sweep((z * s) %*% R, 2, par, "+")
  # log density of the multivariate t proposal
  maha <- rowSums(z^2) * s^2
  logq <- lgamma((df + k)/2) - lgamma(df/2) - k/2*log(df*pi) - logdetV/2 - (df + k)/2*log1p(maha/df)
  nll <- vapply(seq_len(n), function(i){
    p <- th[i, ]
    if (any(p <= lo | p >= up)) return(Inf)
    v <- tryCatch(suppressWarnings(as.numeric(obj$fn(p))), error=function(e) NA_real_)
    if (is.nan(v)) NA_real_ else v
  }, numeric(1))
  failed <- sum(is.na(nll))
  if (failed > 0) warning(failed, " of ", n, " draws could not be evaluated and have weight 0.", call.=FALSE)
  nll[is.na(nll)] <- Inf
  logw <- -nll - logq
  if (all(!is.finite(logw))) stop("all importance weights are 0.")
  mx <- max(logw)
  w <- exp(logw - mx)
  est <- mx + log(mean(w))
  ess <- sum(w)^2/sum(w^2)
  mcse <- sd(w)/(sqrt(n)*mean(w))
  structure(est, method="is", k=k, n=n, ess=ess, se=mcse, failed=failed, laplace=laplace, class="drcmarglik")
}

#' @export
print.drcmarglik <- function(x, digits=max(3, getOption("digits") - 3), ...){
  a <- attributes(x)
  cat("Log marginal likelihood:", format(as.numeric(x), digits=digits), "\n")
  if (a$method == "laplace"){
    cat("Laplace approximation, ", a$k, if (a$k == 1) " parameter\n" else " parameters\n", sep="")
  } else {
    cat("Importance sampling, ", a$n, " draws, effective sample size ", round(a$ess), "\n",
        "Monte Carlo standard error ", format(a$se, digits=2), "; Laplace approximation ", format(a$laplace, digits=digits), "\n", sep="")
    if (a$ess < 100) cat("Warning: small effective sample size; the estimate is unreliable.\n")
  }
  invisible(x)
}
