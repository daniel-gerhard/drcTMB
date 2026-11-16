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
#' \code{method="is"} estimates the marginal likelihood by adaptive importance sampling. The proposal is
#' a multivariate t distribution with \code{df} degrees of freedom, first centred at the posterior mode
#' with the covariance of the fit as its scale matrix; three rounds of n/2 draws move it to the weighted
#' mean and covariance of the draws, which follows skewed or curved posteriors (e.g. a strong
#' correlation of the log BMD with the slope), and a final round of n draws gives the estimate. The
#' proposal of each round is a mixture, with probability 0.2 of the initial proposal, so that the
#' weights stay bounded in its tails. Draws outside the bounds of the parameters have weight 0. The
#' result reports the effective sample size of the weights of the final round and the Monte Carlo
#' standard error of the log marginal likelihood; a small effective sample size (e.g. below 100) means
#' the estimate is unreliable. With random effects, each draw needs an inner optimisation, so that
#' importance sampling is slower.
#'
#' @param object a drmTMB object, fitted with priors for all parameters
#' @param method "laplace" (default) or "is" (importance sampling)
#' @param n number of draws of the final round of importance sampling
#' @param df degrees of freedom of the multivariate t proposal of importance sampling
#'
#' @return the log marginal likelihood, a number of class "drcmarglik" with the attributes method, the
#'   number of parameters k, and for importance sampling n, the effective sample size ess, the Monte Carlo
#'   standard error se, the number of draws whose objective could not be evaluated (failed, given weight
#'   0), and the Laplace approximation laplace.
#' @export
marglik <- function(object, method=c("laplace", "is"), n=5000, df=4){
  if (!inherits(object, "drmTMB")) stop("object must be a drmTMB fit.")
  method <- match.arg(method)
  mlcheck(object)
  la <- mllaplace(object, warn=method == "laplace")
  if (method == "laplace") return(structure(la, method="laplace", k=length(object$opt$par), class="drcmarglik"))
  is <- mlis(object, n, df)
  structure(is$logml, method="is", k=length(object$opt$par), n=n, ess=is$ess, se=is$se, failed=is$failed, laplace=la,
            class="drcmarglik")
}

# stop unless every optimised parameter of a fit has a proper prior, and its Hessian is positive definite
mlcheck <- function(object){
  par <- object$opt$par
  pn <- parnames(object)
  pdat <- object$obj$env$data[c("prwhat", "pridx", "prtrans", "prtype", "prpost", "prhyp")]
  pos <- if (length(pdat$prwhat) > 0) priorpos(pdat, object$fix, names(par)) else numeric(0)
  miss <- setdiff(seq_along(par), pos)
  if (length(miss) > 0){
    mn <- pn[miss]
    corr <- grepl("^corr_par:", mn)
    stop("the marginal likelihood needs a proper prior for every estimated parameter; ",
         if (any(!corr)) paste("no prior for ", paste(sub("^log_", "", mn[!corr]), collapse=", "), sep=""),
         if (any(!corr) && any(corr)) "; ",
         if (any(corr)) "correlations of random effects have no priors, use independent random effects (||)",
         ".", call.=FALSE)
  }
  invisible(TRUE)
}

# the bounds of the optimised parameters of a fit
mlbounds <- function(object){
  k <- length(object$opt$par)
  list(lo=if (is.null(object$lower)) rep(-Inf, k) else object$lower,
       up=if (is.null(object$upper)) rep(Inf, k) else object$upper)
}

# Laplace approximation of the log marginal likelihood at the posterior mode; with warn, a warning for
# posterior modes within 2 standard errors of a bound, where the normal approximation puts mass outside
# the parameter space
mllaplace <- function(object, warn=TRUE){
  if (!is.finite(object$opt$objective)) stop("the objective is not finite at the estimates.")
  if (!isTRUE(object$sdrl$pdHess)) stop("the Hessian at the posterior mode is not positive definite.")
  par <- object$opt$par
  k <- length(par)
  V <- object$sdrl$cov.fixed
  logdetV <- 2*sum(log(diag(chol(V))))
  b <- mlbounds(object)
  near <- pmin(par - b$lo, b$up - par)/sqrt(diag(V)) < 2
  if (warn && any(near)){
    warning("the posterior mode of ", paste(parnames(object)[near], collapse=", "), " is within 2 standard errors of a bound; ",
            "the Laplace approximation may be poor, compare with method=\"is\".", call.=FALSE)
  }
  -object$opt$objective + k/2*log(2*pi) + logdetV/2
}

# log density of the multivariate t distribution with df degrees of freedom, location m and the upper
# Cholesky factor R of its scale matrix, at the rows of x
ldmvt <- function(x, m, R, df){
  k <- length(m)
  z <- backsolve(R, t(x) - m, transpose=TRUE)
  lgamma((df + k)/2) - lgamma(df/2) - k/2*log(df*pi) - sum(log(diag(R))) - (df + k)/2*log1p(colSums(z^2)/df)
}

# n draws (rows) from the multivariate t distribution
rmvt <- function(n, m, R, df){
  k <- length(m)
  z <- matrix(rnorm(n*k), n, k)
  sweep((z * sqrt(df/rchisq(n, df))) %*% R, 2, m, "+")
}

# adaptive importance sampling of the log marginal likelihood. The proposal is a defensive mixture of
# multivariate t distributions with df degrees of freedom: with probability alpha centred at the
# posterior mode with the covariance of the fit as its scale matrix, otherwise at the weighted mean with
# the weighted covariance of the draws of the previous round (the first round uses the mode). iters
# rounds of n/2 draws adapt the proposal, a final round of n draws gives the estimate, its effective
# sample size and Monte Carlo standard error, the number of failed draws, and the draws (rows of the
# optimised parameters) with their normalised weights.
mlis <- function(object, n, df, iters=3, alpha=0.2){
  obj <- object$obj
  par <- object$opt$par
  k <- length(par)
  R0 <- chol(object$sdrl$cov.fixed)
  b <- mlbounds(object)
  # keep the state of the TMB object, as tmbprofile does
  keep <- c("last.par.best", "random.start", "value.best", "last.par", "inner.control", "tracemgc")
  old <- mget(keep, envir=obj$env, ifnotfound=list(NULL))
  on.exit(for (v in keep) assign(v, old[[v]], envir=obj$env))
  obj$env$inner.control$trace <- FALSE
  obj$env$tracemgc <- FALSE
  nllfun <- function(p){
    if (any(p <= b$lo | p >= b$up)) return(Inf)
    v <- tryCatch(suppressWarnings(as.numeric(obj$fn(p))), error=function(e) NA_real_)
    if (is.nan(v)) NA_real_ else v
  }

  m <- par
  R <- R0
  for (it in 0:iters){
    nd <- if (it < iters) max(ceiling(n/2), 10*k) else n
    n0 <- round(alpha*nd)
    th <- rbind(rmvt(n0, par, R0, df), rmvt(nd - n0, m, R, df))
    l0 <- ldmvt(th, par, R0, df)
    l1 <- ldmvt(th, m, R, df)
    mx <- pmax(l0, l1)
    logq <- mx + log(alpha*exp(l0 - mx) + (1 - alpha)*exp(l1 - mx))
    nll <- apply(th, 1, nllfun)
    logw <- -nll - logq
    logw[is.na(logw)] <- -Inf
    if (all(!is.finite(logw))) stop("all importance weights are 0.")
    w <- exp(logw - max(logw))
    if (it < iters){
      # adapt to the weighted draws, unless too few draws carry the weight
      wn <- w/sum(w)
      if (sum(w)^2/sum(w^2) >= 2*k){
        S <- cov.wt(th, wn)$cov
        Rn <- tryCatch(chol(S), error=function(e) NULL)
        if (!is.null(Rn)){
          m <- colSums(th*wn)
          R <- Rn
        }
      }
    }
  }
  failed <- sum(is.na(nll))
  if (failed > 0) warning(failed, " of ", n, " draws could not be evaluated and have weight 0.", call.=FALSE)
  colnames(th) <- parnames(object)
  list(logml=max(logw) + log(mean(w)), ess=sum(w)^2/sum(w^2), se=sd(w)/(sqrt(n)*mean(w)), failed=failed,
       draws=th, weights=w/sum(w))
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
