#' drmTMB Methods
#'
#' @param x,object a drmTMB object
#' @param newdata a data frame for prediction. If NULL, the original data is used.
#' @param random logical, if TRUE, the random effects are included in the predictions, fitted values and residuals (conditional residuals). For simulate, if TRUE, new random effects are drawn for each simulation, otherwise the fitted random effects are kept.
#' @param digits number of significant digits in the printed output
#' @param type type of residuals: "pearson" (default), the response residuals divided by the standard deviation of the fitted distribution and multiplied by the square root of the weights; "response", the observed minus the fitted values, as proportions for a binomial response; or "quantile", randomised quantile residuals (Dunn and Smyth, 1996), which are standard normal under the model. The quantile residuals do not use the weights, and are random for poisson and binomial responses; use set.seed for reproducible values.
#' @param nsim number of response vectors to simulate
#' @param seed an optional seed for the random number generator
#' @param ... anything else
#'
#' @return Estimated model coefficients, variance-covariance matrix, or predictions.
#' coef and vcov include the variance parameters: for each random effect term the log standard deviations
#' (log_sd:b1|group) and, for correlated random effects, the parameters of the correlation matrix
#' (corr_par:b2,b1|group, unconstrained elements of the lower triangle of its Cholesky factor, by row);
#' and the log residual standard deviation (log_sigma, gaussian) or log precision (log_phi, beta).
#' summary returns a summary.drmTMB object with the coefficient table (Estimate, Std. Error), the
#' variance components, log-likelihood, AIC and BIC, sample sizes, a Pearson dispersion statistic for
#' poisson and binomial models, and convergence information. logLik returns the log-likelihood at the
#' estimates, with the number of free coefficients and variance parameters as degrees of freedom.
#' simulate returns a data frame with nsim columns of simulated responses at the estimated
#' parameters; for a binomial response given as cbind(successes, failures), the number of successes.
#' @export
#'
#' @name drmTMB-methods
#' @aliases print.drmTMB summary.drmTMB print.summary.drmTMB logLik.drmTMB nobs.drmTMB fitted.drmTMB residuals.drmTMB predict.drmTMB coef.drmTMB vcov.drmTMB simulate.drmTMB
#' 
#' @rdname drmTMB-methods
print.drmTMB <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  if (!is.null(x$call)){
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("Model:  ", x$model, "\n", sep="")
  cat("Family: ", x$family, ", ", x$link, " link for the asymptotes\n", sep="")

  cat("\nCoefficients:\n")
  print.default(format(coef(x)[seq_len(sum(!x$fix))], digits=digits), print.gap=2L, quote=FALSE)
  printfixed(userfixed(x), digits)
  printre(recov(x), digits)
  printresid(x$family, residpar(x), digits)
  invisible(x)
}

#' @rdname drmTMB-methods
summary.drmTMB <- function(object, ...){
  dat <- object$obj$env$data
  mod <- modelindex(object$model)
  vcs <- recov(object)
  ll <- logLik(object)

  # coefficient table of the free coefficients
  ests <- object$ssdrl
  coefs <- ests[rownames(ests) %in% paste("b", 1:5, sep=""), , drop=FALSE]
  rownames(coefs) <- object$coefnames[!object$fix]

  # direction of the curve from the sign of the steepness b3 at the observations
  npi <- rep(1:5, sapply(object$Xs, ncol))
  direction <- NULL
  if (!(mod %in% 6:8)){
    b <- object$start
    b[!object$fix] <- coefs[, 1]
    s3 <- sign(object$Xs[[3]] %*% b[npi == 3])
    direction <- if (all(s3 > 0)) "increasing" else if (all(s3 < 0)) "decreasing" else "increasing and decreasing"
  }

  # standard deviations at the boundary (singular fit): the likelihood is flat in the log standard
  # deviation, which gives a very large or no standard error, independent of the scale of the parameter
  ses <- ests[rownames(ests) == "theta", 2]
  toff <- 0
  singular <- unlist(lapply(vcs, function(vc){
    np <- length(vc$pars)
    se <- ses[toff + seq_len(np)]
    toff <<- toff + np + np*(np - 1)/2
    nm <- paste(vc$group, paste("b", vc$pars, sep=""))
    nm[!is.finite(se) | se > 100 | vc$sd < 1e-8]
  }))

  # Pearson dispersion for the families without a dispersion parameter, conditional on the
  # random effects; residual degrees of freedom from the number of estimated parameters
  dispersion <- NULL
  if (object$family %in% c("poisson", "binomial")){
    mu <- predict(object)
    if (object$family == "poisson"){
      r2 <- (dat$y - mu)^2 / mu
    } else {
      r2 <- (dat$y - dat$bn*mu)^2 / (dat$bn*mu*(1 - mu))
    }
    rdf <- nobs(object) - attr(ll, "df")
    chisq <- sum(dat$w * r2)
    dispersion <- c(chisq=chisq, df=rdf, ratio=chisq/rdf)
  }

  # number of levels of each grouping factor, once per factor
  ngroups <- sapply(object$reterms, function(tm) length(tm$levels))
  names(ngroups) <- sapply(object$reterms, function(tm) tm$group)
  ngroups <- ngroups[!duplicated(names(ngroups))]

  gr <- object$obj$gr(object$opt$par)
  res <- list(call=object$call, model=object$model, family=object$family, link=object$link,
              direction=direction, b5=any(!object$fix[npi == 5]),
              nobs=nobs(object), ntrials=if (object$family == "binomial") sum(dat$bn), ngroups=ngroups,
              logLik=ll, AIC=AIC(object), BIC=BIC(object),
              coefficients=coefs, fixed=userfixed(object), varcor=vcs, residpar=residpar(object),
              singular=singular, dispersion=dispersion,
              convergence=object$opt$convergence, message=object$opt$message,
              maxgrad=max(abs(gr)), pdHess=object$sdrl$pdHess)
  class(res) <- "summary.drmTMB"
  res
}

#' @rdname drmTMB-methods
print.summary.drmTMB <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  if (!is.null(x$call)){
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("Model:  ", x$model, if (!is.null(x$direction)) paste(",", x$direction, "curve"), "\n", sep="")
  cat("Family: ", x$family, ", ", x$link, " link for the asymptotes\n", sep="")

  cat("\nNumber of observations: ", x$nobs, sep="")
  if (!is.null(x$ntrials)) cat(", number of trials: ", x$ntrials, sep="")
  if (length(x$ngroups) > 0) cat(", groups: ", paste(names(x$ngroups), x$ngroups, sep=" ", collapse="; "), sep="")
  cat("\n")
  cat("Log-likelihood: ", format(as.numeric(x$logLik), digits=digits), " (df = ", attr(x$logLik, "df"), ")",
      ", AIC: ", format(x$AIC, digits=digits), ", BIC: ", format(x$BIC, digits=digits), "\n", sep="")

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits=digits)
  if (x$b5) cat("The asymmetry parameter b5 is on the log scale.\n")
  printfixed(x$fixed, digits)
  printre(x$varcor, digits)
  if (length(x$singular) > 0){
    cat("Standard deviations close to zero (singular fit): ", paste(x$singular, collapse=", "), "\n", sep="")
  }
  printresid(x$family, x$residpar, digits)
  if (!is.null(x$dispersion)){
    cat("\nPearson dispersion: ", format(x$dispersion["ratio"], digits=digits),
        " (chi-square ", format(x$dispersion["chisq"], digits=digits), " on ", x$dispersion["df"], " df)\n", sep="")
  }

  cat("\nOptimiser: ", x$message, ", maximum absolute gradient ", format(x$maxgrad, digits=2), "\n", sep="")
  if (!isTRUE(x$pdHess)) cat("Warning: the Hessian is not positive definite; the standard errors are not reliable.\n")
  invisible(x)
}

#' @rdname drmTMB-methods
logLik.drmTMB <- function(object, ...){
  # degrees of freedom: free coefficients and variance parameters
  structure(-object$opt$objective, df=length(object$opt$par), nobs=nobs(object), class="logLik")
}

#' @rdname drmTMB-methods
nobs.drmTMB <- function(object, ...){
  length(object$obj$env$data$y)
}

#' @rdname drmTMB-methods
fitted.drmTMB <- function(object, random=TRUE, ...){
  predict(object, random=random)
}

#' @rdname drmTMB-methods
residuals.drmTMB <- function(object, type=c("pearson", "response", "quantile"), random=TRUE, ...){
  type <- match.arg(type)
  dat <- object$obj$env$data
  mu <- predict(object, random=random)
  cnt <- as.vector(dat$y)
  n <- as.vector(dat$bn)
  phi <- residpar(object)
  # binomial: on the scale of proportions
  y <- if (object$family == "binomial") cnt/n else cnt

  if (type == "response") return(y - mu)

  if (type == "pearson"){
    v <- switch(object$family,
                gaussian=phi^2,
                poisson=mu,
                binomial=mu*(1 - mu)/n,
                beta=mu*(1 - mu)/(1 + phi))
    return((y - mu)/sqrt(v) * sqrt(as.vector(dat$w)))
  }

  # randomised quantile residuals: normal quantiles of the fitted distribution function at the
  # observations, drawn uniformly between its left and right limits for the discrete families
  if (object$family == "gaussian") return((y - mu)/phi)
  if (object$family == "beta") return(qnorm(pbeta(y, mu*phi, (1 - mu)*phi)))
  if (object$family == "poisson"){
    lo <- ppois(cnt - 1, mu)
    hi <- ppois(cnt, mu)
  } else {
    lo <- pbinom(cnt - 1, n, mu)
    hi <- pbinom(cnt, n, mu)
  }
  qnorm(runif(length(y), lo, hi))
}

#' @rdname drmTMB-methods
predict.drmTMB <- function(object, newdata=NULL, random=TRUE, ...){
  if (is.null(newdata)) newdata <- object$data
  form <- object$formula
  form[2] <- NULL
  mf <- model.frame(form, data=newdata)
  x <- mf[,1]
  Xs <- lapply(object$fform, function(x) model.matrix(x, data=newdata))
  npars <- sapply(Xs, ncol)
  npi <- rep(1:length(npars), npars)
  ests <- object$ssdrl
  rn <- rownames(ests)
  coefs <- ests[rn %in% paste("b", 1:5, sep=""),,drop=FALSE]
  cest <- object$start
  cest[!object$fix] <- coefs[,1]
  pxl <- lapply(1:5, function(i) Xs[[i]] %*% cest[npi == i])
  
  if (random && !is.null(object$rform)){
    # random effects of the fitted levels; levels not in the fitted data get zero
    re <- ests[rn %in% "u", 1]
    uoff <- 0
    for (tm in object$reterms){
      nl <- length(tm$levels)
      np <- length(tm$pars)
      U <- matrix(re[uoff + seq_len(nl*np)], nrow=nl)
      uoff <- uoff + nl*np
      rt <- mkReTrms(findbars(tm$formula), model.frame(subbars(tm$formula), data=newdata))
      Zt <- rt$Ztlist[[tm$name]]
      lid <- match(rownames(Zt), tm$levels)
      Znew <- t(Zt[!is.na(lid), , drop=FALSE])
      for (k in seq_len(np)){
        pxl[[tm$pars[k]]] <- pxl[[tm$pars[k]]] + as.vector(Znew %*% U[lid[!is.na(lid)], k])
      }
    }
  }

  mod <- which(object$model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
  lnk <- which(object$link == c("identity", "log", "logit"))
  fx <- drcmeanR(mod, x, linkfun(lnk, pxl[[1]], inverse=TRUE), linkfun(lnk, pxl[[2]], inverse=TRUE),
                 pxl[[3]], pxl[[4]], pxl[[5]])
  return(as.vector(fx))
}

#' @rdname drmTMB-methods
coef.drmTMB <- function(object, ...){
    est <- object$estimates
    names(est) <- parnames(object)
    return(est)
}

#' @rdname drmTMB-methods
vcov.drmTMB <- function(object, ...){
    vc <- object$sdrl$cov.fixed
    rownames(vc) <- colnames(vc) <- names(coef(object))
    return(vc)
}

#' @rdname drmTMB-methods
simulate.drmTMB <- function(object, nsim=1, seed=NULL, random=TRUE, ...){
  if (!is.null(seed)) set.seed(seed)
  obj <- object$obj
  par <- obj$env$last.par.best
  if (!is.null(object$rform) && !random){
    # keep the fitted random effects: rebuild the objective with simre = 0
    dat <- obj$env$data
    dat$simre <- 0
    obj <- MakeADFun(dat, obj$env$parList(par=par), map=obj$env$map, random="u", DLL="drcTMB_TMBExports", silent=TRUE)
  }
  sims <- lapply(seq_len(nsim), function(i) obj$simulate(par)$y)
  nbad <- sum(sapply(sims, function(y) any(!is.finite(y))))
  if (nbad > 0) warning(nbad, " of ", nsim, " simulations contain non-finite values; simulated random effects may fall outside the valid parameter range (e.g., b4 <= 0).")
  sims <- as.data.frame(sims, col.names=paste("sim", seq_len(nsim), sep="_"))
  attr(sims, "seed") <- seed
  sims
}


# standard deviations, correlation and covariance matrix of each random effect term
recov <- function(object){
  ests <- object$ssdrl
  rn <- rownames(ests)
  thetas <- ests[rn %in% "theta", 1]
  corrs <- ests[rn %in% "corr", 1]
  toff <- 0
  coff <- 0
  lapply(object$reterms, function(tm){
    np <- length(tm$pars)
    sdm <- exp(thetas[toff + seq_len(np)])
    cr <- matrix(corrs[coff + seq_len(np*np)], nrow=np)
    toff <<- toff + np + np*(np-1)/2
    coff <<- coff + np*np
    list(pars=tm$pars, group=tm$group, sd=sdm, corr=cr, cov=diag(sdm, nrow=np) %*% cr %*% diag(sdm, nrow=np))
  })
}


# index of the model, as in the TMB template
modelindex <- function(model){
  which(model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
}

# coefficients fixed by the user, named; the parameters the model does not use are always fixed
userfixed <- function(object){
  npi <- as.integer(sub("^b([1-5]):.*", "\\1", object$coefnames))
  ufix <- object$fix & !(npi %in% unusedpars(modelindex(object$model)))
  setNames(object$start[ufix], object$coefnames[ufix])
}

printfixed <- function(fixed, digits){
  if (length(fixed) > 0){
    cat("Fixed: ", paste(names(fixed), "=", format(fixed, digits=digits), collapse=", "), "\n", sep="")
  }
}

# random effects: one row per parameter of each term, with standard deviation and correlations
printre <- function(vcs, digits){
  if (length(vcs) == 0) return(invisible())
  rows <- lapply(vcs, function(vc){
    np <- length(vc$pars)
    cr <- matrix("", np, max(np - 1, 0))
    for (i in seq_len(np)[-1]) cr[i, seq_len(i - 1)] <- format(vc$corr[i, seq_len(i - 1)], digits=2, nsmall=2)
    cbind(Groups=c(vc$group, rep("", np - 1)), Name=paste("b", vc$pars, sep=""),
          Std.Dev.=format(unname(vc$sd), digits=digits), cr)
  })
  ncr <- max(sapply(rows, ncol)) - 3
  rows <- lapply(rows, function(r) cbind(r, matrix("", nrow(r), ncr - (ncol(r) - 3))))
  tab <- do.call(rbind, rows)
  colnames(tab) <- c("Groups", "Name", "Std.Dev.", if (ncr > 0) c("Corr", rep("", ncr - 1)))
  rownames(tab) <- rep("", nrow(tab))
  cat("\nRandom effects:\n")
  print(tab, quote=FALSE, right=FALSE)
}

# residual standard deviation (gaussian) or precision (beta); NULL for poisson and binomial
residpar <- function(object){
  if (!(object$family %in% c("gaussian", "beta"))) return(NULL)
  unname(exp(object$sdrl$par.fixed["log_sigma"]))
}

printresid <- function(family, residpar, digits){
  if (family == "gaussian") cat("\nResidual standard deviation: ", format(residpar, digits=digits), "\n", sep="")
  if (family == "beta") cat("\nPrecision: ", format(residpar, digits=digits), "\n", sep="")
}

# names of the estimated parameters: the free coefficients; for each random effect term the log
# standard deviations (log_sd:b1|grp) and the parameters of the correlation matrix
# (corr_par:b2,b1|grp, for the pairs of the lower triangle by row, on the unconstrained scale of TMB);
# the log residual standard deviation (gaussian) or log precision (beta)
parnames <- function(object){
  renames <- unlist(lapply(object$reterms, function(tm){
    pn <- paste("b", tm$pars, sep="")
    np <- length(pn)
    pairs <- unlist(lapply(seq_len(np)[-1], function(i) paste(pn[i], pn[seq_len(i - 1)], sep=",")))
    c(paste("log_sd:", pn, "|", tm$group, sep=""),
      if (np > 1) paste("corr_par:", pairs, "|", tm$group, sep=""))
  }))
  c(object$coefnames[!object$fix], renames,
    switch(object$family, gaussian="log_sigma", beta="log_phi"))
}
