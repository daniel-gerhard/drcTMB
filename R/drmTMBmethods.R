#' drmTMB Methods
#'
#' @param x,object a drmTMB object
#' @param newdata a data frame for prediction. If NULL, the original data is used.
#' @param random logical, if TRUE, the random effects are included in the predictions, fitted values and residuals (conditional residuals). For simulate, if TRUE, new random effects are drawn for each simulation, otherwise the fitted random effects are kept.
#'   For prediction intervals, random=TRUE gives the interval of a new observation in the group of each row of newdata, and random=FALSE the interval of a new observation in a new group, including the variance of the random effects. Rows of newdata with a level that is not in the fitted data are treated as new groups for that random effect term.
#' @param se.fit logical, if TRUE, predict returns the standard errors of the predicted means, by the delta method; with random=TRUE they include the uncertainty of the predicted random effects.
#' @param interval type of interval of predict: "none", "confidence" for confidence intervals of the mean, computed on the log scale for poisson and on the logit scale for binomial and beta responses, or "prediction" for prediction intervals of a new observation with weight 1. For the gaussian family without boot, the prediction intervals use the normal distribution, and the variance of the random effects of new groups is added by the delta method. Otherwise they are quantiles of new observations simulated at nsim draws from the normal approximation of the fixed and random effects (or at the bootstrap fits), with new random effects for new groups; the variance parameters are kept at their estimates; for a binomial response as proportions, with the number of trials from the response in newdata (one for a binary response).
#' @param level confidence level of the intervals of predict
#' @param digits number of significant digits in the printed output
#' @param type type of residuals: "pearson" (default), the response residuals divided by the standard deviation of the fitted distribution and multiplied by the square root of the weights; "response", the observed minus the fitted values, as proportions for a binomial response; or "quantile", randomised quantile residuals (Dunn and Smyth, 1996), which are standard normal under the model. The quantile residuals do not use the weights, and are random for poisson and binomial responses; use set.seed for reproducible values.
#' @param nsim number of response vectors to simulate; for predict, the number of parameter draws for simulated prediction intervals
#' @param seed an optional seed for the random number generator
#' @param boot a bootstrap sample of the model from \code{\link{bootstrap}}. predict then gives the standard deviation of the bootstrap predictions as se.fit, percentile confidence intervals, and prediction intervals simulated at the bootstrap fits. Only for population predictions (random=FALSE) of models with random effects.
#' @param ... anything else
#'
#' @return Estimated model coefficients, variance-covariance matrix, or predictions.
#' predict returns a vector of predicted means, or, with se.fit=TRUE or an interval, a data frame with
#' the predicted means (fit), their standard errors (se.fit), and the limits of the intervals (lwr, upr).
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
#' @export
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
#' @export
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
#' @export
logLik.drmTMB <- function(object, ...){
  # degrees of freedom: free coefficients and variance parameters
  structure(-object$opt$objective, df=length(object$opt$par), nobs=nobs(object), class="logLik")
}

#' @rdname drmTMB-methods
#' @export
nobs.drmTMB <- function(object, ...){
  length(object$obj$env$data$y)
}

#' @rdname drmTMB-methods
#' @export
fitted.drmTMB <- function(object, random=TRUE, ...){
  predict(object, random=random)
}

#' @rdname drmTMB-methods
#' @export
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
#' @export
predict.drmTMB <- function(object, newdata=NULL, random=TRUE, se.fit=FALSE, interval=c("none", "confidence", "prediction"), level=0.95, boot=NULL, nsim=1000, seed=NULL, ...){
  interval <- match.arg(interval)
  if (is.null(newdata)) newdata <- object$data
  pd <- preddata(object, newdata, random)
  fit <- predmean(pd, pd$eta)
  if (!se.fit && interval == "none") return(fit)
  if (!is.null(boot)){
    boot <- getboot(object, boot, nsim, seed)
    if (random && length(object$reterms) > 0){
      message("Bootstrap intervals are only available for population predictions, random=FALSE; ",
              "using the delta method and simulation from the estimates.")
      boot <- NULL
    }
  }

  # standard errors and confidence intervals: delta method, or the curves of the bootstrap fits
  probs <- c((1 - level)/2, (1 + level)/2)
  if (is.null(boot)){
    se <- predictse(object, pd)
  } else {
    bc <- bootcurves(object, pd, boot)
    se <- apply(bc, 1, sd, na.rm=TRUE)
  }
  res <- data.frame(fit=fit, se.fit=se)
  z <- qnorm((1 + level)/2)
  if (interval == "confidence"){
    res[c("lwr", "upr")] <- if (is.null(boot)) predci(fit, se, object$family, z)
                            else t(apply(bc, 1, quantile, probs=probs, na.rm=TRUE, names=FALSE))
  }

  if (interval == "prediction"){
    if (object$family == "gaussian" && is.null(boot)){
      # new observation: uncertainty of the mean, residual variance, and the variance of the random
      # effects of the terms without a fitted level (new groups, or all terms for random=FALSE)
      sdp <- sqrt(se^2 + residpar(object)^2 + newgroupvar(object, pd))
      res$lwr <- fit - z*sdp
      res$upr <- fit + z*sdp
    } else {
      # simulated new observations at parameter draws: from the bootstrap fits, or from the normal
      # approximation of the fixed and random effects
      trials <- predtrials(object, newdata)
      if (!is.null(seed)) set.seed(seed)
      draws <- if (is.null(boot)) jointdraws(object, nsim) else bootdraws(object, boot)
      res[c("lwr", "upr")] <- predsim(object, pd, draws, trials, probs)
    }
  }
  res
}

# design of the predictions for newdata: doses, design matrices, linear predictors of b1-b5 (with the
# random effects of the fitted levels if random=TRUE), the random effects design Zp for the TMB
# template, and for each term the rows without a fitted level (newgroup)
preddata <- function(object, newdata, random){
  form <- object$formula
  form[2] <- NULL
  x <- model.frame(form, data=newdata, na.action=na.pass)[, 1]
  n <- length(x)
  Xs <- lapply(object$fterms, function(ft) designmatrix(ft, newdata))
  npi <- rep(1:5, sapply(Xs, ncol))
  ests <- object$ssdrl
  rn <- rownames(ests)
  cest <- object$start
  cest[!object$fix] <- ests[rn %in% paste("b", 1:5, sep=""), 1]
  eta <- sapply(1:5, function(i) as.vector(Xs[[i]] %*% cest[npi == i]))
  eta <- matrix(eta, nrow=n)

  re <- ests[rn %in% "u", 1]
  uoff <- 0
  Zl <- list()
  newgroup <- matrix(TRUE, n, length(object$reterms))
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    np <- length(tm$pars)
    U <- matrix(re[uoff + seq_len(nl*np)], nrow=nl)
    uoff <- uoff + nl*np
    Zt <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(n, nl))
    if (random){
      # random effects of the fitted levels; levels not in the fitted data get zero
      rt <- mkReTrms(findbars(tm$formula), model.frame(subbars(tm$formula), data=newdata))
      Zn <- rt$Ztlist[[tm$name]]
      lid <- match(rownames(Zn), tm$levels)
      ok <- which(!is.na(lid))
      P <- Matrix::sparseMatrix(i=ok, j=lid[ok], x=1, dims=c(nrow(Zn), nl))
      Zt <- t(Zn) %*% P
      newgroup[, t] <- Matrix::rowSums(Zt) == 0
      for (k in seq_len(np)) eta[, tm$pars[k]] <- eta[, tm$pars[k]] + as.vector(Zt %*% U[, k])
    }
    Zl[[t]] <- Zt
  }
  Zp <- if (length(Zl) > 0) do.call(cbind, Zl) else Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(n, 0))
  list(x=x, Xs=Xs, eta=eta, Zp=Zp, Zl=Zl, newgroup=newgroup,
       mod=modelindex(object$model), lnk=which(object$link == c("identity", "log", "logit")))
}

# mean function at the linear predictors eta of b1-b5
predmean <- function(pd, eta){
  as.vector(drcmeanR(pd$mod, pd$x, linkfun(pd$lnk, eta[, 1], inverse=TRUE), linkfun(pd$lnk, eta[, 2], inverse=TRUE),
                     eta[, 3], eta[, 4], eta[, 5]))
}

# standard errors of the predicted means by the delta method of TMB: the mean at the new doses is
# reported by the template, at the estimates; with random effects, sdreport includes the uncertainty
# of the predicted random effects
predictse <- function(object, pd){
  env <- object$obj$env
  dat <- env$data
  dat$xp <- pd$x
  for (i in 1:5) dat[[paste("Xp", i, sep="")]] <- pd$Xs[[i]]
  # TMB reads sparse matrices in triplet form
  dat$Zp <- as(as(pd$Zp, "TsparseMatrix"), "generalMatrix")
  rnd <- !is.null(env$random)
  obj <- MakeADFun(dat, env$parList(object$opt$par), map=env$map, random=if (rnd) "u",
                   DLL="drcTMB_TMBExports", silent=TRUE)
  # the likelihood is unchanged: the Hessian of the fit is reused
  H <- if (isTRUE(object$sdrl$pdHess)) solve(object$sdrl$cov.fixed)
  sdr <- sdreport(obj, par.fixed=object$opt$par, hessian.fixed=H, getReportCovariance=FALSE)
  s <- summary(sdr, "report")
  unname(s[rownames(s) == "mu", 2])
}

# population curves of the bootstrap fits, one column per converged fit
bootcurves <- function(object, pd, boot){
  est <- boot$estimates[!is.na(boot$estimates[, 1]), , drop=FALSE]
  npi <- rep(1:5, sapply(pd$Xs, ncol))
  apply(est, 1, function(b){
    cest <- object$start
    cest[!object$fix] <- b[seq_len(sum(!object$fix))]
    eta <- matrix(sapply(1:5, function(i) as.vector(pd$Xs[[i]] %*% cest[npi == i])), ncol=5)
    predmean(pd, eta)
  })
}

# parameter vectors (as obj$env$par, one per column) of the converged bootstrap fits; the random
# effects are not used for population predictions
bootdraws <- function(object, boot){
  est <- boot$estimates[!is.na(boot$estimates[, 1]), , drop=FALSE]
  env <- object$obj$env
  r <- env$random
  apply(est, 1, function(b){
    par <- numeric(length(env$par))
    if (is.null(r)) par[] <- b else par[-r] <- b
    par
  })
}

# draws from the normal approximation of all parameters: the covariance of the fixed parameters,
# or with random effects the joint precision of the fixed and random effects from sdreport
jointdraws <- function(object, nsim){
  env <- object$obj$env
  r <- env$random
  if (is.null(r)){
    m <- object$opt$par
    R <- tryCatch(chol(object$sdrl$cov.fixed), error=function(e) NULL)
    if (is.null(R)) stop("the covariance matrix of the estimates is not positive definite.")
    return(m + t(R) %*% matrix(rnorm(length(m)*nsim), nrow=length(m)))
  }
  sdr <- sdreport(object$obj, par.fixed=object$opt$par, getJointPrecision=TRUE)
  m <- numeric(length(env$par))
  m[-r] <- sdr$par.fixed
  m[r] <- sdr$par.random
  R <- tryCatch(chol(as.matrix(sdr$jointPrecision)), error=function(e) NULL)
  if (is.null(R)) stop("the joint precision matrix of the parameters is not positive definite.")
  m + backsolve(R, matrix(rnorm(length(m)*nsim), nrow=length(m)))
}

# correlation matrix of a random effect term from its unconstrained parameters, as UNSTRUCTURED_CORR
# in TMB: lower triangle of L filled by row, correlation of L L'
corrmatrix <- function(cpar, np){
  L <- diag(np)
  L[upper.tri(L)] <- cpar
  L <- t(L)
  S <- L %*% t(L)
  S / sqrt(diag(S) %o% diag(S))
}

# number of trials of the new observations of a binomial model: from the response in newdata, or
# one for a binary response
predtrials <- function(object, newdata){
  if (object$family != "binomial") return(NULL)
  y <- tryCatch(model.response(model.frame(object$formula, data=newdata, na.action=na.pass)), error=function(e) NULL)
  if (is.matrix(y) && ncol(y) == 2) return(rowSums(y))
  if (all(object$obj$env$data$bn == 1)) return(rep(1, nrow(newdata)))
  stop("prediction intervals of a binomial model need the number of trials: give the response, e.g. cbind(successes, failures), in newdata.")
}

# quantiles of new observations simulated at each parameter draw (columns of draws): conditional on
# the drawn random effects of the fitted levels for random=TRUE, with new random effects for the
# rows without a fitted level; binomial observations as proportions. The variance parameters are kept
# at their estimates: their normal approximation is poor, e.g. for standard deviations close to zero
predsim <- function(object, pd, draws, trials, probs){
  env <- object$obj$env
  n <- length(pd$x)
  pest <- env$parList(object$opt$par)
  sims <- apply(draws, 2, function(par){
    pl <- env$parList(par=par)
    eta <- matrix(sapply(1:5, function(i) as.vector(pd$Xs[[i]] %*% pl[[paste("b", i, sep="")]])), ncol=5)
    uoff <- 0
    toff <- 0
    for (t in seq_along(object$reterms)){
      tm <- object$reterms[[t]]
      nl <- length(tm$levels)
      np <- length(tm$pars)
      U <- matrix(pl$u[uoff + seq_len(nl*np)], nrow=nl)
      sdt <- exp(pest$theta[toff + seq_len(np)])
      Ct <- corrmatrix(pest$theta[toff + np + seq_len(np*(np - 1)/2)], np)
      uoff <- uoff + nl*np
      toff <- toff + np + np*(np - 1)/2
      eta[, tm$pars] <- eta[, tm$pars] + as.matrix(pd$Zl[[t]] %*% U)
      new <- which(pd$newgroup[, t])
      if (length(new) > 0){
        E <- matrix(rnorm(length(new)*np), ncol=np) %*% chol(Ct) %*% diag(sdt, nrow=np)
        eta[new, tm$pars] <- eta[new, tm$pars] + E
      }
    }
    mu <- predmean(pd, eta)
    phi <- exp(pest$log_sigma)
    # draws outside the valid range of the parameters give NA
    suppressWarnings(switch(object$family,
                            gaussian=rnorm(n, mu, phi),
                            poisson=rpois(n, mu),
                            binomial=rbinom(n, trials, mu)/trials,
                            beta=rbeta(n, mu*phi, (1 - mu)*phi)))
  })
  sims <- matrix(sims, nrow=n)
  bad <- mean(!is.finite(sims))
  if (bad > 0) warning(format(100*bad, digits=2), "% of the simulated observations are not finite and are left out; ",
                       "parameter draws may fall outside the valid parameter range (e.g., b4 <= 0).", call.=FALSE)
  # discrete responses: quantiles of the simulated values
  type <- if (object$family %in% c("poisson", "binomial")) 1 else 7
  t(apply(sims, 1, quantile, probs=probs, na.rm=TRUE, names=FALSE, type=type))
}

# confidence intervals of the means, computed on the log scale for poisson and the logit scale for
# binomial and beta, so that they stay inside the range of the mean
predci <- function(fit, se, family, z){
  tr <- switch(family,
               gaussian=list(f=identity, d=function(m) 1, i=identity),
               poisson=list(f=log, d=function(m) 1/m, i=exp),
               list(f=qlogis, d=function(m) 1/(m*(1 - m)), i=plogis))
  eta <- tr$f(fit)
  se.eta <- se*tr$d(fit)
  ci <- cbind(lwr=tr$i(eta - z*se.eta), upr=tr$i(eta + z*se.eta))
  # means at the boundary, e.g. an asymptote fixed at 0, without uncertainty
  pt <- se == 0 & !is.na(se)
  ci[pt, ] <- fit[pt]
  ci
}

# variance of the mean from the random effects of the terms without a fitted level in each row, by
# the delta method with numerical derivatives of the mean w.r.t. the linear predictors of b1-b5
newgroupvar <- function(object, pd){
  n <- length(pd$x)
  v <- numeric(n)
  if (length(object$reterms) == 0 || !any(pd$newgroup)) return(v)
  g <- sapply(1:5, function(k){
    h <- 1e-5*pmax(1, abs(pd$eta[, k]))
    e1 <- e2 <- pd$eta
    e1[, k] <- e1[, k] + h
    e2[, k] <- e2[, k] - h
    (predmean(pd, e1) - predmean(pd, e2)) / (2*h)
  })
  g <- matrix(g, nrow=n)
  vcs <- recov(object)
  for (t in seq_along(vcs)){
    gt <- g[, vcs[[t]]$pars, drop=FALSE]
    vt <- rowSums((gt %*% vcs[[t]]$cov) * gt)
    v <- v + ifelse(pd$newgroup[, t], vt, 0)
  }
  v
}

#' @rdname drmTMB-methods
#' @export
coef.drmTMB <- function(object, ...){
    est <- object$estimates
    names(est) <- parnames(object)
    return(est)
}

#' @rdname drmTMB-methods
#' @export
vcov.drmTMB <- function(object, ...){
    vc <- object$sdrl$cov.fixed
    rownames(vc) <- colnames(vc) <- names(coef(object))
    return(vc)
}

#' @rdname drmTMB-methods
#' @export
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
