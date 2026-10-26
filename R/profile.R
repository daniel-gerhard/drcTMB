#' Likelihood profiles of a drmTMB model
#'
#' Profiles the (marginal) log-likelihood for single parameters by fixing the parameter at a sequence
#' of values and re-optimising all other parameters, using \code{\link[TMB]{tmbprofile}}. With random
#' effects, the profiled likelihood is the Laplace approximation of the marginal likelihood.
#'
#' The profiles are computed on the scale of \code{coef}: log standard deviations of the random effects
#' (log_sd), parameters of their correlation matrix (corr_par), and the log residual standard deviation
#' (log_sigma) or log precision (log_phi). The confint method transforms the intervals of the log
#' standard deviations and log precision to the natural scale; correlation parameters are only
#' transformed for terms with two random effects, as their parameter determines the correlation
#' only then. If the profile of a log standard deviation does not reach the cutoff at low values,
#' the likelihood is flat towards zero and the lower limit of the standard deviation is 0.
#' Bounds given by \code{lower} and \code{upper} in \code{\link{drmTMB}} are used for the profiled
#' coefficient, but not when re-optimising the other coefficients.
#'
#' @param fitted a drmTMB object
#' @param parm parameters to profile, names or positions in \code{coef(fitted)}. Default are the
#'   estimated coefficients b1-b5, without the variance parameters.
#' @param level the highest confidence level the profiles are computed for: each profile extends until
#'   the deviance from the maximum exceeds the quantile of the chi-square distribution with one degree
#'   of freedom. For confint and plot, the confidence level of the intervals and of the cutoff line.
#' @param trace logical, print the progress of \code{\link[TMB]{tmbprofile}}
#' @param progress logical, show a progress bar over the profiled parameters
#' @param object,x a profile.drmTMB object
#' @param ... further arguments passed to \code{\link[TMB]{tmbprofile}}, e.g. \code{ystep}; ignored by
#'   the confint and plot methods
#'
#' @return profile returns a profile.drmTMB object, a list with a data frame for each parameter with
#'   the parameter values (\code{par}) and the deviance from the maximum (\code{deviance}), i.e. twice
#'   the difference of the negative log-likelihoods. confint returns a matrix of confidence intervals.
#'   plot returns a ggplot object (package ggplot2) with a panel for each parameter, showing the
#'   deviance, the quantile of the chi-square distribution for \code{level} (dashed), and the
#'   estimate (dotted).
#' @export
#'
#' @name profile.drmTMB
#' @aliases profile.drmTMB confint.profile.drmTMB plot.profile.drmTMB
profile.drmTMB <- function(fitted, parm=NULL, level=0.99, trace=FALSE, progress=interactive(), ...){
  est <- coef(fitted)
  pn <- names(est)
  nb <- sum(!fitted$fix)
  if (is.null(parm)) parm <- seq_len(nb)
  if (is.character(parm)){
    id <- match(parm, pn)
    if (any(is.na(id))) stop("unknown parameters in parm: ", paste(parm[is.na(id)], collapse=", "))
    parm <- id
  }
  if (any(parm < 1 | parm > length(est))) stop("parm must select estimated parameters of coef(fitted).")

  # range of the profiled parameter: the bounds of drmTMB (infinite for the variance parameters)
  lower <- if (is.null(fitted$lower)) rep(-Inf, length(est)) else fitted$lower
  upper <- if (is.null(fitted$upper)) rep(Inf, length(est)) else fitted$upper

  nllmin <- fitted$opt$objective
  ytol <- qchisq(level, 1)/2 + 0.1
  pb <- progressbar(length(parm), progress)
  prof <- lapply(parm, function(i){
    p <- TMB::tmbprofile(fitted$obj, name=i, ytol=ytol, parm.range=c(lower[i], upper[i]), trace=trace, ...)
    d <- data.frame(par=p[[1]], deviance=2*(p$value - nllmin))
    pb$tick()
    d[order(d$par), ]
  })
  pb$close()
  names(prof) <- pn[parm]
  attr(prof, "estimates") <- est[parm]
  attr(prof, "corr2") <- corr2names(fitted)
  class(prof) <- "profile.drmTMB"
  prof
}

#' @rdname profile.drmTMB
#' @export
confint.profile.drmTMB <- function(object, parm=NULL, level=0.95, ...){
  if (is.null(parm)) parm <- names(object)
  if (is.numeric(parm)) parm <- names(object)[parm]
  est <- attr(object, "estimates")
  cutoff <- qchisq(level, 1)
  # limits: the parameter values at which the deviance reaches the cutoff, interpolated on each
  # side of the estimate; NA if the profile does not reach the cutoff
  logscale <- attr(object, "logscale")
  if (is.null(logscale)) logscale <- setNames(rep(FALSE, length(object)), names(object))
  limit <- function(d, e, side, nm){
    ds <- d[is.finite(d$deviance) & (if (side == "lower") d$par <= e else d$par >= e), ]
    ds <- ds[order(abs(ds$par - e)), ]
    k <- which(ds$deviance >= cutoff)[1]
    if (is.na(k) || k == 1) return(NA)
    # the square root of the deviance is close to linear in the parameter, on the scale the
    # profile was computed on: spline through the points up to the cutoff (and one beyond), or
    # linear interpolation where the profile is not monotone
    tr <- if (logscale[[nm]]) log else identity
    itr <- if (logscale[[nm]]) exp else identity
    z <- sqrt(pmax(ds$deviance, 0))
    x <- tr(ds$par)
    use <- seq_len(min(k + 1, nrow(ds)))
    if (all(diff(z[use]) > 0)){
      itr(splinefun(z[use], x[use], method="monoH.FC")(sqrt(cutoff)))
    } else {
      itr(approx(z[c(k - 1, k)], x[c(k - 1, k)], xout=sqrt(cutoff))$y)
    }
  }
  ci <- t(sapply(parm, function(nm) c(limit(object[[nm]], est[[nm]], "lower", nm), limit(object[[nm]], est[[nm]], "upper", nm))))
  rn <- parm

  # a flat profile towards small log standard deviations: lower limit 0 for the standard deviation
  logsd <- grepl("^log_sd:", rn)
  ci[logsd & is.na(ci[, 1]), 1] <- -Inf
  if (any(is.na(ci))){
    warning("the profiles of ", paste(rn[apply(is.na(ci), 1, any)], collapse=", "),
            " do not reach the confidence level; compute them with a higher level in profile.", call.=FALSE)
  }

  natscale(ci, rn, attr(object, "corr2"), level)
}

# confidence intervals on the natural scale: standard deviations and precision from their logs,
# correlations of terms with two random effects from their parameter; rows named accordingly
natscale <- function(ci, rn, corr2, level){
  ci <- matrix(ci, ncol=2)
  tolog <- grepl("^log_sd:", rn) | rn %in% c("log_sigma", "log_phi")
  ci[tolog, ] <- exp(ci[tolog, ])
  rn[tolog] <- sub("^log_", "", rn[tolog])
  two <- rn %in% corr2
  ci[two, ] <- ci[two, ] / sqrt(1 + ci[two, ]^2)
  rn[two] <- sub("^corr_par:", "corr:", rn[two])
  dimnames(ci) <- list(rn, paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%"))
  ci
}

# correlation parameters of terms with two random effects, which determine the correlation
corr2names <- function(fit){
  unlist(lapply(fit$reterms, function(tm){
    if (length(tm$pars) == 2) paste("corr_par:b", tm$pars[2], ",b", tm$pars[1], "|", tm$group, sep="")
  }))
}

#' @rdname profile.drmTMB
#' @export
plot.profile.drmTMB <- function(x, level=0.95, ...){
  if (!requireNamespace("ggplot2", quietly=TRUE)) stop("the plot method needs the package ggplot2.")
  nm <- names(x)
  d <- do.call(rbind, lapply(nm, function(n) data.frame(parameter=n, par=x[[n]]$par, deviance=x[[n]]$deviance)))
  d$parameter <- factor(d$parameter, levels=nm)
  est <- data.frame(parameter=factor(nm, levels=nm), par=unname(attr(x, "estimates")))
  ggplot2::ggplot(d, ggplot2::aes(x=par, y=deviance)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept=qchisq(level, 1), linetype="dashed") +
    ggplot2::geom_vline(data=est, ggplot2::aes(xintercept=par), linetype="dotted") +
    ggplot2::facet_wrap(~ parameter, scales="free_x") +
    ggplot2::labs(x=NULL, y="Deviance")
}

# columns of the profile data used in the aesthetics of plot.profile.drmTMB
utils::globalVariables(c("par", "deviance"))

#' Likelihood profiles of effective doses
#'
#' Profiles the (marginal) log-likelihood for effective doses estimated by \code{\link{ed}}. At each
#' value of an effective dose, the model is reparameterised: the effective dose is a function of b3,
#' b4 and b5 that can be solved for b4, so one coefficient of b4 is computed from the effective dose
#' and the other coefficients, which are re-optimised.
#'
#' The profiles are for the effective doses of the population curves, without the random effects;
#' the variation between groups (\code{random=TRUE} in \code{\link{ed}}) and the residual variation are
#' not included. For the models on the log-dose scale and the models without b3, the effective doses
#' are profiled on the log scale.
#'
#' @param fitted an edTMB object
#' @param parm effective doses to profile, names or positions in \code{coef(fitted)}; default all
#' @param level the highest confidence level the profiles are computed for
#' @param maxit maximum number of steps in each direction
#' @param progress logical, show a progress bar over the profiled effective doses
#' @param ... ignored
#'
#' @return a profile.edTMB object, which has the confint and plot methods of
#'   \code{\link{profile.drmTMB}}
#' @export
profile.edTMB <- function(fitted, parm=NULL, level=0.99, maxit=100, progress=interactive(), ...){
  fit <- fitted$fit
  if (is.null(fit)) stop("compute the effective doses again with ed() to profile them.")
  est <- coef(fitted)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown effective doses in parm.")

  obj <- fit$obj
  # keep the state of the TMB object, as tmbprofile does
  keep <- c("last.par.best", "random.start", "value.best", "last.par", "inner.control", "tracemgc")
  old <- mget(keep, envir=obj$env, ifnotfound=list(NULL))
  on.exit(for (v in keep) assign(v, old[[v]], envir=obj$env))
  obj$env$inner.control$trace <- FALSE
  obj$env$tracemgc <- FALSE

  mod <- modelindex(fit$model)
  lnk4 <- b4linkindex(fit)
  logscale <- !(mod %in% c(1, 9))
  free <- which(!fit$fix)
  nf <- length(free)
  par0 <- fit$opt$par
  nllmin <- fit$opt$objective
  cutoff <- qchisq(level, 1)
  Hfull <- tryCatch(solve(fit$sdrl$cov.fixed), error=function(e) NULL)

  pb <- progressbar(length(parm), progress)
  prof <- lapply(parm, function(i){
    L <- fitted$linfct[[i]]
    q <- fitted$q[[i]]
    # the coefficient of b4 computed from the effective dose: free, with the largest weight
    w4 <- L[4, free]
    if (all(w4 == 0)) stop("the effective dose ", names(est)[i], " does not depend on an estimated coefficient of b4.")
    kp <- which.max(abs(w4))
    k <- free[kp]
    if (L[3, k] != 0 || L[5, k] != 0) stop("the coefficient ", fit$coefnames[k], " of b4 must not enter b3 or b5 in linfct.")

    # full parameter vector from the reduced one (without coefficient k) at effective dose ed
    pvfull <- function(par){
      pv <- fit$start
      pv[free] <- par[seq_len(nf)]
      pv
    }
    # b4 on the scale of its coefficients
    b4fun <- function(ed, b3, b5) linkfun(lnk4, edb4(mod, q, ed, b3, b5))
    expand <- function(x, ed){
      par <- append(x, 0, after=kp - 1)
      pv <- pvfull(par)
      lpv <- as.vector(L %*% pv)
      b4 <- b4fun(ed, lpv[3], lpv[5])
      par[kp] <- (b4 - sum(L[4, -k] * pv[-k])) / L[4, k]
      list(par=par, lpv=lpv)
    }
    # derivatives of the full parameter vector w.r.t. the reduced one: coefficient k depends on the
    # other coefficients of b4 and, through the effective dose, on b3 and b5
    jacobian <- function(e, ed){
      h <- 1e-6
      d3 <- (b4fun(ed, e$lpv[3] + h, e$lpv[5]) - b4fun(ed, e$lpv[3] - h, e$lpv[5])) / (2*h)
      d5 <- (b4fun(ed, e$lpv[3], e$lpv[5] + h) - b4fun(ed, e$lpv[3], e$lpv[5] - h)) / (2*h)
      J <- diag(length(par0))[, -kp, drop=FALSE]
      J[kp, seq_len(nf - 1)] <- ((d3*L[3, free] + d5*L[5, free] - L[4, free]) / L[4, k])[-kp]
      J
    }
    fn <- function(x, ed) obj$fn(expand(x, ed)$par)
    gr <- function(x, ed){
      e <- expand(x, ed)
      as.vector(obj$gr(e$par) %*% jacobian(e, ed))
    }
    # the Hessian at the estimates, transformed to the reduced parameters, as a fixed approximation:
    # far fewer iterations than quasi-Newton steps from scratch
    he <- if (is.null(Hfull)) NULL else function(x, ed){
      J <- jacobian(expand(x, ed), ed)
      t(J) %*% Hfull %*% J
    }

    # steps of a fifth of the standard error, on the log scale for positive effective doses
    e0 <- est[[i]]
    se <- sqrt(fitted$variance[i, i])
    tr <- if (logscale) log else identity
    itr <- if (logscale) exp else identity
    step <- (if (logscale) se/e0 else se) / 5
    side <- function(dir){
      x <- par0[-kp]
      out <- NULL
      for (it in seq_len(maxit)){
        ed <- itr(tr(e0) + dir*it*step)
        o <- suppressWarnings(nlminb(x, fn, gr, he, ed=ed))
        x <- o$par
        out <- rbind(out, data.frame(par=ed, deviance=2*(o$objective - nllmin)))
        if (out$deviance[nrow(out)] > cutoff) break
      }
      out
    }
    d <- rbind(side(-1), data.frame(par=e0, deviance=0), side(1))
    pb$tick()
    d[order(d$par), ]
  })
  pb$close()
  names(prof) <- names(est)[parm]
  attr(prof, "estimates") <- est[parm]
  attr(prof, "logscale") <- setNames(rep(logscale, length(parm)), names(est)[parm])
  class(prof) <- c("profile.edTMB", "profile.drmTMB")
  prof
}

# effective doses at level q of F (residual zero), as drced in the TMB template: b4 times a function
# h of q, b3 and b5, or b4 plus h for logistic and gompertz
edh <- function(mod, q, b3, b5){
  switch(mod,
         -log(q^(-1/exp(b5)) - 1)/b3,       # logistic: ed = b4 + h
         (q^(-1/exp(b5)) - 1)^(-1/b3),      # loglogistic: ed = b4 h
         (-log(q))^(-1/b3),                 # weibull1
         (-log(1 - q))^(1/b3),              # weibull2
         exp(qnorm(q)/b3),                  # lognormal
         q/(1 - q),                         # michaelismenten
         -log(1 - q),                       # asymptoticregression
         -log(1 - q),                       # exponentialdecay
         -log(-log(q))/b3)                  # gompertz
}

# effective dose from b3, b4 and b5
edfun <- function(mod, q, b3, b4, b5){
  h <- edh(mod, q, b3, b5)
  if (mod %in% c(1, 9)) b4 + h else b4 * h
}

# b4 of a curve with effective dose ed, given b3 and b5: the inverse of edfun
edb4 <- function(mod, q, ed, b3, b5){
  h <- edh(mod, q, b3, b5)
  if (mod %in% c(1, 9)) ed - h else ed / h
}

#' Confidence intervals for drmTMB models and effective doses
#'
#' Profile likelihood intervals (default), Wald intervals, or parametric bootstrap intervals for the
#' coefficients and variance parameters of a drmTMB model, or for effective doses estimated by
#' \code{\link{ed}}.
#'
#' Standard deviations of the random effects, the residual standard deviation and the precision are
#' reported on their natural scale; Wald intervals for them are computed on the log scale. The
#' correlation of a term with two random effects is reported as a correlation; the parameters of
#' larger correlation matrices on the unconstrained scale of \code{coef}.
#'
#' Profile and bootstrap intervals for effective doses are for the population curves. Effective doses computed
#' with \code{random=TRUE} for a model with random effects, or with \code{residual=TRUE}, include
#' further variation in their standard errors; for these, Wald intervals are returned with a message.
#' Bootstrap intervals need many fits of the model; with random effects this can take several minutes.
#' Computing profiles can take a while for models with random effects; see \code{\link{profile.drmTMB}}
#' and \code{\link{profile.edTMB}}, whose results can be reused with \code{confint}.
#'
#' @param object a drmTMB or edTMB object
#' @param parm parameters or effective doses, names or positions in \code{coef(object)}; default all
#' @param level the confidence level
#' @param method "profile" for profile likelihood intervals, "wald" for Wald intervals, or "boot" for
#'   percentile intervals of a parametric bootstrap: responses are simulated from the fitted model
#'   (with new random effects, see \code{\link{simulate.drmTMB}}) and the model is fitted again to each,
#'   starting at the estimates. Fits that do not converge are left out, with a message.
#' @param nsim number of bootstrap simulations
#' @param seed an optional seed for the random number generator of the bootstrap
#' @param boot a bootstrap sample of the model from \code{\link{bootstrap}}, used instead of new
#'   simulations; implies \code{method="boot"}
#' @param ... further arguments passed to \code{\link{profile.drmTMB}}, \code{\link{profile.edTMB}} or
#'   \code{\link{bootstrap}}, e.g. \code{progress=FALSE}
#'
#' @return a matrix with the lower and upper confidence limits
#' @export
#' @name confint.drmTMB
#' @aliases confint.drmTMB confint.edTMB
confint.drmTMB <- function(object, parm=NULL, level=0.95, method=c("profile", "wald", "boot"), nsim=500, seed=NULL, boot=NULL, ...){
  method <- if (!is.null(boot)) "boot" else match.arg(method)
  est <- coef(object)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown parameters in parm.")
  if (method == "profile"){
    message("Computing profile confidence intervals ...")
    return(confint(profile(object, parm=parm, level=level, ...), level=level))
  }
  if (method == "boot") return(confint(getboot(object, boot, nsim, seed, ...), parm=parm, level=level))
  se <- sqrt(diag(vcov(object)))[parm]
  z <- qnorm((1 + level)/2)
  natscale(cbind(est[parm] - z*se, est[parm] + z*se), names(est)[parm], corr2names(object), level)
}

#' @rdname confint.drmTMB
#' @export
confint.edTMB <- function(object, parm=NULL, level=0.95, method=c("profile", "wald", "boot"), nsim=500, seed=NULL, boot=NULL, ...){
  method <- if (!is.null(boot)) "boot" else match.arg(method)
  est <- coef(object)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown effective doses in parm.")
  if (method %in% c("profile", "boot")){
    fit <- object$fit
    if (is.null(fit)){
      message("Wald intervals: compute the effective doses again with ed() for ", method, " intervals.")
      method <- "wald"
    } else if ((isTRUE(object$random) && length(fit$reterms) > 0) || isTRUE(object$residual)){
      message("Wald intervals: ", method, " intervals are only available for population effective doses, ",
              "ed(..., random=FALSE, residual=FALSE).")
      method <- "wald"
    }
  }
  if (method == "profile"){
    message("Computing profile confidence intervals ...")
    return(confint(profile(object, parm=parm, level=level, ...), level=level))
  }
  if (method == "boot"){
    fit <- object$fit
    b <- bootconverged(getboot(fit, boot, nsim, seed, ...))
    mod <- modelindex(fit$model)
    nb <- sum(!fit$fix)
    # effective doses of each bootstrap fit
    eds <- sapply(parm, function(i){
      apply(b[, seq_len(nb), drop=FALSE], 1, function(bf){
        pv <- fit$start
        pv[!fit$fix] <- bf
        lpv <- as.vector(object$linfct[[i]] %*% pv)
        edfun(mod, object$q[[i]], lpv[3], linkfun(b4linkindex(fit), lpv[4], inverse=TRUE), lpv[5])
      })
    })
    ci <- t(apply(matrix(eds, ncol=length(parm)), 2, quantile, probs=c((1 - level)/2, (1 + level)/2), na.rm=TRUE, names=FALSE))
    dimnames(ci) <- list(names(est)[parm], paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%"))
    return(ci)
  }
  se <- sqrt(diag(vcov(object)))[parm]
  z <- qnorm((1 + level)/2)
  ci <- cbind(est[parm] - z*se, est[parm] + z*se)
  dimnames(ci) <- list(names(est)[parm], paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%"))
  ci
}
