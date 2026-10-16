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
profile.drmTMB <- function(fitted, parm=NULL, level=0.99, trace=FALSE, ...){
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
  prof <- lapply(parm, function(i){
    p <- TMB::tmbprofile(fitted$obj, name=i, ytol=ytol, parm.range=c(lower[i], upper[i]), trace=trace, ...)
    d <- data.frame(par=p[[1]], deviance=2*(p$value - nllmin))
    d[order(d$par), ]
  })
  names(prof) <- pn[parm]
  attr(prof, "estimates") <- est[parm]
  # correlation parameters of terms with two random effects, which determine the correlation
  attr(prof, "corr2") <- unlist(lapply(fitted$reterms, function(tm){
    if (length(tm$pars) == 2) paste("corr_par:b", tm$pars[2], ",b", tm$pars[1], "|", tm$group, sep="")
  }))
  class(prof) <- "profile.drmTMB"
  prof
}

#' @rdname profile.drmTMB
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

  # natural scale for standard deviations, precision, and correlations of two random effects
  tolog <- logsd | rn %in% c("log_sigma", "log_phi")
  ci[tolog, ] <- exp(ci[tolog, ])
  rn[tolog] <- sub("^log_", "", rn[tolog])
  two <- rn %in% attr(object, "corr2")
  ci[two, ] <- ci[two, ] / sqrt(1 + ci[two, ]^2)
  rn[two] <- sub("^corr_par:", "corr:", rn[two])
  ci <- matrix(ci, ncol=2, dimnames=list(rn, NULL))
  colnames(ci) <- paste(format(100*c((1 - level)/2, (1 + level)/2), trim=TRUE, scientific=FALSE, digits=3), "%")
  ci
}

#' @rdname profile.drmTMB
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
#' @param ... ignored
#'
#' @return a profile.edTMB object, which has the confint and plot methods of
#'   \code{\link{profile.drmTMB}}
#' @export
profile.edTMB <- function(fitted, parm=NULL, level=0.99, maxit=100, ...){
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
  logscale <- !(mod %in% c(1, 9))
  free <- which(!fit$fix)
  nf <- length(free)
  par0 <- fit$opt$par
  nllmin <- fit$opt$objective
  cutoff <- qchisq(level, 1)
  Hfull <- tryCatch(solve(fit$sdrl$cov.fixed), error=function(e) NULL)

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
    b4fun <- function(ed, b3, b5) edb4(mod, q, ed, b3, b5)
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
    d[order(d$par), ]
  })
  names(prof) <- names(est)[parm]
  attr(prof, "estimates") <- est[parm]
  attr(prof, "logscale") <- setNames(rep(logscale, length(parm)), names(est)[parm])
  class(prof) <- c("profile.edTMB", "profile.drmTMB")
  prof
}

# b4 of a curve with effective dose ed at level q of F, given b3 and b5 (residual zero):
# the inverse of drced in the TMB template, which is b4 times or plus a function of q, b3 and b5
edb4 <- function(mod, q, ed, b3, b5){
  h <- switch(mod,
              -log(q^(-1/exp(b5)) - 1)/b3,       # logistic: ed = b4 + h
              (q^(-1/exp(b5)) - 1)^(-1/b3),      # loglogistic: ed = b4 h
              (-log(q))^(-1/b3),                 # weibull1
              (-log(1 - q))^(1/b3),              # weibull2
              exp(qnorm(q)/b3),                  # lognormal
              q/(1 - q),                         # michaelismenten
              -log(1 - q),                       # asymptoticregression
              -log(1 - q),                       # exponentialdecay
              -log(-log(q))/b3)                  # gompertz
  if (mod %in% c(1, 9)) ed - h else ed / h
}
