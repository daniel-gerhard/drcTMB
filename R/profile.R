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
#'   of freedom.
#' @param trace logical, print the progress of \code{\link[TMB]{tmbprofile}}
#' @param object,x a profile.drmTMB object
#' @param ... further arguments passed to \code{\link[TMB]{tmbprofile}}, e.g. \code{ystep}; ignored by
#'   the confint and plot methods
#'
#' @return profile returns a profile.drmTMB object, a list with a data frame for each parameter with
#'   the parameter values (\code{par}) and the deviance from the maximum (\code{deviance}), i.e. twice
#'   the difference of the negative log-likelihoods. confint returns a matrix of confidence intervals.
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
  limit <- function(d, e, side){
    ds <- d[is.finite(d$deviance) & (if (side == "lower") d$par <= e else d$par >= e), ]
    ds <- ds[order(abs(ds$par - e)), ]
    k <- which(ds$deviance >= cutoff)[1]
    if (is.na(k) || k == 1) return(NA)
    approx(ds$deviance[c(k - 1, k)], ds$par[c(k - 1, k)], xout=cutoff)$y
  }
  ci <- t(sapply(parm, function(nm) c(limit(object[[nm]], est[[nm]], "lower"), limit(object[[nm]], est[[nm]], "upper"))))
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
  nm <- names(x)
  op <- par(mfrow=n2mfrow(length(nm)))
  on.exit(par(op))
  for (i in seq_along(nm)){
    d <- x[[i]]
    plot(d$par, d$deviance, type="l", xlab=nm[i], ylab="deviance", ...)
    abline(h=qchisq(level, 1), lty=2)
    abline(v=attr(x, "estimates")[i], lty=3)
  }
  invisible(x)
}
