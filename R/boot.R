#' Parametric bootstrap of a drmTMB model
#'
#' Simulates responses from the fitted model, with new random effects (see
#' \code{\link{simulate.drmTMB}}), and fits the model again to each simulated response, with the iteration limits of the fit (\code{control} in \code{\link{drmTMB}}), starting at the
#' estimates. The bootstrap sample can be used for confidence intervals of the model parameters and of
#' effective doses, see \code{\link{confint.drmTMB}}, without fitting the models again.
#'
#' @param object a drmTMB object
#' @param nsim number of simulations
#' @param seed an optional seed for the random number generator
#' @param progress logical, show a progress bar
#' @param x a bootstrap.drmTMB object
#' @param digits number of significant digits in the printed output
#' @param parm parameters, names or positions in \code{coef}; default all
#' @param level the confidence level
#' @param ... ignored
#'
#' @return bootstrap returns a bootstrap.drmTMB object, a list with the estimates of \code{coef(object)}
#'   for each simulation (\code{estimates}, a matrix with NA rows for fits that did not converge), the
#'   fitted model (\code{fit}), and \code{nsim} and \code{seed}. confint returns percentile intervals
#'   of the parameters, on the scale of \code{\link{confint.drmTMB}}.
#' @export
#' @aliases bootstrap print.bootstrap.drmTMB confint.bootstrap.drmTMB
bootstrap <- function(object, nsim=500, seed=NULL, progress=interactive()){
  if (!inherits(object, "drmTMB")) stop("object must be a drmTMB object.")
  sims <- simulate(object, nsim=nsim, seed=seed)
  env <- object$obj$env
  rnd <- !is.null(env$random)
  pars <- env$parList(object$opt$par)
  # iteration limits of the fit
  control <- if (is.null(object$control)) drmTMBcontrol() else object$control
  ctrl <- list(iter.max=control$iter.max, eval.max=control$eval.max)
  pb <- progressbar(nsim, progress)
  est <- t(sapply(sims, function(y){
    on.exit(pb$tick())
    dat <- env$data
    dat$y <- y
    obj <- MakeADFun(dat, pars, map=env$map, random=if (rnd) "u", DLL=env$DLL, silent=TRUE)
    opt <- suppressWarnings(try(nlminb(obj$par, obj$fn, obj$gr, if (!rnd) obj$he,
                                       lower=object$lower, upper=object$upper, control=ctrl), silent=TRUE))
    if (inherits(opt, "try-error") || opt$convergence != 0 || !is.finite(opt$objective)) return(rep(NA, length(object$opt$par)))
    opt$par
  }))
  pb$close()
  est <- matrix(est, nrow=nsim)
  colnames(est) <- names(coef(object))
  nfail <- sum(is.na(est[, 1]))
  if (nfail > 0) message(nfail, " of ", nsim, " bootstrap fits did not converge and are left out.")
  res <- list(estimates=est, fit=object, nsim=nsim, seed=seed)
  class(res) <- "bootstrap.drmTMB"
  res
}

#' @rdname bootstrap
#' @export
print.bootstrap.drmTMB <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  est <- x$estimates
  ok <- !is.na(est[, 1])
  cat("Parametric bootstrap: ", sum(ok), " of ", x$nsim, " fits converged\n\n", sep="")
  orig <- coef(x$fit)
  tab <- cbind(Estimate=orig, Mean=colMeans(est[ok, , drop=FALSE]),
               Bias=colMeans(est[ok, , drop=FALSE]) - orig, "Std. Error"=apply(est[ok, , drop=FALSE], 2, sd))
  print(tab, digits=digits)
  invisible(x)
}

#' @rdname bootstrap
#' @export
confint.bootstrap.drmTMB <- function(object, parm=NULL, level=0.95, ...){
  est <- coef(object$fit)
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown parameters in parm.")
  ci <- t(apply(bootconverged(object)[, parm, drop=FALSE], 2, quantile, probs=c((1 - level)/2, (1 + level)/2), na.rm=TRUE, names=FALSE))
  natscale(ci, names(est)[parm], corr2names(object$fit), level)
}

# bootstrap sample for confint: the given one, checked to belong to the fit, or a new one
getboot <- function(fit, boot, nsim, seed, progress=interactive(), ...){
  if (is.null(boot)){
    message("Computing bootstrap confidence intervals (", nsim, " simulations) ...")
    return(bootstrap(fit, nsim=nsim, seed=seed, progress=progress))
  }
  if (!inherits(boot, "bootstrap.drmTMB")) stop("boot must be the result of bootstrap().")
  if (!identical(coef(boot$fit), coef(fit)) || !identical(boot$fit$call, fit$call)) stop("boot is a bootstrap sample of a different model.")
  boot
}

# estimates of the converged bootstrap fits, at least two for standard errors and intervals
bootconverged <- function(boot){
  est <- boot$estimates[!is.na(boot$estimates[, 1]), , drop=FALSE]
  if (nrow(est) < 2) stop("only ", nrow(est), " of ", boot$nsim, " bootstrap fits converged; at least 2 are needed for bootstrap intervals.")
  est
}

# text progress bar of package utils, or nothing
progressbar <- function(n, show){
  if (!show || n < 1) return(list(tick=function() invisible(), close=function() invisible()))
  pb <- txtProgressBar(min=0, max=n, style=3)
  i <- 0
  list(tick=function(){
         i <<- i + 1
         setTxtProgressBar(pb, i)
       },
       close=function() close(pb))
}
