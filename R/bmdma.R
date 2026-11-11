#' Model-averaged benchmark dose
#'
#' Benchmark doses (see \code{\link{bmd}}) averaged over several dose-response models fitted to the
#' same data, with AIC or BIC weights, and bootstrap confidence intervals: the BMDL is the lower limit of the
#' two-sided 90\% interval, \code{confint(bmdma(...))}.
#'
#' @param ... drmTMB objects, fitted to the same response, doses and weights; the models and families
#'   may differ (e.g. gaussian and lognormal, whose AICs are comparable). Arguments are named by the
#'   argument names, or by the expressions.
#' @param bmr,def,newdata,scale,nodes the benchmark response and the definition, the covariates of the
#'   curves, and the BMD of the marginal or conditional curve with random effects, as in
#'   \code{\link{bmd}}; newdata cannot give fitted groups of random effects, and \code{def="sd"} with
#'   \code{average="curve"} needs gaussian families
#' @param average "curve" (default) for the BMD of the model-averaged curve, the weighted mean of the
#'   curves of the models, or "bmd" for the weighted mean of the BMDs of the models
#' @param weights the information criterion of the model weights, "AIC" (default) or "BIC"
#' @param nsim number of bootstrap samples; 0 for no bootstrap
#' @param seed an optional seed for the random number generator
#' @param progress logical, show a progress bar
#' @param x,object a bmdma object
#' @param digits number of significant digits in the printed output
#' @param parm BMDs, names or positions in \code{coef}; default all
#' @param level the confidence level; 0.90 by default, for the BMDL and BMDU
#' @param method "boot" for bootstrap percentile intervals, or "buckland" for Wald intervals with the
#'   standard error of Buckland et al. (1997), for \code{average="bmd"} only
#'
#' @details Model k gets the weight exp(-IC_k/2), relative to the sum over the models, with the
#' information criterion IC_k, AIC_k = -2 logLik_k + 2 df_k or BIC_k = -2 logLik_k + log(n) df_k, with
#' the number of observations n (\code{nobs}), as in Ritz, Gerhard and Hothorn (2013). With
#' \code{average="curve"}, the BMD is that of the averaged curve, the weighted mean of the curves (the
#' marginal mean curves with \code{scale="marginal"} and random effects), with the benchmark response
#' relative to the averaged response at dose zero, as in EFSA (2017); with \code{def="sd"}, the standard
#' deviation is the weighted mean of the residual standard deviations. With \code{average="bmd"}, it is
#' the weighted mean of the BMDs of the models; models whose curve does not reach the benchmark response
#' are left out, with the weights of the others rescaled.
#'
#' The bootstrap samples the responses from the mixture of the models: for each sample, a model is drawn
#' with its weight and a response simulated from it, with new random effects (see
#' \code{\link{simulate.drmTMB}}). All models are fitted to each sample, starting at their estimates,
#' and the weights and the averaged BMD are computed again. Models whose fit fails are left out of a
#' sample, with the weights of the others rescaled; samples without any fit are left out. Models of the
#' lognormal family get the weight 0 in samples with non-positive responses (e.g. simulated from a
#' gaussian model), whose likelihood is zero under the lognormal family. The
#' confidence intervals are percentile intervals of the averaged BMDs.
#'
#' @references Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model selection: an integral
#'   part of inference. Biometrics 53, 603-618.
#'
#'   Ritz, C., Gerhard, D. and Hothorn, L. A. (2013). A unified framework for benchmark dose estimation
#'   applied to mixed models and model averaging. Statistics in Biopharmaceutical Research 5, 79-90.
#'
#'   EFSA Scientific Committee (2017). Update: use of the benchmark dose approach in risk assessment.
#'   EFSA Journal 15(1), 4658.
#'
#' @return bmdma returns a bmdma object, a list with the averaged BMDs (\code{estimates}), the table of
#'   the models (\code{models}: model, family, logLik, df, AIC or BIC, weight), the \code{\link{bmd}} results
#'   of the models (\code{bmds}), the averaged BMDs of the bootstrap samples (\code{boot}, with NA rows
#'   for samples left out), their weights (\code{bootweights}, NA for failed fits) and the models the
#'   samples were simulated from (\code{source}), and the arguments. coef returns the averaged BMDs and
#'   confint their confidence intervals.
#' @export
#' @aliases bmdma print.bmdma coef.bmdma confint.bmdma
bmdma <- function(..., bmr=0.1, def=c("excess", "additional", "relative", "sd"), newdata=NULL,
                  scale=c("marginal", "conditional"), average=c("curve", "bmd"), weights=c("AIC", "BIC"), nodes=NULL,
                  nsim=500, seed=NULL, progress=interactive()){
  def <- match.arg(def)
  scale <- match.arg(scale)
  average <- match.arg(average)
  weights <- match.arg(weights)
  fits <- list(...)
  nms <- names(fits)
  exprs <- vapply(as.list(substitute(list(...)))[-1], function(e) paste(deparse(e), collapse=""), character(1))
  if (is.null(nms)) nms <- exprs
  nms[nms == ""] <- exprs[nms == ""]
  names(fits) <- make.unique(nms)
  K <- length(fits)
  if (K < 2) stop("give at least two fitted models.")
  if (!all(vapply(fits, inherits, logical(1), "drmTMB"))) stop("all models must be drmTMB objects.")
  dat <- function(f) lapply(f$obj$env$data[c("x", "y", "w", "bn")], function(v) unname(as.vector(v)))
  for (k in 2:K) if (!isTRUE(all.equal(dat(fits[[1]]), dat(fits[[k]])))){
    stop("the models must be fitted to the same doses, responses and weights: ", names(fits)[k], " differs from ", names(fits)[1], ".")
  }
  nconv <- vapply(fits, function(f) f$opt$convergence != 0, logical(1))
  if (any(nconv)) warning("the fits of ", paste(names(fits)[nconv], collapse=", "), " did not converge; their AICs and BMDs may be unreliable.", call.=FALSE)
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim < 0) stop("nsim must be a non-negative number.")
  if (def == "sd" && average == "curve" && any(vapply(fits, function(f) f$family != "gaussian", logical(1)))){
    stop("def=\"sd\" with average=\"curve\" needs gaussian families: the standard deviations of the families differ; use average=\"bmd\".")
  }

  # BMDs of each model
  bmds <- lapply(names(fits), function(k){
    tryCatch(bmd(fits[[k]], bmr=bmr, def=def, newdata=newdata, scale=scale, nodes=nodes),
             error=function(e) stop("model ", k, ": ", conditionMessage(e), call.=FALSE))
  })
  names(bmds) <- names(fits)
  if (any(vapply(bmds, function(b) any(b$conditional), logical(1)))){
    stop("BMDs of fitted groups of random effects are not available for model averaging; give newdata without fitted groups.")
  }
  edn <- names(bmds[[1]]$estimates)
  nbmd <- length(edn)
  inc <- vapply(bmds, function(b) vapply(b$bmd, function(i) i$inc, numeric(1)), numeric(nbmd))
  inc <- matrix(inc, nbmd)
  if (average == "curve" && any(apply(inc, 1, function(s) length(unique(s)) > 1))){
    stop("the curves of the models change in different directions; the averaged curve has no BMD.")
  }

  # AIC or BIC weights
  ll <- vapply(fits, function(f) -f$opt$objective, numeric(1))
  df <- vapply(fits, function(f) length(f$opt$par), numeric(1))
  pen <- if (weights == "AIC") 2 else log(nobs(fits[[1]]))
  ic <- -2*ll + pen*df
  w <- maweights(ic)
  models <- data.frame(model=vapply(fits, function(f) f$model, character(1)),
                       family=vapply(fits, function(f) f$family, character(1)),
                       logLik=ll, df=df, IC=ic, weight=w, row.names=names(fits))
  names(models)[5] <- weights

  # averaged BMDs, from the BMDs of the models, or by the root of the averaged curve starting at the
  # averaged BMDs
  est <- matrix(vapply(bmds, function(b) unname(b$estimates), numeric(nbmd)), nbmd)
  bavg <- apply(est, 1, function(d) bmdmean(d, w))
  dpos <- mean(fits[[1]]$obj$env$data$x[fits[[1]]$obj$env$data$x > 0])
  x0 <- ifelse(is.na(bavg), dpos, bavg)
  if (average == "bmd"){
    ma <- bavg
  } else {
    cv <- lapply(bmds, function(b) lapply(seq_len(nbmd), function(i) bmdcurve(b, i, b$fit$opt$par)))
    ma <- vapply(seq_len(nbmd), function(i) macurve(lapply(cv, `[[`, i), w, def, bmds[[1]]$q[[i]], x0[i]), numeric(1))
  }
  names(ma) <- edn
  if (anyNA(ma)) warning("the averaged curve does not reach the benchmark response for ", paste(edn[is.na(ma)], collapse=", "), "; the BMDs are NA.", call.=FALSE)
  x0 <- ifelse(is.na(ma), x0, ma)

  # bootstrap: responses from the mixture of the models, all models fitted again
  boot <- matrix(NA_real_, nsim, nbmd, dimnames=list(NULL, edn))
  bw <- matrix(NA_real_, nsim, K, dimnames=list(NULL, names(fits)))
  src <- integer(nsim)
  fail <- outs <- matrix(FALSE, nsim, K)
  if (nsim > 0){
    if (!is.null(seed)) set.seed(seed)
    pars <- lapply(fits, function(f) f$obj$env$parList(f$opt$par))
    ctrls <- lapply(fits, function(f){
      control <- if (is.null(f$control)) drmTMBcontrol() else f$control
      list(iter.max=control$iter.max, eval.max=control$eval.max)
    })
    pb <- progressbar(nsim, progress)
    for (r in seq_len(nsim)){
      src[r] <- sample.int(K, 1, prob=w)
      y <- suppressWarnings(simulate(fits[[src[r]]], nsim=1)[[1]])
      pb$tick()
      if (any(!is.finite(y))) next
      # responses outside the support of a model (non-positive for the lognormal family): weight 0
      supp <- vapply(fits, function(f) f$family != "lognormal" || all(y[f$obj$env$data$w > 0] > 0), logical(1))
      opts <- lapply(seq_len(K), function(k){
        if (!supp[k]) return(NULL)
        opt <- refitopt(fits[[k]], y, pars[[k]], ctrls[[k]])
        if (is.null(opt) || opt$convergence != 0 || !is.finite(opt$objective)) NULL else opt
      })
      ok <- !vapply(opts, is.null, logical(1))
      outs[r, ] <- !supp
      fail[r, ] <- supp & !ok
      if (!any(ok)) next
      icr <- vapply(opts[ok], function(o) 2*o$objective + pen*length(o$par), numeric(1))
      bw[r, !supp] <- 0
      bw[r, ok] <- wr <- maweights(icr)
      boot[r, ] <- vapply(seq_len(nbmd), function(i){
        cvr <- lapply(which(ok), function(k) bmdcurve(bmds[[k]], i, opts[[k]]$par))
        if (average == "bmd") bmdmean(vapply(cvr, function(cr) cr$bmd(), numeric(1)), wr)
        else macurve(cvr, wr, def, bmds[[1]]$q[[i]], x0[i])
      }, numeric(1))
    }
    pb$close()
    nout <- colSums(outs)
    if (any(nout > 0)) message("bootstrap samples outside the support of the model, with weight 0: ", paste(names(fits)[nout > 0], nout[nout > 0], collapse=", "), " of ", nsim, ".")
    nfail <- colSums(fail)
    if (any(nfail > 0)) message("bootstrap fits that failed and were left out of their samples: ", paste(names(fits)[nfail > 0], nfail[nfail > 0], collapse=", "), " of ", nsim, ".")
    nna <- sum(is.na(boot[, 1]))
    if (nna > 0) message(nna, " of ", nsim, " bootstrap samples gave no averaged BMD and are left out.")
  }

  res <- list(estimates=ma, models=models, bmds=bmds, boot=boot, bootweights=bw, source=src,
              bmr=bmds[[1]]$q, def=def, scale=scale, average=average, weights=weights, nsim=nsim, seed=seed)
  class(res) <- "bmdma"
  res
}

#' @rdname bmdma
#' @export
print.bmdma <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("Model-averaged benchmark doses, ", c(excess="excess risk", additional="additional risk", relative="relative deviation", sd="standard deviations")[[x$def]],
      ", averaged ", if (x$average == "curve") "curves" else "BMDs", if (any(vapply(x$bmds, function(b) length(b$fit$reterms) > 0, logical(1)))) paste(",", x$scale, "BMD"), ":\n\n", sep="")
  tab <- x$models
  est <- matrix(vapply(x$bmds, function(b) unname(b$estimates), numeric(length(x$estimates))), length(x$estimates))
  bm <- t(est)
  colnames(bm) <- if (ncol(bm) == 1) "BMD" else paste("BMD", names(x$estimates))
  print(cbind(tab, bm), digits=digits)
  cat("\n")
  out <- data.frame(BMD=x$estimates, row.names=names(x$estimates))
  ci <- tryCatch(confint(x), error=function(e) NULL)
  if (!is.null(ci)){
    out$BMDL <- ci[, 1]
    out$BMDU <- ci[, 2]
  }
  print(out, digits=digits)
  if (!is.null(ci)) cat("\n", if (x$nsim > 0) paste("Bootstrap percentile intervals (", sum(!is.na(x$boot[, 1])), " samples)", sep="") else "Buckland intervals", ", level 90%\n", sep="")
  invisible(x)
}

#' @rdname bmdma
#' @export
coef.bmdma <- function(object, ...){
  object$estimates
}

#' @rdname bmdma
#' @export
confint.bmdma <- function(object, parm=NULL, level=0.90, method=c("boot", "buckland"), ...){
  if (missing(method)) method <- if (object$nsim > 0) "boot" else "buckland"
  method <- match.arg(method)
  est <- object$estimates
  if (is.null(parm)) parm <- seq_along(est)
  if (is.character(parm)) parm <- match(parm, names(est))
  if (any(is.na(parm))) stop("unknown BMDs in parm.")
  probs <- c((1 - level)/2, (1 + level)/2)
  if (method == "boot"){
    if (object$nsim == 0) stop("no bootstrap samples: use bmdma with nsim > 0, or method=\"buckland\".")
    b <- object$boot[!is.na(object$boot[, 1]), , drop=FALSE]
    if (nrow(b) < 2) stop("fewer than 2 bootstrap samples gave an averaged BMD.")
    ci <- t(apply(b[, parm, drop=FALSE], 2, quantile, probs=probs, na.rm=TRUE, names=FALSE))
  } else {
    if (object$average != "bmd") stop("Buckland intervals are only available for average=\"bmd\".")
    w <- object$models$weight
    se <- vapply(parm, function(i){
      d <- vapply(object$bmds, function(b) unname(b$estimates[i]), numeric(1))
      s <- vapply(object$bmds, function(b) sqrt(b$variance[i, i]), numeric(1))
      ok <- !is.na(d)
      sum(w[ok] * sqrt(s[ok]^2 + (d[ok] - est[i])^2)) / sum(w[ok])
    }, numeric(1))
    z <- qnorm(probs[2])
    ci <- cbind(est[parm] - z*se, est[parm] + z*se)
  }
  dimnames(ci) <- list(names(est)[parm], paste(format(100*probs, trim=TRUE, scientific=FALSE, digits=3), "%"))
  ci
}

# weights of the models with information criteria ic
maweights <- function(ic){
  w <- exp(-(ic - min(ic))/2)
  w / sum(w)
}

# weighted mean of the BMDs d with weights w, leaving out NA BMDs
bmdmean <- function(d, w){
  ok <- !is.na(d)
  if (!any(ok)) NA else sum(w[ok]*d[ok]) / sum(w[ok])
}

# curve of BMD i of the bmd result b at the fixed parameters par (as opt$par of the fit): the mean
# (marginal mean if averaged over random effects) as a function of the dose, the BMD of the curve as a
# function, the residual standard deviation for def="sd", and the direction
bmdcurve <- function(b, i, par){
  fit <- b$fit
  info <- b$bmd[[i]]
  pv <- fit$start
  pv[!fit$fix] <- par[seq_len(sum(!fit$fix))]
  lpv <- as.vector(b$linfct[[i]] %*% pv)
  bp <- bmdpars(fit, par)
  lnk <- linkindex(fit$link)
  if (length(info$wq) > 0){
    B <- linkpars(lnk, t(lpv + info$A %*% recholesky(fit, bp$theta) %*% info$E))
    wq <- info$wq
  } else {
    B <- matrix(linkpars(lnk, lpv), 1)
    wq <- 1
  }
  mod <- modelindex(fit$model)
  sm <- smoothwidth(fit$smooth)
  list(mean=function(d) sum(wq * drcmeanR(mod, d, B[, 1], B[, 2], B[, 3], B[, 4], B[, 5], fit$alpha, sm)),
       bmd=function() bmdvalue(fit, info, lpv, bp$theta, bp$ls),
       sigma=info$sfac*exp(bp$ls), inc=info$inc)
}

# BMD of the curve averaged with weights w over the curves cv (from bmdcurve), with the definition def
# and the benchmark response bmr, starting the root search at x0
macurve <- function(cv, w, def, bmr, x0){
  inc <- cv[[1]]$inc
  m <- function(d) sum(w * vapply(cv, function(cr) cr$mean(d), numeric(1)))
  sigma <- sum(w * vapply(cv, function(cr) cr$sigma, numeric(1)))
  target <- bmdtarget(def, bmr, inc, m(0), sigma)
  bmdroot(function(d) inc*(m(d) - target), x0)
}
