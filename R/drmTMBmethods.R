#' drmTMB Methods
#'
#' @param x,object a drmTMB object
#' @param newdata a data frame for prediction. If NULL, the original data is used.
#' @param random logical, if TRUE, the random effects are included in the predictions. For simulate, if TRUE, new random effects are drawn for each simulation, otherwise the fitted random effects are kept.
#' @param digits number of significant digits in the printed output
#' @param nsim number of response vectors to simulate
#' @param seed an optional seed for the random number generator
#' @param ... anything else
#'
#' @return Estimated model coefficients, variance-covariance matrix, or predictions.
#' simulate returns a data frame with nsim columns of simulated responses at the estimated
#' parameters; for a binomial response given as cbind(successes, failures), the number of successes.
#' @export
#'
#' @name drmTMB-methods
#' @aliases print.drmTMB predict.drmTMB coef.drmTMB vcov.drmTMB simulate.drmTMB
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

  # coefficients fixed by the user; the parameters the model does not use are always fixed
  mod <- which(x$model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
  npi <- as.integer(sub("^b([1-5]):.*", "\\1", x$coefnames))
  ufix <- x$fix & !(npi %in% unusedpars(mod))
  if (any(ufix)){
    cat("Fixed: ", paste(x$coefnames[ufix], "=", format(x$start[ufix], digits=digits), collapse=", "), "\n", sep="")
  }

  # random effects: one row per parameter of each term, with standard deviation and correlations
  vcs <- recov(x)
  if (length(vcs) > 0){
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

  # residual standard deviation (gaussian) or precision (beta)
  lsig <- x$sdrl$par.fixed["log_sigma"]
  if (x$family == "gaussian") cat("\nResidual standard deviation: ", format(exp(lsig), digits=digits), "\n", sep="")
  if (x$family == "beta") cat("\nPrecision: ", format(exp(lsig), digits=digits), "\n", sep="")
  invisible(x)
}

#' @rdname drmTMB-methods
predict.drmTMB <- function(object, newdata=NULL, random=TRUE){
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
    names(est)[seq_len(sum(!object$fix))] <- object$coefnames[!object$fix]
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
