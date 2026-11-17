#' Bayesian model-averaged benchmark dose
#'
#' The posterior distribution of the benchmark dose averaged over dose-response models fitted in the
#' BMD parameterisation with priors, as in the Bayesian benchmark dose modelling of EFSA (2022): the
#' models are weighted by their posterior probabilities, from their marginal likelihoods
#' (\code{\link{marglik}}), and the BMD and its credible interval are the median and quantiles of the
#' mixture of the posterior distributions of the BMD of the models.
#'
#' @param ... drmTMB objects fitted in the BMD parameterisation (\code{bmdpar} of \code{\link{drmTMB}})
#'   with the same benchmark response and definition, the same prior for the BMD
#'   (\code{prior=list(bmd=...)}) and proper priors for all estimated parameters, to the same doses,
#'   responses and weights; the models and families may differ. Arguments are named by the argument
#'   names, or by the expressions.
#' @param prior.weights prior probabilities of the models, in the order of \code{...}; equal by default.
#' @param method "laplace" (default) for the Laplace approximation of the marginal likelihoods and a
#'   normal posterior distribution of the log BMD of each model, or "is" for importance sampling of both
#'   (see \code{\link{marglik}})
#' @param n,df number of draws of the final round and degrees of freedom of the proposal of importance
#'   sampling (see \code{\link{marglik}})
#' @param x,object a bmdbayes object
#' @param digits number of significant digits in the printed output
#' @param level the probability of the credible interval; 0.90 by default, for the BMDL and BMDU
#' @param parm not used; there is one BMD
#'
#' @details The posterior probability of model k is proportional to its prior probability times its
#' marginal likelihood p(y | model k). With \code{method="laplace"}, the posterior distribution of the
#' log BMD of each model is normal with the posterior mode of log_bmd as the mean and its standard error
#' (the curvature of the log posterior) as the standard deviation, so that the median BMD of a model is
#' its posterior mode; the model-averaged posterior distribution of the log BMD is the mixture of these
#' normal distributions with the posterior probabilities of the models as weights. With
#' \code{method="is"}, the importance sampling draws of each model, weighted by their importance weights
#' and the posterior probability of the model, give both the marginal likelihoods and the mixture,
#' without the normal approximation. The BMD is the median of the mixture, and the BMDL and BMDU the
#' quantiles (1 - level)/2 and (1 + level)/2, the 5\% and 95\% quantiles of a 90\% credible interval.
#'
#' With random effects, the BMD of the BMD parameterisation is that of the curve of the intercepts, with
#' the random effects at zero (the conditional BMD of \code{\link{bmd}}).
#'
#' A model whose marginal likelihood cannot be computed (e.g. the fit did not converge, or its Hessian is
#' not positive definite) gets the weight 0, with a warning. \code{print} notes the criteria of EFSA
#' (2022) for a credible interval that is too wide, with the 90\% interval: BMD/BMDL above 20, BMDU/BMDL
#' above 50, and a BMDL below a tenth of the lowest dose above zero.
#'
#' @references EFSA Scientific Committee (2022). Guidance on the use of the benchmark dose approach in
#'   risk assessment. EFSA Journal 20(10), 7584.
#'
#' @return bmdbayes returns a bmdbayes object, a list with the model-averaged BMD (\code{estimates}), the
#'   table of the models (\code{models}: model, family, log marginal likelihood, prior and posterior
#'   probabilities, and with importance sampling the effective sample size), the posterior distribution of
#'   the log BMD of each model (\code{posterior}: the mean and standard deviation of the normal
#'   distribution, or the draws and their weights), the lowest dose above zero, and the arguments. coef
#'   returns the BMD, confint the credible interval, and \code{bmdbayesmodels} the BMDs and credible
#'   intervals of the models.
#' @export
#' @aliases bmdbayes print.bmdbayes coef.bmdbayes confint.bmdbayes
bmdbayes <- function(..., prior.weights=NULL, method=c("laplace", "is"), n=5000, df=4){
  method <- match.arg(method)
  fits <- namefits(list(...), as.list(substitute(list(...)))[-1])
  K <- length(fits)
  if (K < 1) stop("give at least one fitted model.")
  samedata(fits)
  if (is.null(prior.weights)) prior.weights <- rep(1, K)
  if (!is.numeric(prior.weights) || length(prior.weights) != K || any(!is.finite(prior.weights)) || any(prior.weights < 0) || sum(prior.weights) <= 0){
    stop("prior.weights must be ", K, " non-negative numbers, not all 0.")
  }
  prior.weights <- prior.weights/sum(prior.weights)

  # the same BMD (benchmark response, definition) and the same prior for it in all models
  for (k in names(fits)){
    if (is.null(fits[[k]]$bmdpar)) stop("model ", k, " is not fitted in the BMD parameterisation; use drmTMB with bmdpar.")
  }
  bp <- lapply(fits, function(f) f$bmdpar[c("bmr", "def")])
  for (k in seq_len(K)[-1]) if (!identical(bp[[k]], bp[[1]])){
    stop("the models must have the same benchmark response and definition (bmdpar): ", names(fits)[k], " differs from ", names(fits)[1], ".")
  }
  bprior <- vapply(fits, function(f){
    p <- f$prior
    if (is.null(p) || !("bmd" %in% p$parameter)) NA_character_ else p$prior[p$parameter == "bmd"]
  }, character(1))
  if (anyNA(bprior)) stop("the models need a prior for the BMD, prior=list(bmd=...): missing for ", paste(names(fits)[is.na(bprior)], collapse=", "), ".")
  if (length(unique(bprior)) > 1){
    stop("the models must have the same prior for the BMD, otherwise their posterior probabilities compare the priors: ",
         paste(names(fits), bprior, sep=": ", collapse=", "), ".")
  }
  for (k in names(fits)){
    tryCatch(mlcheck(fits[[k]]), error=function(e) stop("model ", k, ": ", conditionMessage(e), call.=FALSE))
  }

  # marginal likelihoods and posterior distributions of the log BMD; models that fail get weight 0
  post <- lapply(names(fits), function(k){
    f <- fits[[k]]
    i <- which(parnames(f) == "log_bmd")
    tryCatch({
      if (f$opt$convergence != 0) stop("the fit did not converge")
      if (method == "laplace"){
        la <- withCallingHandlers(mllaplace(f), warning=function(wn){
          warning("model ", k, ": ", conditionMessage(wn), call.=FALSE)
          invokeRestart("muffleWarning")
        })
        list(logml=la, mean=f$opt$par[[i]], sd=sqrt(f$sdrl$cov.fixed[i, i]))
      } else {
        s <- mlis(f, n, df)
        list(logml=s$logml, ess=s$ess, se=s$se, draws=unname(s$draws[, i]), weights=s$weights)
      }
    }, error=function(e){
      warning("model ", k, " gets the weight 0: ", conditionMessage(e), ".", call.=FALSE)
      NULL
    })
  })
  names(post) <- names(fits)
  ok <- !vapply(post, is.null, logical(1))
  if (!any(ok & prior.weights > 0)) stop("none of the models with a positive prior probability gives a marginal likelihood.")
  logml <- vapply(post, function(p) if (is.null(p)) NA_real_ else p$logml, numeric(1))
  lw <- log(prior.weights) + logml
  lw[!ok] <- -Inf
  w <- exp(lw - max(lw))
  w <- w/sum(w)

  models <- data.frame(model=vapply(fits, function(f) f$model, character(1)),
                       family=vapply(fits, function(f) f$family, character(1)),
                       logML=logml, prior=prior.weights, weight=w, row.names=names(fits))
  if (method == "is") models$ess <- vapply(post, function(p) if (is.null(p)) NA_real_ else p$ess, numeric(1))
  x <- fits[[1]]$obj$env$data$x
  res <- list(posterior=post, models=models, method=method, bmr=bp[[1]]$bmr, def=bp[[1]]$def, prior=bprior[[1]],
              lowdose=min(x[x > 0]))
  class(res) <- "bmdbayes"
  res$estimates <- c(BMD=exp(bmdbayesq(res, 0.5)))
  res
}

# quantiles p of the posterior distribution of the log BMD: of the mixture over the models with
# weights w, or of model k
bmdbayesq <- function(object, p, k=NULL){
  post <- object$posterior
  w <- object$models$weight
  if (!is.null(k)){
    w <- as.numeric(seq_along(post) == k)
    if (is.null(post[[k]])) return(rep(NA_real_, length(p)))
  }
  use <- which(w > 0 & !vapply(post, is.null, logical(1)))
  if (object$method == "laplace"){
    mu <- vapply(post[use], function(q) q$mean, numeric(1))
    s <- vapply(post[use], function(q) q$sd, numeric(1))
    ww <- w[use]
    cdf <- function(z) sum(ww*pnorm(z, mu, s))
    lo <- min(mu - 10*s)
    hi <- max(mu + 10*s)
    vapply(p, function(pp) uniroot(function(z) cdf(z) - pp, c(lo, hi), tol=1e-10)$root, numeric(1))
  } else {
    d <- unlist(lapply(post[use], function(q) q$draws))
    dw <- unlist(lapply(seq_along(use), function(j) w[use[j]]*post[[use[j]]]$weights))
    keep <- dw > 0
    wquantile(d[keep], dw[keep], p)
  }
}

# weighted quantiles p of x with weights w: the smallest x whose cumulative weight reaches p
wquantile <- function(x, w, p){
  o <- order(x)
  x <- x[o]
  cw <- cumsum(w[o])/sum(w)
  vapply(p, function(pp) x[min(which(cw >= pp - 1e-12))], numeric(1))
}

#' @rdname bmdbayes
#' @export
coef.bmdbayes <- function(object, ...){
  object$estimates
}

#' @rdname bmdbayes
#' @export
confint.bmdbayes <- function(object, parm=NULL, level=0.90, ...){
  probs <- c((1 - level)/2, (1 + level)/2)
  ci <- matrix(exp(bmdbayesq(object, probs)), 1)
  dimnames(ci) <- list("BMD", paste(format(100*probs, trim=TRUE, scientific=FALSE, digits=3), "%"))
  ci
}

#' @rdname bmdbayes
#' @export
bmdbayesmodels <- function(object, level=0.90){
  if (!inherits(object, "bmdbayes")) stop("object must be a bmdbayes object.")
  probs <- c(0.5, (1 - level)/2, (1 + level)/2)
  q <- t(vapply(seq_along(object$posterior), function(k) exp(bmdbayesq(object, probs, k)), numeric(3)))
  data.frame(BMD=q[, 1], BMDL=q[, 2], BMDU=q[, 3], row.names=rownames(object$models))
}

#' @rdname bmdbayes
#' @export
print.bmdbayes <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("Bayesian model-averaged benchmark dose, ", deflabel(x$def),
      " ", format(x$bmr), ", prior ", x$prior, " for the BMD\n",
      if (x$method == "laplace") "Laplace approximation" else "Importance sampling", "\n\n", sep="")
  print(cbind(x$models, bmdbayesmodels(x)), digits=digits)
  ci <- confint(x)
  out <- data.frame(BMD=x$estimates, BMDL=ci[, 1], BMDU=ci[, 2], row.names="model average")
  cat("\n")
  print(out, digits=digits)
  cat("\nPosterior median and 90% credible interval\n")
  notes <- c(if (out$BMD/out$BMDL > 20) "BMD/BMDL above 20",
             if (out$BMDU/out$BMDL > 50) "BMDU/BMDL above 50",
             if (out$BMDL < x$lowdose/10) "BMDL below a tenth of the lowest dose")
  if (length(notes) > 0) cat("Note (EFSA 2022): ", paste(notes, collapse=", "), "; the credible interval is wide.\n", sep="")
  invisible(x)
}
