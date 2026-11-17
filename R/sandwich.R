#' Robust (sandwich) covariance of a drmTMB fit
#'
#' vcovrobust computes the sandwich covariance V M V of all estimated parameters (the coefficients
#' and the variance parameters, as \code{vcov}), with the inverse Hessian V of the fit and the outer
#' products M of the scores of independent clusters, the gradients of their contributions to the
#' log-likelihood. It is robust to a misspecified variance (e.g. heteroscedasticity, overdispersion)
#' and, with clusters, to correlation within the clusters. The result can be given as the argument
#' vcov. of \code{\link{summary.drmTMB}}, \code{\link{confint.drmTMB}} (Wald intervals),
#' \code{\link{predict.drmTMB}}, \code{\link{ed}}, \code{\link{peak}} and \code{\link{bmd}}, e.g.
#' \code{ed(fit, 50, vcov.=vcovrobust)} or \code{ed(fit, 50, vcov.=function(f) vcovrobust(f, cluster=~plate))}.
#'
#' Without random effects, the scores of the observations are the gradients of their contributions to
#' the log-likelihood, by automatic differentiation. With random effects, the clusters must contain
#' whole groups: the score of a cluster is the gradient of the Laplace approximation of its marginal
#' log-likelihood, from a fit of the model with the weights of all other observations set to zero (as
#' in glmmTMB). The groups of all random effect terms, the series of a correlation structure
#' (\code{\link{correlation}}) and the series of event data must be nested in the clusters; by default,
#' the clusters are the groups of the random effect term that all others are nested in (e.g. lab for
#' 1 | lab/plate), or the series, and otherwise the observations.
#'
#' The sandwich covariance is not available for fits with priors, whose prior densities are not
#' contributions of independent clusters. Small numbers of clusters give too small standard errors;
#' type="HC1" multiplies the covariance by G/(G - 1) (n - 1)/(n - p), for G clusters, n observations with
#' positive weights and p estimated parameters, as type="HC1" of \code{sandwich::vcovCL}.
#'
#' estfun and bread methods for the generics of package sandwich are provided, with the scores of the
#' default clusters as estfun and n V as bread for n clusters, so that \code{sandwich::sandwich(fit)}
#' equals \code{vcovrobust(fit)}, and \code{sandwich::vcovCL} can aggregate the scores of the
#' observations of a fit without random effects or correlation structures. With random effects, estfun has
#' a row for each cluster.
#'
#' @param object,x a drmTMB object
#' @param cluster the clusters: a one-sided formula, e.g. ~ plate, evaluated in the data of the fit, or a
#'   vector with an entry for each observation of the fit (the rows of the data without missing values).
#'   NULL for the default clusters, see Details.
#' @param type "HC0" (default) without, or "HC1" with the small-sample adjustment.
#' @param ... anything else
#'
#' @return vcovrobust and bread return a covariance matrix of the estimated parameters, named as
#'   \code{vcov}; estfun returns a matrix of scores with a row for each cluster (or observation) and a
#'   column for each parameter.
#'
#' @examples
#' \dontrun{
#' fit <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic")
#' vcovrobust(fit)
#' summary(fit, vcov.=vcovrobust)
#' ed(fit, c(10, 50), vcov.=vcovrobust)
#' }
#' @export
vcovrobust <- function(object, cluster=NULL, type=c("HC0", "HC1")){
  type <- match.arg(type)
  S <- scores(object, cluster)
  G <- nrow(S)
  if (G < 2) stop("the sandwich covariance needs at least 2 clusters.")
  V <- mlsdr(object)$cov.fixed
  vc <- V %*% crossprod(S) %*% V
  if (type == "HC1"){
    n <- sum(object$obj$env$data$w > 0)
    p <- ncol(S)
    if (n <= p) stop("type=\"HC1\" needs more observations than parameters.")
    vc <- vc * G/(G - 1) * (n - 1)/(n - p)
  }
  rownames(vc) <- colnames(vc) <- names(coef(object))
  vc
}

#' @rdname vcovrobust
#' @exportS3Method sandwich::estfun
estfun.drmTMB <- function(x, cluster=NULL, ...){
  S <- scores(x, cluster)
  colnames(S) <- names(coef(x))
  S
}

#' @rdname vcovrobust
#' @exportS3Method sandwich::bread
bread.drmTMB <- function(x, cluster=NULL, ...){
  V <- mlsdr(x)$cov.fixed
  rownames(V) <- colnames(V) <- names(coef(x))
  nlevels(scoreclusters(x, cluster)) * V
}

# scores (gradients of the log-likelihood) of the clusters of a fit, a row for each cluster with a
# positive weight
scores <- function(object, cluster=NULL){
  if (!is.null(object$prior)) stop("the sandwich covariance is not available for fits with priors.")
  cl <- scoreclusters(object, cluster)
  env <- object$obj$env
  dat <- env$data
  dat$scores <- 1
  if (is.null(env$random)){
    # contributions of the observations to the negative log-likelihood and their Jacobian
    sobj <- MakeADFun(dat, env$parList(object$opt$par), map=env$map, ADreport=TRUE, DLL=env$DLL, silent=TRUE)
    S <- -rowsum(sobj$gr(object$opt$par), cl, reorder=FALSE)
  } else {
    # Laplace approximation of the marginal likelihood of each cluster: the fit with the weights of the
    # other observations set to zero; the random effects of the other groups only have their normal
    # density, which integrates to one
    par <- env$last.par.best
    par[-env$random] <- object$opt$par
    par[env$random] <- object$sdrl$par.random
    pl <- env$parList(object$opt$par, par)
    w <- dat$w
    res <- lapply(levels(cl), function(g){
      dat$w <- ifelse(cl == g, w, 0)
      cobj <- MakeADFun(dat, pl, map=env$map, random="u", DLL=env$DLL, silent=TRUE)
      setinner(cobj, fitcontrol(object))
      nll <- cobj$fn(object$opt$par)
      list(nll=nll, gr=cobj$gr(object$opt$par))
    })
    total <- sum(vapply(res, function(r) r$nll, numeric(1)))
    if (!isTRUE(all.equal(total, object$opt$objective, tolerance=1e-6))){
      warning("the negative log-likelihoods of the clusters sum to ", signif(total, 8), ", not to that of the fit, ",
              signif(object$opt$objective, 8), "; check that the random effect groups are nested in the clusters.", call.=FALSE)
    }
    S <- -do.call(rbind, lapply(res, function(r) r$gr))
    rownames(S) <- levels(cl)
  }
  # clusters without observations of positive weight
  S[tapply(dat$w, cl, function(v) any(v > 0))[rownames(S)], , drop=FALSE]
}

# clusters of the observations of a fit: from cluster (formula or vector), or by default the coarsest
# of the factors that must be nested in them (random effect groups, series), or the observations; the
# factors must be nested in the clusters
scoreclusters <- function(object, cluster=NULL){
  n <- length(object$obj$env$data$y)
  nested <- nestedfactors(object)
  if (is.null(cluster)){
    if (length(nested) == 0) return(factor(seq_len(n)))
    top <- Filter(function(f) all(vapply(nested, function(g) isnested(g, f), logical(1))), nested)
    if (length(top) == 0) stop("the random effect groups (or series) are not nested in one another; give the clusters in cluster.")
    return(top[[which.min(vapply(top, nlevels, integer(1)))]])
  }
  if (inherits(cluster, "formula")){
    if (length(cluster) != 2) stop("cluster must be a one-sided formula, e.g. ~ plate, or a vector.")
    cl <- groupfactor(cluster[[2]], object$data)
  } else {
    if (length(cluster) != n) stop("cluster must have an entry for each of the ", n, " observations of the fit.")
    cl <- factor(cluster)
  }
  if (anyNA(cl)) stop("cluster has missing values.")
  cl <- droplevels(cl)
  for (nm in names(nested)){
    if (!isnested(nested[[nm]], cl)) stop("the ", nm, " must be nested in the clusters.")
  }
  cl
}

# factors of the observations that must be nested in the clusters of the scores: the groups of the
# random effect terms, the series of a correlation structure and of event data
nestedfactors <- function(object){
  dat <- object$obj$env$data
  f <- list()
  for (tm in object$reterms){
    f[[paste("groups of", tm$group)]] <- droplevels(groupfactor(str2lang(tm$group), object$data))
  }
  if (length(dat$clen) > 0 && dat$ctype > 0){
    s <- integer(length(dat$y))
    s[dat$cord + 1] <- rep(seq_along(dat$clen), dat$clen)
    f[["series of the correlation structure"]] <- factor(s)
  }
  if (identical(object$family, "event")) f[["series of the event data"]] <- droplevels(factor(object$series))
  f
}

# whether factor g is nested in factor f: each level of g within one level of f
isnested <- function(g, f) all(tapply(as.integer(f), g, function(v) length(unique(v)) == 1))

# the fit with the covariance vcov. of its estimated parameters: a matrix, or a function of the fit
# that returns one (e.g. vcovrobust), for the Wald standard errors and intervals of summary, confint,
# predict, ed and bmd; the fit unchanged for vcov. = NULL
vcovfit <- function(object, vcov.){
  if (is.null(vcov.)) return(object)
  V <- if (is.function(vcov.)) vcov.(object) else vcov.
  np <- length(object$opt$par)
  V <- tryCatch(as.matrix(V), error=function(e) NULL)
  if (!is.numeric(V) || !identical(dim(V), c(np, np)) || anyNA(V)){
    stop("vcov. must be a function of the fit or a ", np, " x ", np, " covariance matrix of all estimated parameters, as vcov() of the fit.")
  }
  if (!is.null(rownames(V)) && !identical(rownames(V), names(coef(object)))) stop("the names of vcov. are not those of vcov() of the fit.")
  object$sdrml <- mlsdr(object)
  dimnames(V) <- dimnames(object$sdrl$cov.fixed)
  object$sdrl$cov.fixed <- V
  object$vcovuser <- TRUE
  object
}

# sdreport of a fit with the inverse Hessian as cov.fixed, also for a fit with a covariance vcov.
mlsdr <- function(object) if (is.null(object$sdrml)) object$sdrl else object$sdrml

# joint covariance of all parameters (fixed and random, in the order of the parameters of the TMB
# object) for a covariance V of the fixed parameters given as vcov., from the joint precision P of
# sdreport at the Hessian of the fit: the random effects u(theta) + e with the gradient
# A = -P_rr^-1 P_rf of the predicted random effects and the prediction error e of covariance P_rr^-1
# (TMB's joint precision has the inverse of this covariance for V the inverse Hessian); V may be
# singular, e.g. for fewer clusters than parameters
jointcovuser <- function(object, P){
  r <- object$obj$env$random
  V <- object$sdrl$cov.fixed
  Prr <- P[r, r]
  A <- -as.matrix(Matrix::solve(Prr, P[r, -r]))
  C <- matrix(0, nrow(P), ncol(P))
  C[-r, -r] <- V
  C[r, -r] <- A %*% V
  C[-r, r] <- t(C[r, -r])
  C[r, r] <- as.matrix(Matrix::solve(Prr)) + A %*% V %*% t(A)
  C
}

# draws of a normal distribution with mean m (a column for each draw) and a covariance C that may be
# singular (a covariance vcov., see vcovfit)
psddraws <- function(m, C, nsim){
  e <- eigen(C, symmetric=TRUE)
  m + e$vectors %*% (sqrt(pmax(e$values, 0)) * matrix(rnorm(length(m)*nsim), nrow=length(m)))
}
