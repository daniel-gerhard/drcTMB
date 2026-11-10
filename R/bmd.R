#' Benchmark dose estimation
#'
#' The benchmark dose (BMD) is the dose at which the response has changed from the response at dose
#' zero by a specified amount, the benchmark response (BMR). Its lower confidence limit, the BMDL, is
#' the reference point of regulatory risk assessment; with the convention of EFSA and the US EPA, the
#' lower limit of a two-sided 90\% interval, \code{confint(bmd(fit), level=0.90)}.
#'
#' @param x a drmTMB object
#' @param bmr the benchmark response, see \code{def}; a vector gives a BMD for each value
#' @param def the definition of the benchmark response, see Details: "excess" (extra) risk or
#'   "additional" (added) risk for quantal responses, "relative" deviation, or "sd", the change in
#'   standard deviations of the response (gaussian and lognormal families, models without random effects)
#' @param linfct a matrix of linear functions to identify the parameters of the dose-response curve,
#'   or a list of such matrices, as in \code{\link{ed}}; each matrix is paired with an entry of bmr
#' @param newdata a data frame with the covariates of the curves, as an alternative to linfct, as in
#'   \code{\link{ed}}; rows with a level of a grouping factor of the fitted data give the BMD of the
#'   curve of that group
#' @param scale with random effects, "marginal" (default) for the BMD of the marginal (population-averaged)
#'   mean curve, the mean response over the distribution of the random effects, or "conditional" for
#'   the BMD of the typical curve, with the random effects at zero
#' @param nodes number of Gauss-Hermite quadrature nodes per random coefficient for the marginal mean;
#'   the default decreases with the number of random coefficients, from 30 for one to 5 for five or more
#'
#' @details With the response mu0 at dose zero, the BMD is the dose at which the curve reaches the
#' response
#' \itemize{
#'   \item excess (extra) risk: mu0 + bmr (1 - mu0) for an increasing curve, and mu0 (1 - bmr) for a
#'     decreasing curve, i.e. the change relative to the proportion that can still change;
#'   \item additional (added) risk: mu0 + bmr, or mu0 - bmr for a decreasing curve;
#'   \item relative deviation: mu0 (1 + bmr), or mu0 (1 - bmr) for a decreasing curve;
#'   \item sd: mu0 + bmr sigma, or mu0 - bmr sigma for a decreasing curve, with the residual standard
#'     deviation sigma (for weights of 1); for the lognormal family, mu0 exp(bmr sigma), or
#'     mu0 exp(-bmr sigma), with the standard deviation sigma of log(y), as in BMDS.
#' }
#' The definitions excess and additional are for responses that are proportions (binomial,
#' beta-binomial and beta families). For excess risk, the BMD of a curve from zero to one is the
#' effective dose \code{ed(x, 100*bmr)}. The response at dose zero, and sigma, are estimated, and their
#' uncertainty is included in the standard errors and confidence intervals. Curves that do not reach
#' the response give NA, with a warning.
#'
#' For the lognormal family, mu0 is the median at dose zero; as the mean is the median times a constant,
#' the relative BMD is the same for the median and the mean.
#'
#' With random effects, the BMD depends on the curve it refers to. The conditional BMD is that of the
#' typical curve, with all random effects at zero, as the effective doses of \code{\link{ed}} with
#' \code{random=FALSE}. The marginal BMD is that of the mean response of the population of groups,
#' the mean curve over the distribution of the random effects, found by Gauss-Hermite quadrature.
#' Both are the same if the random effects enter the mean linearly, i.e. on b1 and b2 with the
#' identity link; otherwise they differ, and the marginal BMD also depends on the variances of the
#' random effects, whose uncertainty is included in its standard errors. The standard errors are for
#' the uncertainty of the estimates, not the variation between groups. Rows of newdata with fitted
#' groups give the BMDs of the curves of these groups, with their predicted random effects; with
#' \code{scale="marginal"}, the terms without a fitted group are averaged over.
#'
#' Confidence intervals: see \code{\link{confint.edTMB}}. Profile likelihood intervals are available for
#' BMDs without averaging over random effects, except for the models logistic and gompertz and the
#' smoothed threshold of nec, whose response at dose zero depends on b4; otherwise Wald intervals are
#' returned, with a message, and bootstrap intervals are available.
#'
#' @return an edTMB object, as from \code{\link{ed}}, with the BMDs and their covariance matrix
#' @export
bmd <- function(x, bmr=0.1, def=c("excess", "additional", "relative", "sd"), linfct=NULL, newdata=NULL,
                scale=c("marginal", "conditional"), nodes=NULL){
  def <- match.arg(def)
  scale <- match.arg(scale)
  mod <- modelindex(x$model)
  m <- drcmodel(mod)
  if (is.na(m$ed) || m$ed %in% "hormesis") stop("benchmark doses are not available for the ", x$model, " model.")
  if (def %in% c("excess", "additional") && !(x$family %in% c("binomial", "betabinomial", "beta"))){
    stop("the ", def, " risk is for responses that are proportions (binomial, betabinomial or beta family); use def=\"relative\" or \"sd\".")
  }
  terms <- x$reterms
  if (def == "sd"){
    if (!(x$family %in% c("gaussian", "lognormal"))) stop("def=\"sd\" is only available for the gaussian and lognormal families.")
    if (length(terms) > 0) stop("def=\"sd\" is not available for models with random effects: the standard deviation of the response depends on the random effects.")
  }
  if (!is.numeric(bmr) || length(bmr) == 0 || any(!is.finite(bmr)) || any(bmr <= 0)) stop("bmr must be positive.")
  if (def %in% c("excess", "additional") && any(bmr >= 1)) stop("the ", def, " risk bmr must be below 1.")

  # curves and benchmark responses, paired as in ed
  rows <- NULL
  if (!is.null(newdata)){
    if (!is.null(linfct)) stop("give either linfct or newdata, not both.")
    nl <- edlinfct(x, newdata)
    rows <- rep(seq_along(nl), length(bmr))
    linfct <- rep(nl, length(bmr))
    bmr <- rep(bmr, each=length(nl))
  }
  if (is.null(linfct)) linfct <- diag(5)
  if (!is.list(linfct)) linfct <- list(linfct)
  if (length(bmr) != length(linfct)){
    if (length(bmr) == 1) bmr <- rep(bmr, length(linfct))
    if (length(linfct) == 1) linfct <- rep(linfct, length(bmr))
    if (length(bmr) != length(linfct)) stop("bmr and linfct must have the same length or one of them must be of length 1.")
  }
  if (is.null(names(linfct))) names(linfct) <- paste("lf", seq_along(linfct), sep="")
  ncf <- length(x$start)
  for (L in linfct){
    if (!is.matrix(L) || nrow(L) != 5 || ncol(L) != ncf){
      stop("each matrix in linfct must have 5 rows (b1-b5) and one column for each of the ", ncf,
           " model coefficients: ", paste(x$coefnames, collapse=", "), ".",
           if (is.null(newdata) && ncf > 5) " Alternatively, give the covariate values of the curves in newdata.")
    }
  }

  pf <- x$sdrl$par.fixed
  nb <- sum(!x$fix)
  pv <- x$start
  pv[!x$fix] <- pf[seq_len(nb)]
  ith <- which(names(pf) == "theta")
  ils <- which(names(pf) == "log_sigma")
  theta <- unname(pf[ith])
  ls <- if (length(ils) > 0) unname(pf[ils]) else 0
  uhat <- if (length(terms) > 0) x$sdrl$par.random else numeric(0)
  nu <- length(uhat)
  nbmd <- length(linfct)
  tnp <- vapply(terms, function(tm) length(tm$pars), integer(1))
  ncoef <- sum(tnp)

  # random effects of each BMD: the weights W of the fitted random effects of the groups in newdata,
  # and the terms averaged over (marginal BMDs), with the loads of their random coefficients
  W <- rep(list(matrix(0, 5, nu)), nbmd)
  integ <- matrix(scale == "marginal", nbmd, length(terms))
  loads <- lapply(terms, function(tm) reloads(tm, matrix(as.numeric(tm$cnms == "(Intercept)"), nbmd, length(tm$cnms), byrow=TRUE)))
  if (!is.null(newdata) && length(terms) > 0){
    rl <- relevels(x, newdata, skip=TRUE)
    W <- lapply(rows, function(r) reweights(x, rl, r))
    integ <- integ & rl$newgroup[rows, , drop=FALSE]
    if (scale == "marginal"){
      loads <- lapply(seq_along(terms), function(t) reloads(terms[[t]], rl$X[[t]][rows, , drop=FALSE]))
      if (any(is.na(unlist(loads)[unlist(lapply(seq_along(terms), function(t) rep(integ[, t], length(terms[[t]]$pars))))]))){
        stop("newdata has missing values in the covariates of the random slopes.")
      }
    }
  } else if (scale == "marginal" && any(unlist(lapply(terms, function(tm) tm$cnms)) != "(Intercept)")){
    stop("the marginal BMD depends on the covariates of the random slopes: give their values in newdata, or use scale=\"conditional\".")
  }
  cond <- vapply(W, function(w) any(w != 0), logical(1))

  lnk <- linkindex(x$link)
  sfac <- if (def == "sd") sqrt(varscale(x)) else 1
  edl <- numeric(nbmd)
  grad <- matrix(0, nbmd, length(pf) + nu)
  info <- vector("list", nbmd)
  for (i in seq_len(nbmd)){
    lpv <- as.vector(linfct[[i]] %*% pv) + as.vector(W[[i]] %*% uhat)
    # loads of the random coefficients averaged over, and the quadrature nodes
    A <- matrix(0, 5, ncoef)
    coff <- 0
    for (t in seq_along(terms)){
      if (integ[i, t]) for (k in seq_len(tnp[t])) A[terms[[t]]$pars[k], coff + k] <- loads[[t]][i, k]
      coff <- coff + tnp[t]
    }
    # random effects only on b1 and b2 with the identity link enter the mean linearly: the marginal BMD
    # is the conditional one
    if (all(A[setdiff(1:5, which(lnk[1:2] == 1)), ] == 0)) A[] <- 0
    act <- which(colSums(A != 0) > 0)
    quad <- bmdnodes(ncoef, act, nodes)
    # lognormal: bmr standard deviations of log(y)
    info[[i]] <- list(def=if (def == "sd" && x$family == "lognormal") "logsd" else def, bmr=bmr[i], inc=bmdinc(x, lpv), A=A, E=quad$E, wq=quad$w, sfac=sfac)
    edl[i] <- bmdvalue(x, info[[i]], lpv, theta, ls)
    if (is.na(edl[i])){
      grad[i, ] <- NA
      next
    }
    x0 <- if (length(quad$w) > 0) edl[i] else if (m$ed %in% "polynomial") edl[i] else 0
    eobj <- MakeADFun(list(model="drcBMD", mod=mod, lnk=as.integer(lnk), alpha=alphavalue(x$alpha), smooth=smoothwidth(x$smooth),
                           def=match(info[[i]]$def, c("excess", "additional", "relative", "sd", "logsd")), bmr=bmr[i], inc=as.integer(info[[i]]$inc),
                           sfac=sfac, x0=x0, tnp=as.integer(tnp), A=A, E=quad$E, wq=quad$w),
                      list(b=as.vector(linfct[[i]] %*% pv), u=as.vector(W[[i]] %*% uhat), theta=theta, log_sigma=ls),
                      ADreport=TRUE, DLL="drcTMB_TMBExports", silent=TRUE)
    edl[i] <- eobj$fn()
    egrad <- eobj$gr()
    g <- numeric(length(pf) + nu)
    g[seq_len(nb)] <- (egrad[1:5] %*% linfct[[i]])[!x$fix]
    g[ith] <- egrad[10 + seq_along(ith)]
    if (length(ils) > 0) g[ils] <- egrad[11 + length(theta)]
    g[length(pf) + seq_len(nu)] <- egrad[6:10] %*% W[[i]]
    grad[i, ] <- g
  }

  cvc <- bmdcov(x, any(cond))
  edvar <- grad %*% cvc %*% t(grad)
  edn <- paste(bmr, names(linfct), sep=":")
  if (anyNA(edl)) warning("the curves do not reach the benchmark response for ", paste(edn[is.na(edl)], collapse=", "), "; the BMDs are NA.", call.=FALSE)
  names(edl) <- edn
  colnames(edvar) <- rownames(edvar) <- edn
  integrated <- vapply(info, function(b) length(b$wq) > 0, logical(1))
  res <- list(estimates=edl, variance=edvar, fit=x, linfct=linfct, q=setNames(bmr, edn), type="bmd", random=FALSE, residual=FALSE,
              uweights=W, conditional=setNames(cond, edn), newvar=setNames(rep(FALSE, nbmd), edn),
              integrated=setNames(integrated, edn), bmd=info, def=def, scale=scale)
  class(res) <- "edTMB"
  res
}

# direction of the curve with linear predictors lpv: 1 increasing, -1 decreasing in the dose
bmdinc <- function(fit, lpv){
  m <- drcmodel(fit$model)
  b <- linkpars(linkindex(fit$link), lpv)
  s <- sign(b[2] - b[1]) * (if (m$b3 %in% "steepness" && b[3] < 0) -1 else 1)
  if (s == 0) 1 else s
}

# response at the BMD, as bmdtarget in the TMB template
bmdtarget <- function(def, bmr, inc, mu0, sigma){
  switch(def,
         excess=if (inc > 0) mu0 + bmr*(1 - mu0) else mu0*(1 - bmr),
         additional=mu0 + inc*bmr,
         relative=mu0*(1 + inc*bmr),
         sd=mu0 + inc*bmr*sigma,
         logsd=mu0*exp(inc*bmr*sigma))
}

# level of F at the BMD of the curve with linear predictors lpv (without averaging over random
# effects), with log residual standard deviation ls
bmdlevel <- function(fit, info, lpv, ls){
  b <- linkpars(linkindex(fit$link), lpv)
  (bmdtarget(info$def, info$bmr, info$inc, bmdmu0(fit, b), info$sfac*exp(ls)) - b[1]) / (b[2] - b[1])
}

# response at dose zero of the curve with parameters b (b1-b5): the asymptote at dose zero, which does
# not depend on b4, except for the models with a location b4 and the smoothed threshold of nec
bmdmu0 <- function(fit, b){
  m <- drcmodel(fit$model)
  if (m$b4 %in% "location" || smoothwidth(fit$smooth) > 0){
    return(drcmeanR(modelindex(fit$model), 0, b[1], b[2], b[3], b[4], b[5], fit$alpha, smoothwidth(fit$smooth)))
  }
  if (m$logdose && b[3] < 0) b[2] else b[1]
}

# BMD of the curve with linear predictors lpv, the variance parameters theta of the random effects and
# the log residual standard deviation ls; marginal BMDs by the root of the marginal mean; NA if the
# curve does not reach the benchmark response
bmdvalue <- function(fit, info, lpv, theta, ls){
  mod <- modelindex(fit$model)
  b <- linkpars(linkindex(fit$link), lpv)
  q <- bmdlevel(fit, info, lpv, ls)
  xc <- if (is.finite(q) && q > 0 && q < 1) edfun(mod, q, b[3], b[4], b[5], smooth=smoothwidth(fit$smooth)) else NA
  if (length(info$wq) == 0) return(if (isTRUE(xc > 0)) xc else NA)

  # marginal mean at dose d, over the linear predictors lpv + Z of the quadrature nodes
  lnk <- linkindex(fit$link)
  Z <- info$A %*% recholesky(fit, theta) %*% info$E
  B <- linkpars(lnk, t(lpv + Z))
  sm <- smoothwidth(fit$smooth)
  mm <- function(d) sum(info$wq * drcmeanR(mod, d, B[, 1], B[, 2], B[, 3], B[, 4], B[, 5], fit$alpha, sm))
  target <- bmdtarget(info$def, info$bmr, info$inc, mm(0), info$sfac*exp(ls))
  g <- function(d) info$inc * (mm(d) - target)
  # bracket the root, starting from the BMD of the curve lpv, or from the mean positive dose
  x0 <- if (isTRUE(xc > 0)) xc else mean(fit$obj$env$data$x[fit$obj$env$data$x > 0])
  lo <- hi <- x0
  for (k in 1:60){
    if (isTRUE(g(lo) < 0)) break
    lo <- lo/2
  }
  for (k in 1:60){
    if (isTRUE(g(hi) > 0)) break
    hi <- hi*2
  }
  if (!isTRUE(g(lo) < 0) || !isTRUE(g(hi) > 0)) return(NA)
  uniroot(g, c(lo, hi), tol=1e-12*x0)$root
}

# Cholesky factor of the block-diagonal covariance matrix of the random coefficients of all terms, from
# their parameters theta, as in the drcBMD template
recholesky <- function(fit, theta){
  tnp <- vapply(fit$reterms, function(tm) length(tm$pars), integer(1))
  C <- matrix(0, sum(tnp), sum(tnp))
  toff <- 0
  coff <- 0
  for (np in tnp){
    L <- diag(np)
    L[upper.tri(L)] <- theta[toff + np + seq_len(np*(np - 1)/2)]
    L <- t(L)
    L <- L / sqrt(rowSums(L^2))
    C[coff + seq_len(np), coff + seq_len(np)] <- exp(theta[toff + seq_len(np)]) * L
    toff <- toff + np + np*(np - 1)/2
    coff <- coff + np
  }
  C
}

# Gauss-Hermite nodes and weights for the standard normal distribution: a product grid over the random
# coefficients act of ncoef, with nodes per coefficient (default by their number); nodes with
# negligible weights are dropped; no nodes without random coefficients to average over
bmdnodes <- function(ncoef, act, nodes=NULL){
  d <- length(act)
  if (d == 0) return(list(E=matrix(0, ncoef, 0), w=numeric(0)))
  k <- if (is.null(nodes)) c(30, 15, 9, 7, 5)[min(d, 5)] else nodes
  gh <- ghnormal(k)
  grid <- as.matrix(expand.grid(rep(list(gh$x), d)))
  w <- apply(as.matrix(expand.grid(rep(list(gh$w), d))), 1, prod)
  keep <- w > 1e-14*max(w)
  E <- matrix(0, ncoef, sum(keep))
  E[act, ] <- t(grid[keep, , drop=FALSE])
  list(E=E, w=w[keep]/sum(w[keep]))
}

# Gauss-Hermite nodes and weights for the standard normal distribution, by the Golub-Welsch algorithm
ghnormal <- function(k){
  if (k == 1) return(list(x=0, w=1))
  J <- matrix(0, k, k)
  off <- sqrt(seq_len(k - 1))
  J[cbind(1:(k - 1), 2:k)] <- off
  J[cbind(2:k, 1:(k - 1))] <- off
  e <- eigen(J, symmetric=TRUE)
  list(x=e$values, w=e$vectors[1, ]^2)
}

# covariance of all fixed parameters (coefficients, variance parameters of the random effects, residual
# parameter) and the random effects: jointly for BMDs of fitted groups, otherwise of the fixed parameters
bmdcov <- function(object, cond){
  nu <- length(object$sdrl$par.random)
  if (cond) return(jointcov(object, all=TRUE))
  as.matrix(bdiag(fixedcov(object), matrix(0, nu, nu)))
}

# theta and log residual parameter in a parameter vector par of the fit: the fixed parameters, or all
# parameters with the random effects
bmdpars <- function(fit, par){
  env <- fit$obj$env
  nm <- if (length(par) == length(env$par)) names(env$par) else names(fit$opt$par)
  ls <- par[nm == "log_sigma"]
  list(theta=unname(par[nm == "theta"]), ls=if (length(ls) > 0) unname(ls[1]) else 0)
}

# BMDs parm of an edTMB object with profile likelihood intervals: not averaged over random effects, and
# the response at dose zero does not depend on b4
bmdprofile <- function(object, parm){
  fit <- object$fit
  m <- drcmodel(fit$model)
  !any(object$integrated[parm]) && !(m$b4 %in% "location") && smoothwidth(fit$smooth) == 0
}
