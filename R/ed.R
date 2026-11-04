#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level in percent, default is 50. The effective dose ED p is the dose at which the response has changed by p percent of the difference between the asymptote at low doses (the control) and the asymptote at high doses; for the hormesis models, the dose above the peak of the stimulation, see Details. For \code{type="absolute"}, values of the response.
#' @param linfct a matrix of linear functions to identify the parameters of the dose-response curve, with 5 rows (b1-b5) and one column for each model coefficient, or a list of such matrices; default is the identity matrix, for models without covariates in fform. Each matrix is paired with an entry of respLev.
#' @param newdata a data frame with the values of the covariates in fform, one row for each dose-response curve, as an alternative to linfct; the dose is not needed. With grouping variables of the random effects, the curves of fitted groups, see Details. The effective doses are computed at each level of respLev for each row, named by the covariate and group values (or by the row names of newdata, if set).
#' @param random logical, if TRUE, the effective doses of rows of newdata with fitted groups are conditional on their random effects, and for the other random effect terms, the variance-covariance matrix of the random effects is included in the standard error estimation. If FALSE, the effective doses are for the population curves.
#' @param residual logical, if TRUE, the residual variance is included in the standard error estimation (gaussian family only)
#' @param type "relative" (default) for effective doses relative to the asymptotes, or "absolute" for the doses at which the curves reach the responses in respLev, on the scale of the response (proportions for binomial, beta-binomial and beta). Absolute effective doses depend on the asymptotes, and on their random effects; responses outside the range of a curve give NA, with a warning.
#'
#' @details Without random effects, or with \code{random=FALSE}, the effective doses are those of the
#' population curves, with all random effects at zero, and their standard errors include the
#' uncertainty of the fixed effects only. With \code{random=TRUE}, the standard errors also include the
#' variance of the random effects, the variation of the effective dose between groups.
#'
#' If a row of newdata gives a level of a grouping factor that is in the fitted data, the effective
#' dose of that row is for the curve of the group, with its predicted random effects (the conditional
#' modes), as in \code{\link{predict.drmTMB}}. These are shrunken towards the population curve, more so
#' for groups with few observations. Their standard errors include the uncertainty of the fixed
#' effects and the prediction error of the random effects, from the joint precision matrix of
#' \code{\link[TMB]{sdreport}}. Random effect terms without a grouping variable in newdata, or with a
#' level that is not in the fitted data, are treated as above. Profile intervals are not available
#' for the effective doses of fitted groups, see \code{\link{confint.edTMB}}.
#'
#' For the hormesis models braincousens and cedergreen, the effective dose ED p is the dose above the
#' peak of the stimulation at which the response has decreased from the response at dose zero (b2) by
#' p percent of the difference to the lower asymptote b1, as in drc; with \code{type="absolute"}, the
#' dose above the peak at which the curve reaches the response. It has no closed form: it is found
#' numerically, and its derivatives, for the standard errors, follow from the implicit function
#' theorem. The models need a decreasing curve, with b3 < 0. With \code{respLev=0}, the effective
#' dose is the dose above the peak at which the response returns to the response at dose zero, the
#' limit of the stimulation (the zero equivalent point); see \code{\link{peak}} for the peak itself.
#'
#' @return a edTMB object
#' @export
#'
ed <- function(x, respLev=50, linfct=NULL, newdata=NULL, random=TRUE, residual=FALSE, type=c("relative", "absolute")){
  type <- match.arg(type)
  if (is.na(drcmodel(x$model)$ed)) stop("the ", x$model, " model has no effective doses.")
  edcore(x, respLev, linfct, newdata, random, residual, type)
}

#' Peak of hormesis and substrate inhibition curves
#'
#' The dose at which a curve of the hormesis models braincousens and cedergreen, or of the substrate
#' inhibition model haldane, reaches its maximum, and the maximal response, as \code{MAX} in drc. For the
#' hormesis models the peak has no closed form: it is found numerically, and its derivatives, for the
#' standard errors, follow from the implicit function theorem, as for the effective doses of
#' \code{\link{ed}}. Curves without a stimulation (b5 <= 0, or a curve that does not decrease at high
#' doses) have no peak, and give NA with a warning. For haldane, the peak is at the dose sqrt(b4/b5), the
#' optimal substrate concentration sqrt(Km Ki); without inhibition (b5 = 0) there is no peak.
#'
#' The result is an edTMB object, as from \code{\link{ed}}, with the dose and the response at the
#' peak of each curve, named "dose" and "response", and their joint covariance matrix, e.g. for
#' \code{glht} of multcomp. \code{\link{confint.edTMB}} gives Wald and bootstrap intervals for both,
#' and profile likelihood intervals for the dose; for the response, profile intervals are not
#' available and Wald intervals are returned, with a message.
#'
#' @param x a drmTMB object of the braincousens, cedergreen or haldane model
#' @param linfct,newdata,random the curves, as in \code{\link{ed}}
#'
#' @return an edTMB object
#' @export
peak <- function(x, linfct=NULL, newdata=NULL, random=TRUE){
  m <- drcmodel(x$model)
  if (is.null(m$hormesis) && is.null(m$peak)) stop("peak is only available for the hormesis and substrate inhibition models ", modelswith(function(m) !is.null(m$hormesis) || !is.null(m$peak)), ".")
  # the dose (1) and the response (2) at the peak of every curve
  if (is.null(newdata)){
    if (is.null(linfct)) linfct <- diag(5)
    if (!is.list(linfct)) linfct <- list(linfct)
    if (is.null(names(linfct))) names(linfct) <- paste("lf", seq_along(linfct), sep="")
    edcore(x, rep(1:2, each=length(linfct)), rep(linfct, 2), NULL, random, FALSE, "peak")
  } else {
    edcore(x, 1:2, linfct, newdata, random, FALSE, "peak")
  }
}

# effective doses (type "relative" or "absolute", respLev the levels) or peaks of hormesis curves
# (type "peak", respLev 1 for the dose and 2 for the response at the peak) of the curves in linfct
# or newdata, with their covariance by the delta method
edcore <- function(x, respLev, linfct, newdata, random, residual, type){
  absolute <- type == "absolute"
  peak <- type == "peak"
  if (!is.null(newdata)){
    if (!is.null(linfct)) stop("give either linfct or newdata, not both.")
    # every level of respLev for every curve
    nl <- edlinfct(x, newdata)
    rows <- rep(seq_along(nl), length(respLev))
    linfct <- rep(nl, length(respLev))
    respLev <- rep(respLev, each=length(nl))
  }
  if (is.null(linfct)) linfct <- diag(5)
  if (!is.list(linfct)) linfct <- list(linfct)
  if (length(respLev) != length(linfct)){
    if (length(respLev) == 1){
      respLev <- rep(respLev, length(linfct))
    }
    if (length(linfct) == 1){
      linfct <- rep(linfct, length(respLev))
    } 
    if (length(respLev) != length(linfct)){
      stop("respLev and linfct must have the same length or one of them must be of length 1.")
    }
  }
  if (is.null(names(linfct))) names(linfct) <- paste("lf", 1:length(linfct), sep="")
  ncf <- length(x$start)
  for (L in linfct){
    if (!is.matrix(L) || nrow(L) != 5 || ncol(L) != ncf){
      stop("each matrix in linfct must have 5 rows (b1-b5) and one column for each of the ", ncf,
           " model coefficients: ", paste(x$coefnames, collapse=", "), ".",
           if (is.null(newdata) && ncf > 5) " Alternatively, give the covariate values of the curves in newdata.")
    }
  }

  if (residual == TRUE && x$family != "gaussian"){
    stop("residual=TRUE is only available for the gaussian family.")
  }
  rvar <- if (residual) residpar(x)^2 else 0

  pv <- x$start
  pv[!x$fix] <- x$sdrl$par.fixed[seq_len(sum(!x$fix))]
  nb <- sum(!x$fix)
  terms <- x$reterms
  uhat <- if (length(terms) > 0) x$sdrl$par.random else numeric(0)
  nu <- length(uhat)

  # random effects of each effective dose: the weights W of the fitted random effects on b1-b5, for
  # the rows of newdata with a fitted group (random=TRUE), and the terms whose variance is added
  # for a new group (random=TRUE, no fitted level)
  ned <- length(linfct)
  W <- rep(list(matrix(0, 5, nu)), ned)
  newterm <- matrix(random, ned, length(terms))
  # columns of each term for a new group (intercept, covariates of slopes), for each ED
  loads <- lapply(terms, function(tm) reloads(tm, matrix(as.numeric(tm$cnms == "(Intercept)"), ned, length(tm$cnms), byrow=TRUE)))
  if (!is.null(newdata) && random && length(terms) > 0){
    rl <- relevels(x, newdata, skip=TRUE)
    W <- lapply(rows, function(r) reweights(x, rl, r))
    newterm <- rl$newgroup[rows, , drop=FALSE]
    loads <- lapply(seq_along(terms), function(t) reloads(terms[[t]], rl$X[[t]][rows, , drop=FALSE]))
    if (any(is.na(unlist(loads)))) stop("newdata has missing values in the covariates of the random slopes.")
  } else if (random && any(unlist(lapply(terms, function(tm) tm$cnms)) != "(Intercept)")){
    stop("the variance of the effective doses depends on the covariates of the random slopes: give their values in newdata, or use random=FALSE.")
  }
  cond <- vapply(W, function(w) any(w != 0), logical(1))

  mod <- modelindex(x$model)
  m <- drcmodel(mod)
  lnk <- linkindex(x$link)

  # gradient of each ED w.r.t. the free coefficients, the fitted random effects, the random effects
  # of a new group for each term, and the residual
  edl <- numeric(ned)
  ql <- numeric(ned)
  grad <- matrix(0, ned, nb + nu + sum(vapply(terms, function(tm) length(tm$pars), integer(1))) + 1)

  for (i in seq_len(ned)){
    lpv <- as.vector(linfct[[i]] %*% pv)
    u <- as.vector(W[[i]] %*% uhat)
    # parameters b1-b5 of the curve
    b <- linkpars(lnk, lpv + u)

    # level of the increasing function F of the model at the ED: the control is b1 for an
    # increasing curve (b3 > 0) and for the models without b3, and b2 for a decreasing curve;
    # for absolute effective doses, the response, and the level of F follows from the asymptotes;
    # for peaks, the target
    if (peak){
      q <- respLev[i]
    } else if (absolute){
      q <- respLev[i]
      qf <- edq(x, q, lpv + u, TRUE)
      if (!is.finite(qf) || qf <= 0 || qf >= 1){
        edl[i] <- NA
        grad[i, ] <- NA
        ql[i] <- q
        next
      }
    } else {
      q <- respLev[i]/100
      if (m$b3 %in% "steepness" && b[3] < 0) q <- 1 - q
    }
    ql[i] <- q

    # hormesis models: log dose of the effective dose above the peak, or of the peak, found numerically
    s0 <- 0
    if (m$ed %in% "hormesis"){
      if (b[3] >= 0) stop(if (peak) "the peaks" else "the effective doses", " of the ", x$model, " model need a decreasing curve, with b3 < 0.")
      s0 <- if (peak) hormpeak(mod, b[3], b[4], hormc(x, lpv + u), x$alpha)
            else hormroot(mod, edq(x, q, lpv + u, absolute), b[3], b[4], hormc(x, lpv + u), x$alpha)
      if (is.na(s0)){
        edl[i] <- NA
        grad[i, ] <- NA
        next
      }
    }
    # multistage model: log of the smallest positive root of the polynomial
    # closed-form peaks: NA without a peak
    if (peak && !is.null(m$peak) && is.na(m$peak$dose(b))){
      edl[i] <- NA
      grad[i, ] <- NA
      next
    }
    if (m$ed %in% "polynomial"){
      s0 <- log(msroot(edq(x, q, lpv + u, absolute), b[3:5]))
      if (is.na(s0)){
        edl[i] <- NA
        grad[i, ] <- NA
        next
      }
    }

    # ED and its gradient w.r.t. (b1-b5, u1-u5, r) from the drcED TMB template
    eobj <- MakeADFun(list(model="drcED", mod=mod, lnk=as.integer(lnk), absolute=as.integer(absolute), q=q,
                           alpha=alphavalue(x$alpha), smooth=smoothwidth(x$smooth), s0=s0, target=as.integer(if (peak) q else 0)),
                      list(b=lpv, u=u, r=0),
                      ADreport=TRUE, DLL="drcTMB_TMBExports", silent=TRUE)
    edl[i] <- eobj$fn()
    egrad <- eobj$gr()
    gu <- egrad[6:10]
    gnew <- unlist(lapply(seq_along(terms), function(t) gu[terms[[t]]$pars] * loads[[t]][i, ] * newterm[i, t]))
    grad[i, ] <- c((egrad[1:5] %*% linfct[[i]])[!x$fix], gu %*% W[[i]], gnew, egrad[11])
  }

  # covariance of the free coefficients and the fitted random effects: jointly for effective doses of
  # fitted groups, otherwise of the coefficients only
  cbu <- if (any(cond)) jointcov(x) else as.matrix(bdiag(fixedcov(x)[seq_len(nb), seq_len(nb), drop=FALSE], matrix(0, nu, nu)))
  cvc <- as.matrix(bdiag(c(list(cbu), lapply(recov(x), function(vc) vc$cov), list(matrix(rvar)))))
  edvar <- grad %*% cvc %*% t(grad)
  edn <- paste(if (absolute) paste("y=", respLev, sep="") else if (peak) c("dose", "response")[respLev] else respLev, names(linfct), sep=":")
  if (anyNA(edl)){
    if (peak) warning("the curves ", paste(unique(names(linfct)[is.na(edl)]), collapse=", "), " have no peak; the estimates are NA.", call.=FALSE)
    else warning("the curves do not reach the response for ", paste(edn[is.na(edl)], collapse=", "), "; the effective doses are NA.", call.=FALSE)
  }
  names(edl) <- edn
  colnames(edvar) <- rownames(edvar) <- edn
  # the fit, linear functions, random effect weights and levels of F are kept for profile and confint
  names(ql) <- names(edl)
  res <- list(estimates=edl, variance=edvar, fit=x, linfct=linfct, q=ql, type=type, random=random, residual=residual,
              uweights=W, conditional=setNames(cond, names(edl)),
              newvar=setNames(apply(newterm, 1, any) & length(terms) > 0, names(edl)))
  class(res) <- "edTMB"
  return(res)
}

# level of F at an effective dose with level q (absolute=FALSE) or response q (absolute=TRUE), for the
# parameters lpv (b1-b5 on the scale of the linear predictors) of the curve
edq <- function(fit, q, lpv, absolute){
  if (!absolute) return(q)
  b <- linkpars(linkindex(fit$link), lpv)
  (q - b[1]) / (b[2] - b[1])
}

# hormesis models: the size of the hormesis effect b5 relative to the difference of the asymptotes,
# for the parameters lpv (b1-b5 on the scale of the linear predictors) of the curve
hormc <- function(fit, lpv){
  b <- linkpars(linkindex(fit$link), lpv)
  b[5] / (b[2] - b[1])
}

# log dose of the effective dose of a hormesis model with steepness b3 < 0, b4,
# relative hormesis effect c5 (see hormc) and alpha: the dose above the peak at which
# F(x) (1 + c5 g(x)) equals q, with the loglogistic F. The last sign change on a grid of log doses,
# from where F is close to 1 up to a dose at which the curve is below q, is refined by uniroot. NA
# if the curve does not reach q.
hormroot <- function(mod, q, b3, b4, c5, alpha){
  if (!is.finite(b3) || b3 >= 0 || !is.finite(q) || !is.finite(c5) || !is.finite(b4) || b4 <= 0) return(NA)
  G <- function(s) plogis(b3*(s - log(b4))) * (1 + c5*hormg(mod, exp(s), alpha)) - q
  s <- hormgrid(G, b3, b4)
  if (is.null(s)) return(NA)
  gs <- G(s)
  i <- which(gs[-length(gs)] >= 0 & gs[-1] < 0)
  if (length(i) == 0) return(NA)
  i <- max(i)
  uniroot(G, s[i + 0:1], tol=1e-12)$root
}

# effective dose of the multistage model at level q of F, for the coefficients beta of the polynomial:
# the smallest positive root of beta1 x + beta2 x^2 + beta3 x^3 = -log(1 - q); NA without one
msroot <- function(q, beta){
  if (!is.finite(q) || q <= 0 || q >= 1 || any(!is.finite(beta))) return(NA)
  r <- polyroot(c(log(1 - q), beta))
  r <- Re(r)[abs(Im(r)) < 1e-6*pmax(1, abs(Re(r))) & Re(r) > 0]
  if (length(r) == 0) NA else min(r)
}

# grid of 2001 log doses for a hormesis curve with b3 < 0 and b4, from where F is close to 1 to a dose
# at which G is negative, doubling the distance from b4 up to at most exp(700); NULL if G stays
# non-negative
hormgrid <- function(G, b3, b4){
  s4 <- log(b4)
  w <- 1/abs(b3)
  d <- 5*w
  while (isTRUE(G(s4 + d) >= 0) && s4 + d < 700) d <- 2*d
  su <- min(s4 + d, 700)
  if (!isTRUE(G(su) < 0)) return(NULL)
  seq(s4 - 30*w, su, length.out=2001)
}

# log dose of the peak of a hormesis curve with b3 < 0, b4, relative hormesis effect
# c5 (see hormc) and alpha: the maximum of F(x) (1 + c5 g(x)) on a grid of log doses, up to where the
# curve is back below the response at dose zero, refined by optimise. NA without a stimulation
hormpeak <- function(mod, b3, b4, c5, alpha){
  if (!is.finite(b3) || b3 >= 0 || !is.finite(c5) || c5 <= 0 || !is.finite(b4) || b4 <= 0) return(NA)
  H <- function(s) plogis(b3*(s - log(b4))) * (1 + c5*hormg(mod, exp(s), alpha))
  s <- hormgrid(function(s) H(s) - 1, b3, b4)
  if (is.null(s)) return(NA)
  i <- which.max(H(s))
  if (length(i) == 0 || i == 1 || i == length(s)) return(NA)
  optimise(H, s[c(i - 1, i + 1)], maximum=TRUE, tol=1e-12)$maximum
}

# weights of the fitted random effects u on b1-b5 for row r of the random effect design rl of newdata
reweights <- function(object, rl, r){
  w <- matrix(0, 5, length(object$sdrl$par.random))
  uoff <- 0
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    for (k in seq_along(tm$pars)) w[tm$pars[k], uoff + (k - 1)*nl + seq_len(nl)] <- as.vector(rl$Zl[[t]][[tm$zcol[k]]][r, ])
    uoff <- uoff + nl*length(tm$pars)
  }
  w
}

# joint covariance of the free coefficients b1-b5 and the prediction errors of the random effects,
# from the joint precision matrix of sdreport; the likelihood is unchanged, so the Hessian of the fit
# is reused
jointcov <- function(object){
  env <- object$obj$env
  H <- if (isTRUE(object$sdrl$pdHess)) solve(object$sdrl$cov.fixed)
  sdr <- sdreport(object$obj, par.fixed=object$opt$par, hessian.fixed=H, getJointPrecision=TRUE)
  id <- c(seq_len(sum(!object$fix)), env$random)
  C <- tryCatch(solve(sdr$jointPrecision), error=function(e) NULL)
  if (is.null(C)) stop("the joint precision matrix of the coefficients and random effects is singular; check the convergence of the fit.")
  as.matrix(C)[id, id, drop=FALSE]
}

# linear functions of the curves in the rows of newdata: the rows of the design matrices of b1-b5,
# named by the values of the covariates in fform and the grouping variables in newdata, or by the
# row names of newdata if set
edlinfct <- function(x, newdata){
  if (!is.data.frame(newdata) || nrow(newdata) == 0) stop("newdata must be a data frame with at least one row.")
  vars <- unique(unlist(lapply(x$fterms, function(ft) all.vars(ft$terms))))
  miss <- setdiff(vars, names(newdata))
  if (length(miss) > 0) stop("newdata needs the covariates of fform: ", paste(miss, collapse=", "), ".")
  Xs <- lapply(x$fterms, function(ft) designmatrix(ft, newdata))
  if (any(sapply(Xs, anyNA))) stop("newdata has missing values in the covariates of fform.")
  npi <- rep(1:5, sapply(Xs, ncol))
  L <- lapply(seq_len(nrow(newdata)), function(r){
    L <- matrix(0, 5, length(npi))
    for (i in 1:5) L[i, npi == i] <- Xs[[i]][r, ]
    L
  })
  # grouping variables of the random effects in newdata
  gvars <- unique(unlist(lapply(x$reterms, function(tm) all.vars(tm$formula))))
  vars <- c(vars, intersect(setdiff(gvars, vars), names(newdata)))
  if (.row_names_info(newdata) > 0){
    names(L) <- rownames(newdata)
  } else if (length(vars) > 0){
    names(L) <- do.call(paste, c(lapply(vars, function(v) paste(v, newdata[[v]], sep="=")), sep=","))
  } else {
    names(L) <- seq_len(nrow(newdata))
  }
  L
}
