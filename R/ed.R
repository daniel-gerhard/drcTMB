#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level in percent, default is 50. The effective dose ED p is the dose at which the response has changed by p percent of the difference between the asymptote at low doses (the control) and the asymptote at high doses.
#' @param linfct a matrix of linear functions to identify the parameters of the dose-response curve, with 5 rows (b1-b5) and one column for each model coefficient, or a list of such matrices; default is the identity matrix, for models without covariates in fform. Each matrix is paired with an entry of respLev.
#' @param newdata a data frame with the values of the covariates in fform, one row for each dose-response curve, as an alternative to linfct; the dose is not needed. With grouping variables of the random effects, the curves of fitted groups, see Details. The effective doses are computed at each level of respLev for each row, named by the covariate and group values (or by the row names of newdata, if set).
#' @param random logical, if TRUE, the effective doses of rows of newdata with fitted groups are conditional on their random effects, and for the other random effect terms, the variance-covariance matrix of the random effects is included in the standard error estimation. If FALSE, the effective doses are for the population curves.
#' @param residual logical, if TRUE, the residual variance is included in the standard error estimation (gaussian family only)
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
#' @return a edTMB object
#' @export
#'
ed <- function(x, respLev=50, linfct=NULL, newdata=NULL, random=TRUE, residual=FALSE){
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
  rvar <- if (residual) exp(x$sdrl$par.fixed[["log_sigma"]])^2 else 0

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
  if (!is.null(newdata) && random && length(terms) > 0){
    rl <- relevels(x, newdata, skip=TRUE)
    W <- lapply(rows, function(r) reweights(x, rl, r))
    newterm <- rl$newgroup[rows, , drop=FALSE]
  }
  cond <- vapply(W, function(w) any(w != 0), logical(1))

  mod <- modelindex(x$model)
  lnk <- which(x$link == c("identity", "log", "logit"))

  # gradient of each ED w.r.t. the free coefficients, the fitted random effects, the random effects
  # of a new group for each term, and the residual
  edl <- numeric(ned)
  ql <- numeric(ned)
  grad <- matrix(0, ned, nb + nu + sum(vapply(terms, function(tm) length(tm$pars), integer(1))) + 1)

  for (i in seq_len(ned)){
    lpv <- as.vector(linfct[[i]] %*% pv)
    u <- as.vector(W[[i]] %*% uhat)

    # level of the increasing function F of the model at the ED: the control is b1 for an
    # increasing curve (b3 > 0) and for the models without b3, and b2 for a decreasing curve
    q <- respLev[i]/100
    if (!(mod %in% c(6, 7, 8)) && lpv[3] + u[3] < 0) q <- 1 - q
    ql[i] <- q

    # ED and its gradient w.r.t. (b1-b5, u1-u5, r) from the drcED TMB template
    eobj <- MakeADFun(list(model="drcED", mod=mod, lnk=lnk, lnk4=b4linkindex(x), q=q),
                      list(b=lpv, u=u, r=0),
                      ADreport=TRUE, DLL="drcTMB_TMBExports", silent=TRUE)
    edl[i] <- eobj$fn()
    egrad <- eobj$gr()
    gu <- egrad[6:10]
    gnew <- unlist(lapply(seq_along(terms), function(t) gu[terms[[t]]$pars] * newterm[i, t]))
    grad[i, ] <- c((egrad[1:5] %*% linfct[[i]])[!x$fix], gu %*% W[[i]], gnew, egrad[11])
  }

  # covariance of the free coefficients and the fitted random effects: jointly for effective doses of
  # fitted groups, otherwise of the coefficients only
  cbu <- if (any(cond)) jointcov(x) else as.matrix(bdiag(x$sdrl$cov.fixed[seq_len(nb), seq_len(nb), drop=FALSE], matrix(0, nu, nu)))
  cvc <- as.matrix(bdiag(c(list(cbu), lapply(recov(x), function(vc) vc$cov), list(matrix(rvar)))))
  edvar <- grad %*% cvc %*% t(grad)
  names(edl) <- paste(respLev, names(linfct), sep=":")
  colnames(edvar) <- rownames(edvar) <- paste(respLev, names(linfct), sep=":")
  # the fit, linear functions, random effect weights and levels of F are kept for profile and confint
  names(ql) <- names(edl)
  res <- list(estimates=edl, variance=edvar, fit=x, linfct=linfct, q=ql, random=random, residual=residual,
              uweights=W, conditional=setNames(cond, names(edl)),
              newvar=setNames(apply(newterm, 1, any) & length(terms) > 0, names(edl)))
  class(res) <- "edTMB"
  return(res)
}

# weights of the fitted random effects u on b1-b5 for row r of the random effect design rl of newdata
reweights <- function(object, rl, r){
  w <- matrix(0, 5, length(object$sdrl$par.random))
  uoff <- 0
  for (t in seq_along(object$reterms)){
    tm <- object$reterms[[t]]
    nl <- length(tm$levels)
    z <- as.vector(rl$Zl[[t]][r, ])
    for (k in seq_along(tm$pars)) w[tm$pars[k], uoff + (k - 1)*nl + seq_len(nl)] <- z
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
  as.matrix(solve(sdr$jointPrecision))[id, id, drop=FALSE]
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
