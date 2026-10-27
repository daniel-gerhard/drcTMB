#' Effective dose estimation
#'
#' @param x a drmTMB object
#' @param respLev effective dose level in percent, default is 50. The effective dose ED p is the dose at which the response has changed by p percent of the difference between the asymptote at low doses (the control) and the asymptote at high doses.
#' @param linfct a matrix of linear functions to identify the parameters of the dose-response curve, with 5 rows (b1-b5) and one column for each model coefficient, or a list of such matrices; default is the identity matrix, for models without covariates in fform. Each matrix is paired with an entry of respLev.
#' @param newdata a data frame with the values of the covariates in fform, one row for each dose-response curve, as an alternative to linfct; the dose is not needed. The effective doses are computed at each level of respLev for each row, named by the covariate values (or by the row names of newdata, if set).
#' @param random logical, if TRUE, the variance-covariance matrix of the random effects is included in the standard error estimation
#' @param residual logical, if TRUE, the residual variance is included in the standard error estimation (gaussian family only)
#'
#' @return a edTMB object
#' @export
#'
ed <- function(x, respLev=50, linfct=NULL, newdata=NULL, random=TRUE, residual=FALSE){
  if (!is.null(newdata)){
    if (!is.null(linfct)) stop("give either linfct or newdata, not both.")
    # every level of respLev for every curve
    nl <- edlinfct(x, newdata)
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

  est <- x$sdrl$par.fixed
  covest <- x$sdrl$cov.fixed
  rn <- names(est)
  id <- rn %in% paste("b", 1:5, sep="")
  est <- est[id,drop=FALSE]
  covest <- covest[id,id,drop=FALSE]
  
  pv <- numeric(length=length(x$start))
  pv[!x$fix] <- est
  pv[x$fix] <- x$start[x$fix]
  pvar <- matrix(0, nrow=length(x$start), ncol=length(x$start))
  pvar[!x$fix, !x$fix] <- covest

  ests <- x$ssdrl

  # random effects: sum of the covariance matrices of all terms
  rvc <- matrix(0, nrow=5, ncol=5)
  if (!is.null(x$rform) & random == TRUE){
    for (vc in recov(x)){
      rvc[vc$pars, vc$pars] <- rvc[vc$pars, vc$pars] + vc$cov
    }
  }

  if (residual == TRUE && x$family != "gaussian"){
    stop("residual=TRUE is only available for the gaussian family.")
  }
  if (residual == TRUE){
      rvar <- exp(x$sdrl$par.fixed["log_sigma"])^2
  } else {
      rvar <- 0
  }

  mod <- which(x$model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
  lnk <- which(x$link == c("identity", "log", "logit"))

  # gradient of each ED w.r.t. (parameters, u1-u5, r)
  edl <- numeric(length=length(linfct))
  ql <- numeric(length=length(linfct))
  grad <- matrix(0, nrow=length(linfct), ncol=length(pv) + 6)

  for (i in seq_along(linfct)){
    lpv <- as.vector(linfct[[i]] %*% pv)

    # level of the increasing function F of the model at the ED: the control is b1 for an
    # increasing curve (b3 > 0) and for the models without b3, and b2 for a decreasing curve
    q <- respLev[i]/100
    if (!(mod %in% c(6, 7, 8)) && lpv[3] < 0) q <- 1 - q
    ql[i] <- q

    # ED and its gradient w.r.t. (b1-b5, u1-u5, r) from the drcED TMB template
    eobj <- MakeADFun(list(model="drcED", mod=mod, lnk=lnk, lnk4=b4linkindex(x), q=q),
                      list(b=lpv, u=rep(0, 5), r=0),
                      ADreport=TRUE, DLL="drcTMB_TMBExports", silent=TRUE)
    edl[i] <- eobj$fn()
    egrad <- eobj$gr()
    grad[i, ] <- c(egrad[1:5] %*% linfct[[i]], egrad[6:11])
  }

  cvc <- as.matrix(bdiag(pvar, rvc, rvar))
  edvar <- grad %*% cvc %*% t(grad)
  names(edl) <- paste(respLev, names(linfct), sep=":")
  colnames(edvar) <- rownames(edvar) <- paste(respLev, names(linfct), sep=":")
  # the fit, linear functions and levels of F are kept for profile
  names(ql) <- names(edl)
  res <- list(estimates=edl, variance=edvar, fit=x, linfct=linfct, q=ql, random=random, residual=residual)
  class(res) <- "edTMB"
  return(res)
}


# linear functions of the curves in the rows of newdata: the rows of the design matrices of b1-b5,
# named by the values of the covariates in fform, or by the row names of newdata if set
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
  if (.row_names_info(newdata) > 0){
    names(L) <- rownames(newdata)
  } else if (length(vars) > 0){
    names(L) <- do.call(paste, c(lapply(vars, function(v) paste(v, newdata[[v]], sep="=")), sep=","))
  } else {
    names(L) <- seq_len(nrow(newdata))
  }
  L
}
