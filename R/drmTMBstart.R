#' Starting values for a drmTMB model
#'
#' Computes starting values for the fixed-effect parameters of a dose-response model,
#' which can be inspected, modified, and passed to the \code{start} argument of \code{\link{drmTMB}}.
#'
#' Method "lm" takes the asymptotes from the mean response at the lowest and highest dose,
#' transforms the response to the relative scale between the asymptotes, and estimates
#' b3 and b4 by linear regression on a linearised version of the model function.
#' The asymmetry parameter b5 starts at 0. With covariates in \code{fform}, the pooled
#' starting value of each parameter is projected onto its design matrix, giving the
#' intercept or each group the pooled value and contrasts a value of zero.
#'
#' Method "grid" evaluates the negative log-likelihood (without random effects) at all
#' combinations of a grid of parameter values and returns the best grid point; parameters
#' not included in the grid start at the values of method "lm".
#' Without bounds in \code{gsl} and \code{gsu}, default ranges are used for b3 (the sign
#' of the "lm" starting value, with 1/4 to 4 times its steepness, log-spaced), b4 (the dose
#' range extended at both ends; log-spaced for models on the log-dose scale), and b5
#' (-2 to 2, logistic and loglogistic only). Default ranges are only used for parameters
#' with a single coefficient in \code{fform}; the asymptotes b1 and b2 are only included
#' in the grid with user-specified bounds.
#'
#' @param form formula with name of response and dose scale
#' @param fform fixed effects for the model parameters: a formula of the form b2 + b4 ~ Stress, or a list of such formulas, e.g. list(b2 ~ Stress, b4 ~ Stress + Time). Parameters not listed have a single intercept (~ 1). A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~ 1), is also accepted.
#' @param data a data.frame object
#' @param start a vector with starting values for each model parameter. Only NA entries are replaced.
#' @param fix a binary vector with same length as start. If TRUE, a model coefficient is fixed to the starting value. Parameters not used by the model (b5 except for logistic and loglogistic, b3 for michaelismenten, asymptoticregression and exponentialdecay) are always fixed, at 0 unless given in start.
#' @param weights a vector of weights for each observation.
#' @param family Specify the distributional assumptions. One of "gaussian", "binomial", "beta", and "poisson".
#' @param model Specify the nonlinear model. One of "logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz".
#' @param link Specify the link function for estimating asymptote parameters. One of "identity", "log", "logit".
#' @param method one of "lm" or "grid".
#' @param gsl a vector of lower bounds for the grid search, same length as start; NA (or NULL) for the default range.
#' @param gsu a vector of upper bounds for the grid search, same length as start; NA (or NULL) for the default range.
#' @param gsn number of grid points for each parameter, a single number or a vector with the same length as start (default 10). A parameter with gsn = 1 stays at its "lm" starting value.
#'
#' @return a named vector of starting values
#' @export
#'
drmTMBstart <- function(form, fform=NULL, data, start=NULL, fix=NULL, weights=NULL, family="gaussian", model="logistic", link="identity", method="lm", gsl=NULL, gsu=NULL, gsn=NULL){
  dd <- drmTMBdata(form, fform=fform, data=data, weights=weights, family=family, model=model, link=link)
  if (is.null(start)) start <- rep(NA, sum(dd$npars))
  if (is.null(fix)) fix <- rep(FALSE, sum(dd$npars))
  fix[dd$npi %in% unusedpars(dd$mod)] <- TRUE
  sv <- drmTMBstartvals(dd, start, fix, method=method, gsl=gsl, gsu=gsu, gsn=gsn)
  names(sv) <- coefnames(dd$Xs)
  sv
}


# model data shared by drmTMB and drmTMBstart
drmTMBdata <- function(form, fform=NULL, rform=NULL, data, weights=NULL, family="gaussian", model="logistic", link="identity"){
  # response
  mf <- model.frame(form, data=data)
  # fixed effects
  fform <- fforms(fform)
  Xs <- lapply(fform, function(x) model.matrix(x, data=data))
  npars <- sapply(Xs, ncol)
  if (is.null(weights)) weights <- rep(1, nrow(data))

  # model
  mod <- which(model == c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal", "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz"))
  lnk <- which(link == c("identity", "log", "logit"))
  fam <- which(family == c("gaussian", "poisson", "binomial", "beta"))
  if (length(mod) != 1) stop("Unknown model: ", model)

  # response
  y <- mr <- model.response(mf, "any")
  if (family == "binomial" & ncol(as.matrix(mr)) == 2){
    bn <- mr[,1] + mr[,2]
    y <- mr[,1]
  } else if (family == "binomial"){
    # binary response: one trial per observation
    bn <- rep(1, length(y))
  } else {
    bn <- rep(0, length(y))
  }

  # random effects
  re <- reterms(rform, data)
  rpars <- unlist(lapply(re$terms, function(tm) tm$pars))
  unused <- intersect(unusedpars(mod), rpars)
  if (length(unused) > 0) stop("The ", model, " model does not use ", paste("b", unused, sep="", collapse=", "), "; remove it from rform.")

  list(y=y, x=mf[,2], w=weights, bn=bn, fam=fam, mod=mod, lnk=lnk, re=re,
       Xs=Xs, fform=fform, npars=npars, npi=rep(1:5, npars),
       family=family, model=model, link=link)
}


# starting values for all fixed-effect parameters; only NA entries of start are replaced
drmTMBstartvals <- function(dd, start, fix, method="lm", gsl=NULL, gsu=NULL, gsn=NULL){
  method <- match.arg(method, c("lm", "grid"))
  sv <- lmstart(dd, start, fix)
  if (method == "grid") sv <- gridstart(dd, sv, fix, gsl, gsu, gsn)
  sv
}


# starting values from linearised model functions
lmstart <- function(dd, start, fix){
  npi <- dd$npi
  mod <- dd$mod
  known1 <- fix[npi == 1][1] | !is.na(start[npi == 1][1])
  known2 <- fix[npi == 2][1] | !is.na(start[npi == 2][1])

  # response on the mean scale, aggregated by dose
  if (dd$family == "binomial" & any(dd$bn > 0)){
    yd <- tapply(dd$w * dd$y, dd$x, sum) / tapply(dd$w * dd$bn, dd$x, sum)
  } else {
    yd <- tapply(dd$w * dd$y, dd$x, sum) / tapply(dd$w, dd$x, sum)
  }
  xd <- sort(unique(dd$x))
  yd <- as.vector(yd)
  ylow <- yd[1]
  yhigh <- yd[length(yd)]

  # asymptotes on the mean scale: for the models with steepness b3, b1 is the lower and b2 the
  # upper asymptote; for michaelismenten, asymptoticregression and exponentialdecay, b1 is
  # approached at dose zero and b2 at high doses
  hasb3 <- !(mod %in% c(6, 7, 8))
  f1 <- if (hasb3) min(ylow, yhigh) else ylow
  f2 <- if (hasb3) max(ylow, yhigh) else yhigh
  # an asymptote with a given starting value keeps it; the other one is the more distant end of the curve
  if (known1) f1 <- linkfun(dd$lnk, start[npi == 1][1], inverse=TRUE)
  if (known2) f2 <- linkfun(dd$lnk, start[npi == 2][1], inverse=TRUE)
  if (known1 & !known2) f2 <- if (abs(ylow - f1) > abs(yhigh - f1)) ylow else yhigh
  if (known2 & !known1) f1 <- if (abs(ylow - f2) > abs(yhigh - f2)) ylow else yhigh
  # observed ends are not the asymptotes: extend the range by 5%
  if (!known1) f1 <- f1 + 0.05*(f1 - f2)
  if (!known2) f2 <- f2 - 0.05*(f1 - f2)
  f1 <- linkrange(dd$lnk, f1, dd$family)
  f2 <- linkrange(dd$lnk, f2, dd$family)

  # value of the increasing function F of the model; points at the asymptotes carry no information
  q <- (yd - f1) / (f2 - f1)
  ok <- is.finite(q) & q > 0.01 & q < 0.99
  if (mod %in% c(2, 3, 4, 5, 6)) ok <- ok & xd > 0
  q <- q[ok]
  xo <- xd[ok]

  # linearised model functions: inverse of F, linear in the dose (or log dose) with slope b3
  b3 <- b4 <- NA
  if (length(unique(xo)) >= 2){
    lf <- switch(mod,
                 lm(qlogis(q) ~ xo),                # logistic
                 lm(qlogis(q) ~ log(xo)),           # loglogistic
                 lm(I(-log(-log(q))) ~ log(xo)),    # weibull1
                 lm(log(-log(1 - q)) ~ log(xo)),    # weibull2
                 lm(qnorm(q) ~ log(xo)),            # lognormal
                 lm(I(1/q - 1) ~ 0 + I(1/xo)),      # michaelismenten
                 NULL,                              # asymptoticregression
                 NULL,                              # exponentialdecay
                 lm(I(-log(-log(q))) ~ xo))         # gompertz
    if (!is.null(lf)){
      cf <- coefficients(lf)
      if (mod %in% c(1, 9)){
        b3 <- cf[2]
        b4 <- -cf[1]/cf[2]
      }
      if (mod %in% c(2, 3, 4, 5)){
        b3 <- cf[2]
        b4 <- exp(-cf[1]/cf[2])
      }
      if (mod == 6) b4 <- cf[1]
    }
  }
  # fallback: unit slope in the direction of the curve (b3 > 0: increasing), dose at its midpoint
  if (!is.finite(b3)) b3 <- if (yhigh >= ylow) 1 else -1
  if (!is.finite(b4) | (mod %in% c(2, 3, 4, 5, 6, 7, 8) & b4 <= 0)){
    xpos <- xd[xd > 0]
    b4 <- xpos[which.min(abs((yd[xd > 0] - f1) / (f2 - f1) - 0.5))]
  }
  if (mod %in% c(6, 7, 8)) b3 <- 0

  # asymptoticregression, exponentialdecay: b2 is often not approached within the dose range;
  # for given b4 the model is linear in the asymptotes, so profile the residual sum of squares over b4
  if (mod %in% c(7, 8)){
    xp <- xd[xd > 0]
    cols <- c(!known1, !known2)
    b4grid <- exp(seq(log(min(xp)/10), log(max(xp)*10), length.out=100))
    pfit <- lapply(b4grid, function(b4){
      g <- 1 - exp(-xd/b4)
      off <- (if (known1) f1*(1 - g) else 0) + (if (known2) f2*g else 0)
      if (!any(cols)) return(list(rss=sum((yd - off)^2)))
      fit <- lm.fit(cbind(1 - g, g)[, cols, drop=FALSE], yd - off)
      list(rss=sum(fit$residuals^2), coef=fit$coefficients)
    })
    best <- which.min(sapply(pfit, function(p) p$rss))
    b4 <- b4grid[best]
    fa <- c(f1, f2)
    fa[cols] <- pfit[[best]]$coef
    f1 <- linkrange(dd$lnk, fa[1], dd$family)
    f2 <- linkrange(dd$lnk, fa[2], dd$family)
  }

  pooled <- c(linkfun(dd$lnk, f1), linkfun(dd$lnk, f2), b3, b4, 0)

  # project the pooled value onto each design matrix
  for (i in 1:5){
    id <- which(npi == i)
    X <- dd$Xs[[i]]
    cf <- qr.coef(qr(X), rep(pooled[i], nrow(X)))
    cf[is.na(cf) | abs(cf) < 1e-8 * max(abs(pooled[i]), 1)] <- 0
    fill <- is.na(start[id]) & !fix[id]
    start[id[fill]] <- cf[fill]
  }
  start[is.na(start)] <- 0
  as.vector(start)
}


# best point of a grid of parameter values, evaluated without random effects
gridstart <- function(dd, start, fix, gsl, gsu, gsn){
  np <- length(start)
  if (is.null(gsl)) gsl <- rep(NA, np)
  if (is.null(gsu)) gsu <- rep(NA, np)
  if (is.null(gsn)) gsn <- 10
  if (length(gsl) != np | length(gsu) != np) stop("gsl and gsu must have the same length as start, with NA for default ranges.")
  if (length(gsn) == 1) gsn <- rep(gsn, np)
  if (length(gsn) != np) stop("gsn must be a single number or have the same length as start.")
  gsn[is.na(gsn)] <- 10

  # grid values for each parameter: user bounds, otherwise default ranges
  dflt <- gridrange(dd, start)
  gvals <- lapply(1:np, function(i){
    if (fix[i] | gsn[i] < 2) return(NULL)
    if (!is.na(gsl[i]) & !is.na(gsu[i])) return(seq(gsl[i], gsu[i], length.out=gsn[i]))
    if (is.null(dflt[[i]])) return(NULL)
    r <- dflt[[i]]
    if (!is.na(gsl[i])) r$lower <- gsl[i]
    if (!is.na(gsu[i])) r$upper <- gsu[i]
    if (r$log) sign(r$lower) * exp(seq(log(abs(r$lower)), log(abs(r$upper)), length.out=gsn[i]))
    else seq(r$lower, r$upper, length.out=gsn[i])
  })
  gid <- which(!sapply(gvals, is.null))
  if (length(gid) == 0) return(start)

  obj <- drmTMBobj(dd, start, fix, rform=NULL)
  # position of the grid parameters in obj$par (fixed parameters are mapped out)
  ppos <- cumsum(!fix)[gid]
  grid <- expand.grid(gvals[gid])
  objvals <- apply(grid, 1, function(g){
    par <- obj$par
    par[ppos] <- g
    obj$fn(par)
  })
  start[gid] <- as.numeric(grid[which.min(objvals), ])
  start
}


# default grid ranges for b3, b4 and b5 from the data and the "lm" starting values;
# only for parameters with a single coefficient, and not for the asymptotes
gridrange <- function(dd, start){
  mod <- dd$mod
  npi <- dd$npi
  x <- dd$x
  xp <- x[x > 0]
  rng <- vector("list", length(start))
  single <- function(i) if (dd$npars[i] == 1) which(npi == i) else NULL

  # b3: keep the sign of the starting value, steepness between 1/4 and 4 times
  i3 <- single(3)
  if (!is.null(i3) & mod %in% c(1, 2, 3, 4, 5, 9)){
    b3 <- start[i3]
    if (b3 == 0) b3 <- 1
    rng[[i3]] <- list(lower=b3/4, upper=b3*4, log=TRUE)
  }
  # b4: location within the dose range, extended at both ends
  i4 <- single(4)
  if (!is.null(i4)){
    if (mod %in% c(1, 9)){
      ext <- 0.25*diff(range(x))
      rng[[i4]] <- list(lower=min(x) - ext, upper=max(x) + ext, log=FALSE)
    } else if (length(xp) > 0){
      rng[[i4]] <- list(lower=min(xp)/2, upper=max(xp)*2, log=TRUE)
    }
  }
  # b5: asymmetry on the log scale, only used by logistic and loglogistic
  i5 <- single(5)
  if (!is.null(i5) & mod %in% c(1, 2)){
    rng[[i5]] <- list(lower=-2, upper=2, log=FALSE)
  }
  rng
}


# TMB objective function for given starting values
drmTMBobj <- function(dd, start, fix, rform=NULL){
  npi <- dd$npi
  sb <- lapply(1:5, function(i) start[npi == i])
  re <- if (is.null(rform)) reterms(NULL, length(dd$y)) else dd$re
  tnp <- vapply(re$terms, function(tm) length(tm$pars), integer(1))
  tnl <- vapply(re$terms, function(tm) length(tm$levels), integer(1))
  # log standard deviations and transformed correlations for each term; random effects at zero
  stheta <- unlist(lapply(tnp, function(np) c(rep(0.001, np), rep(0, np*(np-1)/2))))
  up <- rep(0, sum(tnp*tnl))
  dlist <- list(model="drcTMB", y=dd$y, x=dd$x, w=dd$w, bn=dd$bn, fam=dd$fam, mod=dd$mod, lnk=dd$lnk,
                tnp=tnp, tnl=tnl, tpar=as.integer(unlist(lapply(re$terms, function(tm) tm$pars)) - 1),
                simre=1L, Z=re$Z)
  datalist <- c(dlist, dd$Xs)

  # log_sigma: residual sd (gaussian) or precision (beta); fixed for poisson and binomial
  if (dd$family == "gaussian") slsigma <- if (is.null(rform)) 1 else 0
  if (dd$family == "beta") slsigma <- 10
  if (dd$family %in% c("poisson", "binomial")) slsigma <- 0
  parameters <- list(b1=sb[[1]], b2=sb[[2]], b3=sb[[3]], b4=sb[[4]], b5=sb[[5]], u=up, theta=as.numeric(stheta), log_sigma=slsigma)

  # fix parameters to starting values
  lmap <- list()
  if (any(fix)){
    idp <- 1:length(start)
    idp[fix] <- NA
    for (i in 1:5){
      lmap[[paste("b", i, sep="")]] <- as.factor(idp[npi == i])
    }
  }
  if (dd$family %in% c("poisson", "binomial")) lmap$log_sigma <- factor(NA)

  MakeADFun(datalist, parameters, map=lmap, random=if (!is.null(rform)) "u", DLL="drcTMB_TMBExports", silent=TRUE)
}


# right-hand side formulas X1-X5 for b1-b5 from fform: a formula b2 + b4 ~ Stress, a list of
# such formulas, or a list of five right-hand side formulas
fforms <- function(fform){
  res <- rep(list(~ 1), 5)
  names(res) <- paste("X", 1:5, sep="")
  if (is.null(fform)) return(res)
  if (inherits(fform, "formula")) fform <- list(fform)
  if (!is.list(fform) || !all(sapply(fform, inherits, "formula"))) stop("fform must be a formula or a list of formulas.")
  twosided <- sapply(fform, length) == 3
  if (!any(twosided)){
    if (length(fform) != 5) stop("A list of right-hand side formulas in fform needs one formula for each of b1 to b5.")
    res[] <- fform
    return(res)
  }
  if (!all(twosided)) stop("Each formula in fform needs parameters on the left-hand side, e.g. b2 + b4 ~ Stress.")
  done <- integer(0)
  for (f in fform){
    pnames <- all.vars(f[[2]])
    if (!all(pnames %in% paste("b", 1:5, sep=""))) stop("The left-hand side of a formula in fform can only contain b1 to b5.")
    pars <- match(pnames, paste("b", 1:5, sep=""))
    if (any(pars %in% done)) stop("Each parameter can only appear once in fform.")
    done <- c(done, pars)
    rhs <- f
    rhs[2] <- NULL
    res[pars] <- list(rhs)
  }
  res
}

# coefficient names of the form b2:Stressdrought
coefnames <- function(Xs){
  unlist(lapply(1:5, function(i) paste("b", i, ":", colnames(Xs[[i]]), sep="")))
}

# random intercept terms from rform, a formula or a list of formulas of the form
# b1 + b4 ~ 1 | group; nested grouping factors (1 | a/b) give one term per level of nesting
reterms <- function(rform, data){
  if (is.null(rform)){
    n <- if (is.numeric(data)) data else nrow(data)
    return(list(terms=list(), Z=Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(n, 0))))
  }
  if (inherits(rform, "formula")) rform <- list(rform)
  terms <- list()
  Zl <- list()
  for (f in rform){
    if (length(f) != 3) stop("Each random effect formula needs parameters on the left-hand side, e.g. b1 + b4 ~ 1 | group.")
    pnames <- all.vars(f[[2]])
    if (!all(pnames %in% paste("b", 1:5, sep=""))) stop("The left-hand side of a random effect formula can only contain b1 to b5.")
    pars <- sort(match(pnames, paste("b", 1:5, sep="")))
    frh <- f
    frh[2] <- NULL
    rt <- mkReTrms(findbars(frh), model.frame(subbars(frh), data=data))
    if (!all(sapply(rt$cnms, function(cn) identical(cn, "(Intercept)")))) stop("Only random intercepts are supported, e.g. b4 ~ 1 | group.")
    for (j in seq_along(rt$Ztlist)){
      Zt <- rt$Ztlist[[j]]
      # name and formula identify the term again in new data (predict)
      terms[[length(terms) + 1]] <- list(pars=pars, name=names(rt$Ztlist)[j], group=sub("^1 \\| ", "", names(rt$Ztlist)[j]),
                                         levels=rownames(Zt), formula=frh)
      Zl[[length(Zl) + 1]] <- t(Zt)
    }
  }
  list(terms=terms, Z=do.call(cbind, Zl))
}


# dose-response mean function, same as drcmean in the TMB template: f1 + (f2 - f1) * F,
# with F increasing from 0 to 1; f1 is the lower and f2 the upper asymptote for the models with
# steepness f3, and the asymptotes at dose zero and at high doses for models 6 to 8
drcmeanR <- function(mod, x, f1, f2, f3, f4, f5){
  t <- f3 * (if (mod %in% c(2, 3, 4, 5)) log(x) - log(f4) else x - f4)
  Fx <- switch(mod,
               (1 + exp(-t))^(-exp(f5)),  # logistic
               (1 + exp(-t))^(-exp(f5)),  # loglogistic
               exp(-exp(-t)),             # weibull1
               1 - exp(-exp(t)),          # weibull2
               pnorm(t),                  # lognormal
               1 / (1 + f4 / x),          # michaelismenten
               1 - exp(-x / f4),          # asymptoticregression
               1 - exp(-x / f4),          # exponentialdecay
               exp(-exp(-t)))             # gompertz
  f1 + (f2 - f1) * Fx
}

# start or fix as a vector with one entry per coefficient: an unnamed vector is used as it is,
# a vector named by parameter (e.g. c(b1=0, b5=0)) gives the value for all coefficients of that parameter
parvector <- function(v, npi, default, arg){
  out <- rep(default, length(npi))
  if (is.null(v)) return(out)
  if (is.null(names(v))){
    if (length(v) != length(npi)) stop(arg, " must have one entry for each of the ", length(npi), " model coefficients, or be named by parameter, e.g. c(b5=0).")
    return(v)
  }
  pn <- match(names(v), paste("b", 1:5, sep=""))
  if (any(is.na(pn)) | any(duplicated(pn))) stop("the names of ", arg, " must be distinct parameter names b1-b5.")
  for (i in seq_along(v)) out[npi == pn[i]] <- v[[i]]
  out
}

# parameters not used by the model, which cannot be estimated: b3 without steepness
# (michaelismenten, asymptoticregression, exponentialdecay), b5 without asymmetry
unusedpars <- function(mod){
  c(if (mod %in% c(6, 7, 8)) 3, if (!(mod %in% c(1, 2))) 5)
}

# link function (inverse=FALSE) or inverse link for the asymptote parameters
linkfun <- function(lnk, x, inverse=FALSE){
  if (inverse){
    switch(lnk, x, exp(x), plogis(x))
  } else {
    switch(lnk, x, log(x), qlogis(x))
  }
}

# keep asymptotes inside the range of the link function, and inside (0, 1) for proportions
linkrange <- function(lnk, f, family){
  if (family %in% c("binomial", "beta") | lnk == 3) return(min(max(f, 1e-4), 1 - 1e-4))
  if (lnk == 2) return(max(f, 1e-8))
  f
}
