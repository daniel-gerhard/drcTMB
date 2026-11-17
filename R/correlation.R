#' Residual correlation structures and variance functions
#'
#' Correlated residuals within series (e.g. repeated measurements of a subject over time, or plots in a
#' field) and residual standard deviations that depend on a covariate or on the fitted mean, for the
#' gaussian and lognormal families of \code{\link{drmTMB}} (arguments \code{correlation} and
#' \code{variance}), as the corStruct and varFunc classes of nlme.
#'
#' The residuals e (y minus the mean; for the lognormal family log(y) minus the log median) of a series
#' are multivariate normal with covariance D R D: R is the correlation matrix of the structure, and D the
#' diagonal matrix of the residual standard deviations, sigma times the variance function g. Series are
#' independent, and observations without a series are one series.
#'
#' The formula of a correlation structure is \code{~ position | series}: the position of the observations
#' in time or space (one or more numeric variables; Euclidean distances), and an optional grouping
#' factor of the series (several factors with \code{a:b} or \code{a/b}). The correlation of two residuals
#' of a series at distance d:
#' \itemize{
#'   \item \code{cor_ar1}: rho^k for observations k positions apart, -1 < rho < 1; the position is an
#'     integer time (or, with \code{~ 1 | series}, the order of the observations in the data), as corAR1.
#'   \item \code{cor_car1}: phi^d for times d apart, 0 < phi < 1 (continuous-time AR1, corCAR1).
#'   \item \code{cor_cs}: rho for all pairs (compound symmetry, corCompSymm), -1/(m - 1) < rho < 1 for
#'     series of up to m observations; for rho > 0 the same as a random intercept of the series.
#'   \item \code{cor_exp}: exp(-d/range) (corExp); \code{cor_gaus}: exp(-(d/range)^2) (corGaus);
#'     \code{cor_spher}: 1 - 1.5 d/range + 0.5 (d/range)^3 for d < range, 0 beyond (corSpher). With
#'     \code{nugget=TRUE}, the correlation is (1 - nugget) times that, for d > 0.
#' }
#' Positions must be distinct within a series.
#'
#' Variance functions give the residual standard deviation sigma g(v) for a covariate v, the fitted
#' mean (the median for lognormal, \code{~ fitted}, the default) or a variable of the data (e.g.
#' \code{~ dose}):
#' \itemize{
#'   \item \code{var_power}: g = |v|^delta (varPower).
#'   \item \code{var_exp}: g = exp(delta v) (varExp).
#'   \item \code{var_ident}: g = exp(lambda_s) for the strata s of \code{~ 1 | stratum}, with lambda 0 for
#'     the first stratum (varIdent): a separate standard deviation for each stratum.
#' }
#' With the fitted mean as the covariate, the variance depends on the coefficients of the curve, and the
#' likelihood is maximised over both, unlike the iterative estimation of gnls in nlme, whose estimates
#' differ slightly.
#'
#' The parameters are estimated on an unconstrained scale (\code{coef} of the fit): for the correlation,
#' cor_par:rho (atanh(rho) for AR1, the logit of (rho + 1/(m - 1))/(1 + 1/(m - 1)) for compound symmetry),
#' cor_par:phi (logit(phi)), cor_par:range (log(range)) and cor_par:nugget (logit(nugget)); for the
#' variance functions var_par:delta (delta itself) and var_par:<stratum> (lambda, the log ratio of the
#' standard deviation of the stratum to that of the first). \code{summary} prints them on their natural
#' scale, and \code{confint} gives their intervals on that scale, named cor:rho, cor:phi, cor:range,
#' cor:nugget, var:delta and var:ratio:<stratum>. Priors (argument \code{prior} of \code{\link{drmTMB}})
#' are given by these names, on the natural scale, e.g. \code{prior=list("cor:phi"=pert(0, 0.5, 1))}; the
#' marginal likelihood (\code{\link{marglik}}) needs priors for all of them.
#'
#' @param form formula, see Details
#' @param nugget logical, include a nugget effect
#'
#' @return a correlation structure (class drccor) or a variance function (class drcvar)
#'
#' @examples
#' \dontrun{
#' drmTMB(y ~ dose, data=d, model="loglogistic", correlation=cor_car1(~ time | subject))
#' drmTMB(y ~ dose, data=d, model="loglogistic", variance=var_power())
#' }
#' @name correlation
#' @aliases cor_ar1 cor_car1 cor_cs cor_exp cor_gaus cor_spher var_power var_exp var_ident
NULL

newcor <- function(type, form, nugget=FALSE){
  if (!inherits(form, "formula") || length(form) != 2) stop("the formula of a correlation structure must be one-sided, e.g. ~ time | subject.")
  if (!is.logical(nugget) || length(nugget) != 1 || is.na(nugget)) stop("nugget must be TRUE or FALSE.")
  structure(list(type=type, form=form, nugget=nugget), class="drccor")
}

#' @rdname correlation
#' @export
cor_ar1 <- function(form=~ 1) newcor("ar1", form)
#' @rdname correlation
#' @export
cor_car1 <- function(form) newcor("car1", form)
#' @rdname correlation
#' @export
cor_cs <- function(form=~ 1) newcor("cs", form)
#' @rdname correlation
#' @export
cor_exp <- function(form, nugget=FALSE) newcor("exp", form, nugget)
#' @rdname correlation
#' @export
cor_gaus <- function(form, nugget=FALSE) newcor("gaus", form, nugget)
#' @rdname correlation
#' @export
cor_spher <- function(form, nugget=FALSE) newcor("spher", form, nugget)

newvar <- function(type, form){
  if (!inherits(form, "formula") || length(form) != 2) stop("the formula of a variance function must be one-sided, e.g. ~ fitted or ~ dose.")
  structure(list(type=type, form=form), class="drcvar")
}

#' @rdname correlation
#' @export
var_power <- function(form=~ fitted) newvar("power", form)
#' @rdname correlation
#' @export
var_exp <- function(form=~ fitted) newvar("exp", form)
#' @rdname correlation
#' @export
var_ident <- function(form) newvar("ident", form)

#' @export
print.drccor <- function(x, ...){
  cat("Correlation structure: ", cortext(x), "\n", sep="")
  invisible(x)
}

#' @export
print.drcvar <- function(x, ...){
  cat("Variance function: ", vartext(x), "\n", sep="")
  invisible(x)
}

cortext <- function(cs){
  paste(c(ar1="AR1", car1="CAR1", cs="compound symmetry", exp="exponential", gaus="Gaussian", spher="spherical")[[cs$type]],
        if (isTRUE(cs$nugget)) " with nugget", ", ", deparse1(cs$form), sep="")
}

vartext <- function(vf){
  paste(c(power="power", exp="exponential", ident="separate standard deviations")[[vf$type]], ", ", deparse1(vf$form), sep="")
}

# position (lhs) and grouping (rhs) parts of ~ position | group
splitform <- function(form){
  rhs <- form[[2]]
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) list(pos=rhs[[2]], grp=rhs[[3]]) else list(pos=rhs, grp=NULL)
}

# variables of the formulas of a correlation structure and a variance function, for the observations
# with missing values (the keyword fitted is not a variable)
corvars <- function(correlation, variance){
  v <- c(if (!is.null(correlation)) all.vars(correlation$form), if (!is.null(variance)) all.vars(variance$form))
  setdiff(v, "fitted")
}

# grouping factor from the expression g (a factor, a:b or a/b) in data; one group without g
groupfactor <- function(g, data){
  if (is.null(g)) return(factor(rep(1, nrow(data))))
  if (is.call(g) && (identical(g[[1]], as.name(":")) || identical(g[[1]], as.name("/")))){
    return(interaction(groupfactor(g[[2]], data), groupfactor(g[[3]], data), drop=TRUE, lex.order=TRUE))
  }
  factor(eval(g, data, environment()))
}

# positions of the observations: the numeric variables of the expression p (a + b for several), a
# matrix with a column for each
posmatrix <- function(p, data){
  vars <- all.vars(p)
  if (length(vars) == 0) return(NULL)
  terms <- if (is.call(p) && identical(p[[1]], as.name("+"))) as.list(attr(terms(as.formula(call("~", p))), "variables"))[-1] else list(p)
  pos <- sapply(terms, function(e) as.numeric(eval(e, data, environment())))
  matrix(pos, nrow=nrow(data))
}

# data of a correlation structure for the TMB template: type code, nugget, the observations of each
# series (indices from 0, ordered by series and position), the series lengths, the positions, the
# largest integer lag (ar1) and the starting values of the parameters
cordata <- function(correlation, data){
  n <- nrow(data)
  none <- list(ctype=0L, cnug=0L, cord=integer(0), clen=integer(0), cpos=matrix(0, 0, 1), cmaxk=0L)
  if (is.null(correlation)) return(c(none, list(cpar=numeric(0), cnames=character(0))))
  if (!inherits(correlation, "drccor")) stop("correlation must be a correlation structure, e.g. cor_car1(~ time | subject); see ?correlation.")
  sf <- splitform(correlation$form)
  grp <- groupfactor(sf$grp, data)
  pos <- posmatrix(sf$pos, data)
  type <- correlation$type
  if (type %in% c("ar1", "car1", "exp", "gaus", "spher") && is.null(pos) && type != "ar1"){
    stop("the ", cortext(correlation), " structure needs the positions of the observations, e.g. ~ time | subject.")
  }
  if (type %in% c("ar1", "car1") && !is.null(pos) && ncol(pos) > 1) stop("the ", type, " structure needs a single time variable.")
  if (type == "ar1"){
    # integer times, or the order of the observations within each series
    if (is.null(pos)) pos <- matrix(ave(seq_len(n), grp, FUN=seq_along), ncol=1)
    if (any(abs(pos - round(pos)) > 1e-8)) stop("the times of the AR1 structure must be integers; use cor_car1 for continuous times.")
    pos <- round(pos)
  }
  if (is.null(pos)) pos <- matrix(0, n, 1)
  if (anyNA(pos)) stop("the positions of the correlation structure have missing values.")
  o <- order(grp, pos[, 1])
  clen <- as.vector(table(grp))
  clen <- clen[clen > 0]
  # distinct positions within series (a correlation of 1 between residuals)
  if (type != "cs"){
    key <- paste(grp, apply(pos, 1, paste, collapse=","))
    if (anyDuplicated(key)) stop("the positions of the observations must be distinct within each series of the correlation structure.")
  }
  m <- max(clen)
  cmaxk <- if (type == "ar1") as.integer(max(unlist(tapply(pos[, 1], grp, function(p) diff(range(p)))), 0)) else 0L
  # starting values on the unconstrained scale: correlation 0.2, range a third of the median distance
  # within series, nugget 0.1
  if (type == "ar1") cpar <- c("rho"=atanh(0.2))
  if (type == "car1") cpar <- c("phi"=qlogis(0.2))
  if (type == "cs"){
    lo <- if (m > 1) -1/(m - 1) else -1
    cpar <- c("rho"=qlogis((0.2 - lo)/(1 - lo)))
  }
  if (type %in% c("exp", "gaus", "spher")){
    dl <- unlist(lapply(split(seq_len(n), grp), function(i) if (length(i) > 1) as.vector(dist(pos[i, , drop=FALSE]))))
    dl <- dl[dl > 0]
    if (length(dl) == 0) stop("the ", cortext(correlation), " structure needs series with at least two observations.")
    cpar <- c("range"=log(median(dl)/3))
    if (correlation$nugget) cpar <- c(cpar, "nugget"=qlogis(0.1))
  }
  if (m < 2) stop("the series of the correlation structure have one observation each; there is nothing to correlate.")
  list(ctype=match(type, c("ar1", "car1", "cs", "exp", "gaus", "spher")), cnug=as.integer(correlation$nugget),
       cord=as.integer(o - 1), clen=as.integer(clen), cpos=pos, cmaxk=cmaxk, cpar=cpar, cnames=paste("cor_par:", names(cpar), sep=""))
}

# data of a variance function for the TMB template: type code, whether the covariate is the fitted mean,
# the covariate, the strata (from 0) and the starting values of the parameters
vardata <- function(variance, data){
  n <- nrow(data)
  if (is.null(variance)) return(list(vtype=0L, vfit=0L, vcov=numeric(0), vstr=integer(0), vpar=numeric(0), vnames=character(0)))
  if (!inherits(variance, "drcvar")) stop("variance must be a variance function, e.g. var_power(); see ?correlation.")
  sf <- splitform(variance$form)
  type <- variance$type
  if (type == "ident"){
    if (is.null(sf$grp) || !identical(sf$pos, 1)) stop("var_ident needs the strata as ~ 1 | stratum.")
    s <- groupfactor(sf$grp, data)
    s <- droplevels(s)
    if (nlevels(s) < 2) stop("var_ident needs at least two strata.")
    vpar <- setNames(rep(0, nlevels(s) - 1), levels(s)[-1])
    return(list(vtype=3L, vfit=0L, vcov=numeric(n), vstr=as.integer(s) - 1L, vpar=vpar, vnames=paste("var_par:", names(vpar), sep=""),
                vlevels=levels(s)))
  }
  if (!is.null(sf$grp)) stop("the ", type, " variance function has no strata; use ~ covariate.")
  fit <- identical(sf$pos, as.name("fitted")) && !("fitted" %in% names(data))
  v <- if (fit) numeric(n) else as.numeric(eval(sf$pos, data, environment()))
  if (!fit && length(v) != n) stop("the covariate of the variance function must have a value for each observation.")
  if (!fit && type == "power" && any(v == 0)) stop("the covariate of the power variance function is 0 for some observations, which gives a standard deviation of 0; use var_exp or a shifted covariate.")
  list(vtype=match(type, c("power", "exp")), vfit=as.integer(fit), vcov=v, vstr=integer(n), vpar=c("delta"=0), vnames="var_par:delta")
}

# transformations of the parameters of the correlation structure and the variance function from their
# unconstrained scale (coef) to the natural scale, all increasing: a list named by parameter, with the
# natural name and the function
restrans <- function(object){
  cs <- object$correlation
  vf <- object$variance
  out <- list()
  if (!is.null(cs)){
    m <- object$cormax
    lo <- if (!is.null(m) && m > 1) -1/(m - 1) else -1
    if (cs$type == "ar1") out[["cor_par:rho"]] <- list(name="cor:rho", f=tanh)
    if (cs$type == "car1") out[["cor_par:phi"]] <- list(name="cor:phi", f=plogis)
    if (cs$type == "cs") out[["cor_par:rho"]] <- list(name="cor:rho", f=function(p) lo + (1 - lo)*plogis(p))
    if (cs$type %in% c("exp", "gaus", "spher")){
      out[["cor_par:range"]] <- list(name="cor:range", f=exp)
      if (cs$nugget) out[["cor_par:nugget"]] <- list(name="cor:nugget", f=plogis)
    }
  }
  if (!is.null(vf)){
    if (vf$type %in% c("power", "exp")) out[["var_par:delta"]] <- list(name="var:delta", f=identity)
    if (vf$type == "ident") for (s in object$varlevels[-1]) out[[paste("var_par:", s, sep="")]] <- list(name=paste("var:ratio:", s, sep=""), f=exp)
  }
  out
}

# the parameters of the correlation structure and the variance function on their natural scale
resnatural <- function(object){
  rt <- restrans(object)
  if (length(rt) == 0) return(NULL)
  est <- coef(object)[names(rt)]
  setNames(vapply(names(rt), function(nm) rt[[nm]]$f(est[[nm]]), numeric(1)), vapply(rt, function(r) r$name, character(1)))
}

# the correlation structure, variance function and their parameters on the natural scale, for printing
resinfo <- function(object) list(correlation=object$correlation, variance=object$variance, natural=resnatural(object))

printresstruct <- function(object, digits){
  nat <- object$natural
  if (is.null(nat)) return(invisible(NULL))
  if (!is.null(object$correlation)){
    k <- grepl("^cor:", names(nat))
    cat("\nResidual correlation: ", cortext(object$correlation), ": ", paste(sub("^cor:", "", names(nat)[k]), "=", format(nat[k], digits=digits), collapse=", "), "\n", sep="")
  }
  if (!is.null(object$variance)){
    k <- grepl("^var:", names(nat))
    cat(if (is.null(object$correlation)) "\n", "Variance function: ", vartext(object$variance), ": ", paste(sub("^var:", "", names(nat)[k]), "=", format(nat[k], digits=digits), collapse=", "), "\n", sep="")
  }
  invisible(NULL)
}

# variance function g of the residual standard deviation sigma g, at the means mu of the observations of
# the fit (newdata NULL) or of the rows of newdata, with the parameters vpar (by default the estimates)
resg <- function(object, mu, newdata=NULL, vpar=NULL){
  vf <- object$variance
  if (is.null(vf)) return(rep(1, length(mu)))
  env <- object$obj$env
  if (is.null(vpar)) vpar <- env$parList(object$opt$par)$vpar
  dat <- env$data
  sf <- splitform(vf$form)
  if (vf$type == "ident"){
    s <- if (is.null(newdata)) dat$vstr + 1 else match(as.character(groupfactor(sf$grp, newdata)), object$varlevels)
    if (anyNA(s)) stop("newdata has strata of the variance function that are not in the fitted data.")
    return(exp(c(0, vpar))[s])
  }
  v <- if (dat$vfit == 1) mu else if (is.null(newdata)) dat$vcov else as.numeric(eval(sf$pos, newdata, environment()))
  if (length(v) != length(mu)) stop("newdata needs the covariate of the variance function, ", deparse1(sf$pos), ".")
  if (vf$type == "power") abs(v)^vpar[1] else exp(vpar[1]*v)
}

# correlation matrix of the residuals of a series with positions pos (rows), as corvalue of the template
cormatrix <- function(object, pos, cpar=NULL){
  dat <- object$obj$env$data
  if (is.null(cpar)) cpar <- object$obj$env$parList(object$opt$par)$cpar
  m <- nrow(pos)
  d <- as.matrix(dist(pos))
  r <- switch(dat$ctype,
              tanh(cpar[1])^round(d),
              plogis(cpar[1])^d,
              {lo <- if (object$cormax > 1) -1/(object$cormax - 1) else -1; matrix(lo + (1 - lo)*plogis(cpar[1]), m, m)},
              exp(-d/exp(cpar[1])),
              exp(-(d/exp(cpar[1]))^2),
              {z <- d/exp(cpar[1]); ifelse(z < 1, 1 - 1.5*z + 0.5*z^3, 0)})
  if (dat$ctype >= 4 && dat$cnug == 1) r <- r*(1 - plogis(cpar[2]))
  diag(r) <- 1
  r
}

# normalized residuals: standardised residuals z (residuals divided by their standard deviations)
# decorrelated within each series by the inverse Cholesky factor of its correlation matrix
normalize <- function(object, z){
  dat <- object$obj$env$data
  if (dat$ctype == 0) return(z)
  off <- 0
  for (m in dat$clen){
    i <- dat$cord[off + seq_len(m)] + 1
    R <- chol(cormatrix(object, dat$cpos[i, , drop=FALSE]))
    z[i] <- backsolve(R, z[i], transpose=TRUE)
    off <- off + m
  }
  z
}

# parameters of the correlation structure and the variance function that can have priors, on their
# natural scale, for priordata: the prior name, the target in the TMB template (3 cpar, 4 vpar), the
# index, the transformation from the optimised value (see prtrans of the template) and the range
resprior <- function(dd){
  out <- list()
  cd <- dd$cor
  vd <- dd$var
  if (!is.null(cd) && cd$ctype > 0){
    m <- max(cd$clen)
    lo <- if (m > 1) -1/(m - 1) else -1
    for (j in seq_along(cd$cnames)){
      nm <- sub("^cor_par:", "", cd$cnames[j])
      tr <- switch(cd$ctype, list(5L, c(-1, 1)), list(3L, c(0, 1)), list(6L, c(lo, 1)), NULL, NULL, NULL)
      if (nm == "range") tr <- list(2L, c(0, Inf))
      if (nm == "nugget") tr <- list(3L, c(0, 1))
      out[[paste("cor:", nm, sep="")]] <- list(what=3L, idx=j - 1L, trans=tr[[1]], range=tr[[2]], lo=lo)
    }
  }
  if (!is.null(vd) && vd$vtype > 0){
    for (j in seq_along(vd$vnames)){
      nm <- sub("^var_par:", "", vd$vnames[j])
      key <- if (vd$vtype == 3) paste("var:ratio:", nm, sep="") else paste("var:", nm, sep="")
      out[[key]] <- if (vd$vtype == 3) list(what=4L, idx=j - 1L, trans=2L, range=c(0, Inf))
                    else list(what=4L, idx=j - 1L, trans=1L, range=c(-Inf, Inf))
    }
  }
  out
}

# transformation of a prior target from the optimised value (inverse) or to it: linkfun, and tanh (5)
# and the compound symmetry correlation in (lo, 1) (6)
priortrans <- function(trans, x, inverse=FALSE, lo=-1){
  if (trans == 5) return(if (inverse) tanh(x) else atanh(x))
  if (trans == 6) return(if (inverse) lo + (1 - lo)*plogis(x) else qlogis((x - lo)/(1 - lo)))
  linkfun(trans, x, inverse=inverse)
}
