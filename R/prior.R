#' Prior distributions for drmTMB models
#'
#' Prior distributions for the \code{prior} argument of \code{\link{drmTMB}}. With priors, the
#' estimates are the posterior mode (maximum a posteriori), and the standard errors are from the
#' curvature of the log posterior, see \code{\link{drmTMB}}.
#'
#' \itemize{
#'   \item \code{normal(mean, sd)}: normal distribution.
#'   \item \code{student_t(df, location, scale)}: Student t distribution with df degrees of freedom,
#'     shifted by location and scaled by scale.
#'   \item \code{lognormal(meanlog, sdlog)}: lognormal distribution of a positive parameter, with
#'     log(parameter) normal with mean meanlog and standard deviation sdlog.
#'   \item \code{pert(min, mode, max, shape)}: the PERT distribution of EFSA (2022) for benchmark dose
#'     modelling, a beta distribution on (min, max) with mode \code{mode} and the shape parameters
#'     1 + shape (mode - min)/(max - min) and 1 + shape (max - mode)/(max - min); shape = 4 by default.
#'   \item \code{half_normal(scale)}, \code{half_student_t(df, scale)}: normal and Student t
#'     distributions folded at 0, for standard deviations.
#' }
#'
#' A prior whose support is wider than the range of the parameter, e.g. a normal prior for a parameter
#' bounded below by 0, is truncated to that range and renormalised; the normalisation does not change
#' the estimates.
#'
#' @param mean,sd mean and standard deviation of the normal distribution
#' @param df degrees of freedom of the Student t distribution
#' @param location,scale location and scale of the Student t distribution, or scale of the half
#'   normal and half Student t distributions
#' @param meanlog,sdlog mean and standard deviation of the log of a lognormal parameter
#' @param min,mode,max minimum, mode and maximum of the PERT distribution
#' @param shape the shape (weight of the mode) of the PERT distribution
#'
#' @return a drcprior object
#' @name priors
#' @aliases normal student_t lognormal pert half_normal half_student_t
NULL

# prior distribution: type as in the TMB template (see lpdfprior), hyperparameters, support, a central
# value for starting values, and the distribution function for the normalisation of truncated priors
newprior <- function(name, type, hyp, support, center, cdf){
  if (any(!is.finite(hyp))) stop("the parameters of the ", name, " prior must be finite numbers.")
  structure(list(name=name, type=type, hyp=hyp, support=support, center=center, cdf=cdf), class="drcprior")
}

#' @rdname priors
#' @export
normal <- function(mean=0, sd=1){
  if (!isTRUE(sd > 0)) stop("sd of the normal prior must be positive.")
  newprior(sprintf("normal(%g, %g)", mean, sd), 1L, c(mean, sd, 0), c(-Inf, Inf), mean,
           function(z) pnorm(z, mean, sd))
}

#' @rdname priors
#' @export
student_t <- function(df=3, location=0, scale=1){
  if (!isTRUE(df > 0) || !isTRUE(scale > 0)) stop("df and scale of the student_t prior must be positive.")
  newprior(sprintf("student_t(%g, %g, %g)", df, location, scale), 2L, c(location, scale, df), c(-Inf, Inf), location,
           function(z) pt((z - location)/scale, df))
}

#' @rdname priors
#' @export
lognormal <- function(meanlog=0, sdlog=1){
  if (!isTRUE(sdlog > 0)) stop("sdlog of the lognormal prior must be positive.")
  newprior(sprintf("lognormal(%g, %g)", meanlog, sdlog), 3L, c(meanlog, sdlog, 0), c(0, Inf), exp(meanlog),
           function(z) plnorm(z, meanlog, sdlog))
}

#' @rdname priors
#' @export
pert <- function(min, mode, max, shape=4){
  if (!(isTRUE(min < mode) && isTRUE(mode < max))) stop("the pert prior needs min < mode < max.")
  if (!isTRUE(shape > 0)) stop("shape of the pert prior must be positive.")
  a <- 1 + shape*(mode - min)/(max - min)
  b <- 1 + shape*(max - mode)/(max - min)
  p <- newprior(sprintf("pert(%g, %g, %g%s)", min, mode, max, if (shape == 4) "" else sprintf(", shape=%g", shape)), 4L,
                c(min, max, a, b), c(min, max), mode, function(z) pbeta((z - min)/(max - min), a, b))
  p
}

#' @rdname priors
#' @export
half_normal <- function(scale=1){
  if (!isTRUE(scale > 0)) stop("scale of the half_normal prior must be positive.")
  newprior(sprintf("half_normal(%g)", scale), 5L, c(scale, 0, 0), c(0, Inf), scale,
           function(z) pmax(2*pnorm(z, 0, scale) - 1, 0))
}

#' @rdname priors
#' @export
half_student_t <- function(df=3, scale=1){
  if (!isTRUE(df > 0) || !isTRUE(scale > 0)) stop("df and scale of the half_student_t prior must be positive.")
  newprior(sprintf("half_student_t(%g, %g)", df, scale), 6L, c(scale, df, 0), c(0, Inf), scale,
           function(z) pmax(2*pt(z/scale, df) - 1, 0))
}

#' @export
print.drcprior <- function(x, ...){
  cat(x$name, "\n")
  invisible(x)
}

# log density of a prior at z, as lpdfprior in the TMB template (without the normalisation)
lpdfprior <- function(p, z){
  h <- p$hyp
  switch(p$type,
         dnorm(z, h[1], h[2], log=TRUE),
         dt((z - h[1])/h[2], h[3], log=TRUE) - log(h[2]),
         dlnorm(z, h[1], h[2], log=TRUE),
         dbeta((z - h[1])/(h[2] - h[1]), h[3], h[4], log=TRUE) - log(h[2] - h[1]),
         log(2) + dnorm(z, 0, h[1], log=TRUE),
         log(2) + dt(z/h[1], h[2], log=TRUE) - log(h[1]))
}

# names of the variance parameters of a fit (without "log_"), as in coef: sd:b1|grp for the standard
# deviations of the random effects, and sigma, phi or theta for the residual parameter
priorvarnames <- function(terms, family){
  sdn <- unlist(lapply(terms, function(tm) paste("sd:", relabels(tm), "|", tm$group, sep="")))
  list(sd=if (is.null(sdn)) character(0) else sdn,
       resid=switch(family, gaussian=, lognormal="sigma", beta="phi", nbinom2="theta", nbinom1="phi", betabinomial="phi", NULL))
}

# priors for the TMB template from the prior argument of drmTMB: a named list of drcprior objects, named
# by parameter (b1-b5, with a single coefficient, on the scale of the parameter after the inverse link),
# by coefficient (e.g. "b4:(Intercept)", on the scale of the coefficient), "sd" (all standard deviations
# of the random effects) or "sd:b1|grp", and the residual parameter "sigma" (gaussian, lognormal), "phi"
# (beta, nbinom1, betabinomial) or "theta" (nbinom2), on their natural scale. Specific names take
# precedence over b1-b5 and "sd". lower and upper are the bounds of all coefficients (start).
#
# Returns the data of the template (prwhat: 0 coefficient, 1 theta, 2 log_sigma; pridx: index from 0
# in the coefficients b1-b5 (all, including fixed ones) or in theta; prtrans: 1 none, 2 exp, 3 inverse
# logit from the optimised value to the scale of the prior; prtype and prhyp, with the log of the
# probability of the range of the parameter in column 5), the bounds of the optimised values implied
# by the support of the priors, and a table of the priors for printing.
priordata <- function(prior, dd, fix, lower, upper, coefnames){
  empty <- list(prwhat=integer(0), pridx=integer(0), prtrans=integer(0), prtype=integer(0), prhyp=matrix(0, 0, 5))
  if (is.null(prior) || length(prior) == 0) return(list(data=empty, table=NULL))
  if (inherits(prior, "drcprior") || !is.list(prior) || is.null(names(prior)) || any(names(prior) == "")) {
    stop("prior must be a list of priors named by parameter, e.g. list(b1=normal(10, 2), sd=half_normal(1)).")
  }
  if (any(duplicated(names(prior)))) stop("the names of prior must be distinct.")
  if (!all(vapply(prior, inherits, logical(1), "drcprior"))) stop("the entries of prior must be priors, e.g. normal(0, 1) or pert(0, 5, 10); see ?priors.")
  npi <- dd$npi
  terms <- dd$re$terms
  vn <- priorvarnames(terms, dd$family)
  ntheta <- sum(vapply(terms, function(tm){
    np <- length(tm$pars)
    np + np*(np - 1)/2
  }, numeric(1)))
  # position of the log standard deviations in theta
  sdpos <- unlist(lapply(seq_along(terms), function(t){
    off <- sum(vapply(terms[seq_len(t - 1)], function(tm){
      np <- length(tm$pars)
      np + np*(np - 1)/2
    }, numeric(1)))
    off + seq_along(terms[[t]]$pars) - 1
  }))

  entries <- list()
  add <- function(key, what, idx, trans, p, range, label){
    entries[[key]] <<- list(what=what, idx=idx, trans=trans, p=p, range=range, label=label)
  }
  linkrng <- function(l) switch(l, c(-Inf, Inf), c(0, Inf), c(0, 1))
  parname <- paste("b", 1:5, sep="")
  for (nm in names(prior)){
    p <- prior[[nm]]
    if (nm %in% parname){
      k <- match(nm, parname)
      id <- which(npi == k)
      if (all(fix[id])) stop("the prior for ", nm, " is for a fixed or unused parameter.")
      if (length(id) != 1) stop("the prior for ", nm, " needs a single coefficient; give priors for its coefficients by name, e.g. \"", coefnames[id[1]], "\", on the scale of the link.")
      l <- dd$lnk[k]
      rng <- linkfun(l, c(lower[id], upper[id]), inverse=TRUE)
      rng <- c(max(rng[1], linkrng(l)[1]), min(rng[2], linkrng(l)[2]))
      if (!(nm %in% names(entries))) add(nm, 0L, id - 1L, l, p, rng, nm)
    } else if (nm %in% coefnames){
      id <- match(nm, coefnames)
      if (fix[id]) stop("the prior for ", nm, " is for a fixed coefficient.")
      for (key in names(entries)) if (entries[[key]]$what == 0 && entries[[key]]$idx == id - 1L) entries[[key]] <- NULL
      add(nm, 0L, id - 1L, 1L, p, c(lower[id], upper[id]), nm)
    } else if (nm == "sd" || nm %in% vn$sd){
      if (length(vn$sd) == 0) stop("the prior for ", nm, " needs random effects.")
      sel <- if (nm == "sd") seq_along(vn$sd) else match(nm, vn$sd)
      for (j in sel){
        if (nm == "sd" && vn$sd[j] %in% names(prior)) next
        add(vn$sd[j], 1L, as.integer(sdpos[j]), 2L, p, c(0, Inf), vn$sd[j])
      }
    } else if (identical(nm, vn$resid)){
      add(nm, 2L, 0L, 2L, p, c(0, Inf), nm)
    } else {
      stop("unknown name in prior: ", nm, "; use b1-b5, coefficient names (", paste(coefnames[!fix], collapse=", "), ")",
           if (length(vn$sd) > 0) paste(", sd, ", paste(vn$sd, collapse=", "), sep=""),
           if (!is.null(vn$resid)) paste(" or", vn$resid) else "", ".")
    }
  }

  # the range of each parameter within the support of its prior: bounds of the optimised values, and
  # the probability of the range for the normalisation of truncated priors
  ne <- length(entries)
  hyp <- matrix(0, ne, 5)
  bnd <- matrix(NA_real_, ne, 2)
  for (i in seq_len(ne)){
    e <- entries[[i]]
    r <- c(max(e$range[1], e$p$support[1]), min(e$range[2], e$p$support[2]))
    if (!(r[1] < r[2])) stop("the prior ", e$p$name, " of ", e$label, " has no support within the range of the parameter, (", e$range[1], ", ", e$range[2], ").")
    mass <- diff(e$p$cdf(r))
    if (!(mass > 0)) stop("the prior ", e$p$name, " of ", e$label, " has no probability within the range of the parameter, (", r[1], ", ", r[2], ").")
    hyp[i, ] <- c(e$p$hyp, rep(0, 4 - length(e$p$hyp)), log(mass))
    bnd[i, ] <- linkfun(e$trans, r)
    entries[[i]]$center <- linkfun(e$trans, min(max(e$p$center, r[1]), r[2]))
  }
  data <- list(prwhat=vapply(entries, function(e) e$what, integer(1)), pridx=vapply(entries, function(e) as.integer(e$idx), integer(1)),
               prtrans=vapply(entries, function(e) as.integer(e$trans), integer(1)), prtype=vapply(entries, function(e) e$p$type, integer(1)),
               prhyp=hyp)
  data <- lapply(data, unname)
  table <- data.frame(parameter=vapply(entries, function(e) e$label, character(1)),
                      prior=vapply(entries, function(e) e$p$name, character(1)),
                      scale=vapply(entries, function(e) if (e$what == 0 && !(e$label %in% parname)) "link" else "natural", character(1)),
                      stringsAsFactors=FALSE, row.names=NULL)
  list(data=data, lower=bnd[, 1], upper=bnd[, 2], center=vapply(entries, function(e) e$center, numeric(1)), table=table)
}

# the prior data for an objective without random effects (starting values): priors on theta left out
priorsubset <- function(pd, keep){
  lapply(pd, function(v) if (is.matrix(v)) v[keep, , drop=FALSE] else v[keep])
}

# position of each prior in the optimised parameter vector par of a fit with the free coefficients fix
# and the parameter names of obj$par
priorpos <- function(pdata, fix, parnames){
  free <- cumsum(!fix)
  ith <- which(parnames == "theta")
  ils <- which(parnames == "log_sigma")
  vapply(seq_along(pdata$prwhat), function(i){
    switch(pdata$prwhat[i] + 1, free[pdata$pridx[i] + 1], ith[pdata$pridx[i] + 1], ils)
  }, numeric(1))
}
