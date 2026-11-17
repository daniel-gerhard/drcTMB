# BMD parameterisation of drmTMB (argument bmdpar): b4 is replaced by the log benchmark dose of the curve
# of the intercepts, and computed from it and the other parameters in the TMB template (bmdb4).

# check bmdpar, a list with the benchmark response bmr and its definition def (see bmd), for the model,
# family and parameters of dd; returns the list with the role of b4 for the template (see bmdb4)
bmdparcheck <- function(bmdpar, dd, fix){
  if (!is.list(bmdpar) || is.null(bmdpar$bmr) || is.null(bmdpar$def)) stop("bmdpar must be a list with the benchmark response bmr and its definition def, e.g. list(bmr=0.1, def=\"relative\").")
  bmr <- bmdpar$bmr
  def <- match.arg(bmdpar$def, c("excess", "additional", "relative", "sd"))
  if (!is.numeric(bmr) || length(bmr) != 1 || !is.finite(bmr) || bmr <= 0) stop("bmr of bmdpar must be a single positive number.")
  if (def %in% c("excess", "additional") && bmr >= 1) stop("the ", def, " risk bmr must be below 1.")
  if (def %in% c("excess", "additional") && !(dd$family %in% c("binomial", "betabinomial", "beta"))){
    stop("the ", def, " risk is for responses that are proportions (binomial, betabinomial or beta family); use def=\"relative\" or \"sd\".")
  }
  if (def == "sd" && !(dd$family %in% c("gaussian", "lognormal"))) stop("def=\"sd\" is only available for the gaussian and lognormal families.")
  if (def == "sd" && length(dd$re$terms) > 0) stop("def=\"sd\" is not available for models with random effects.")
  if (def == "sd" && isTRUE(dd$var$vtype > 0)) stop("def=\"sd\" is not available for models with a variance function.")
  m <- drcmodel(dd$mod)
  role <- if (m$ed %in% "polynomial") 4 else if (m$ed %in% "unbounded") 3 else if (m$ed %in% "closed" && m$b4 %in% "scale") 1
          else if (m$ed %in% "closed" && m$b4 %in% "threshold") 2 else NA
  if (is.na(role)) stop("the BMD parameterisation is not available for the ", dd$model, " model", if (m$b4 %in% "location") ": its response at dose zero depends on b4" else "", ".")
  if (m$b4 %in% "threshold" && smoothwidth(dd$smooth) > 0) stop("the BMD parameterisation is not available for the smoothed threshold of the nec model: its response at dose zero depends on b4.")
  npi <- dd$npi
  if (dd$npars[4] != 1 || fix[npi == 4]) stop("the BMD parameterisation needs b4 with a single, estimated coefficient.")
  if (any(dd$npars[-4] != 1)) stop("the BMD parameterisation needs b1-b3 and b5 without covariates in fform.")
  rpars <- unlist(lapply(dd$re$terms, function(tm) tm$pars))
  if (4 %in% rpars) stop("the BMD parameterisation is not available with random effects on b4.")
  # multistage of degree 2 (b5 fixed at 0): the quadratic exponential model of EFSA (2022), with b3 the
  # technical parameter log d and both coefficients of the polynomial computed from the BMD (role 5)
  if (role == 4 && all(fix[npi == 5])){
    if (fix[npi == 3]) stop("the BMD parameterisation of the multistage model of degree 2 needs an estimated b3.")
    if (3 %in% rpars) stop("the BMD parameterisation of the multistage model of degree 2 is not available with random effects on b3.")
    return(list(bmr=bmr, def=def, role=5, b4pos=FALSE))
  }
  # b4 bounded below by 0 by default (multistage coefficient, nec threshold, with the identity link): the
  # computed b4 is kept non-negative
  list(bmr=bmr, def=def, role=role, b4pos=4 %in% m$lower0 && dd$lnk[4] == 1)
}

# stop for methods that need b4 itself, for fits in the BMD parameterisation
nobmdpar <- function(x, what){
  if (!is.null(x$bmdpar)) stop(what, " is not available for fits in the BMD parameterisation (bmdpar); use the fit without bmdpar.", call.=FALSE)
}

# description of the BMD parameterisation of a fit, for print and summary
printbmdpar <- function(bp, est, digits){
  if (is.null(bp)) return(invisible(NULL))
  if (isTRUE(bp$role == 5)) cat("log_d is the technical parameter of the quadratic exponential model: the polynomial is b ((1 - c) x/dmax + c (x/dmax)^2) with c = exp(-exp(log_d)) and the highest dose dmax\n")
  cat("log_bmd is the log benchmark dose (", bp$def, ", bmr = ", format(bp$bmr), ", ", if (bp$inc > 0) "increasing" else "decreasing",
      " curve) of the curve of the intercepts: BMD = ", format(exp(est), digits=digits), "\n", sep="")
  invisible(NULL)
}
