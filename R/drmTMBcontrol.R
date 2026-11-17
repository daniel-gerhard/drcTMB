#' Set control parameters for drcTMB
#'
#' The optimiser minimises the negative log-likelihood (with random effects, its Laplace approximation)
#' over the fixed parameters, with the gradient from TMB, and for nlminb without random effects also
#' the Hessian. The same settings are used for the bootstrap fits of \code{\link{bootstrap}},
#' \code{\link{anova.drmTMB}} and \code{\link{bmdma}}.
#'
#' optim does not use the Hessian, and its methods other than "L-BFGS-B" do not allow bounds; models with
#' bounds (e.g. the default lower bound 0 of b4 for some models, see \code{\link{drmTMB}}, or user
#' bounds lower and upper) need method "L-BFGS-B". Parameter values at which the negative log-likelihood
#' is not finite (e.g. a negative b4 of the BMD parameterisation) are given a large finite value for optim.
#'
#' With restart=TRUE, the model is fitted again starting at the estimates, which restarts the
#' approximations of the optimiser (e.g. after a false convergence, a stop at a kink of the
#' likelihood, or in a flat likelihood of a variance), and the better fit is kept: the restart if its
#' negative log-likelihood is lower, or as low and it converged. This doubles the time of the fits at
#' least, also of the bootstrap fits.
#'
#' @param optimiser the optimiser: "nlminb" (default, \code{\link{nlminb}}) or "optim" (\code{\link{optim}})
#' @param method the method of optim: "L-BFGS-B" (default), "BFGS", "CG" or "Nelder-Mead"
#' @param iter.max maximum number of iterations of the optimiser (iter.max of nlminb, maxit of optim)
#' @param eval.max maximum number of evaluations of the objective function of nlminb
#' @param optcontrol a list of further control parameters of the optimiser, the argument control of
#'   \code{\link{nlminb}} (e.g. rel.tol, x.tol, trace) or \code{\link{optim}} (e.g. reltol, factr, pgtol,
#'   trace, parscale); its iter.max, eval.max or maxit replace those above. The tolerances of optim are
#'   factr=1e5 (L-BFGS-B) and reltol=1e-10 (the other methods) unless given, close to rel.tol=1e-10 of nlminb.
#' @param inner a list of control parameters of the inner optimisation of the random effects for the
#'   Laplace approximation, see the argument inner.control of \code{\link[TMB]{MakeADFun}}, e.g.
#'   list(maxit=2000, tol10=1e-6)
#' @param restart if TRUE, the fit is restarted once at the estimates and the better fit is kept, see
#'   Details
#' @param smallsample if TRUE (default), gaussian and lognormal models without random effects use the unbiased
#'   residual variance with n - p residual degrees of freedom and t quantiles, as lm, nls and drc; if
#'   FALSE, the maximum likelihood estimate and normal quantiles. See \code{\link{drmTMB}}.
#'
#' @return a list with control parameters
#' @export
#'
#' @examples
#' drmTMBcontrol(optimiser="optim", optcontrol=list(factr=1e5), restart=TRUE)

drmTMBcontrol <- function(optimiser=c("nlminb", "optim"), method=c("L-BFGS-B", "BFGS", "CG", "Nelder-Mead"),
                          iter.max=1000, eval.max=2000, optcontrol=list(), inner=list(), restart=FALSE, smallsample=TRUE){
  optimiser <- match.arg(optimiser)
  method <- match.arg(method)
  if (!is.list(optcontrol) || (length(optcontrol) > 0 && is.null(names(optcontrol)))) stop("optcontrol must be a named list.")
  known <- if (optimiser == "nlminb"){
    c("eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol", "scale.init", "diff.g")
  } else {
    c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT",
      "warn.1d.NelderMead", "type", "lmm", "factr", "pgtol")
  }
  bad <- setdiff(names(optcontrol), known)
  if (length(bad) > 0) stop("unknown control parameters of ", optimiser, " in optcontrol: ", paste(bad, collapse=", "), ".")
  if (!is.null(optcontrol$fnscale)) stop("fnscale of optim cannot be changed: the objective is minimised.")
  if (!is.list(inner) || (length(inner) > 0 && is.null(names(inner)))) stop("inner must be a named list.")
  if (!is.logical(restart) || length(restart) != 1 || is.na(restart)) stop("restart must be TRUE or FALSE.")
  list(optimiser=optimiser, method=method, iter.max=iter.max, eval.max=eval.max, optcontrol=optcontrol, inner=inner,
       restart=restart, smallsample=smallsample)
}

# control parameters of a fit, completed with the defaults (fits from older versions of drcTMB)
fitcontrol <- function(fit){
  ctl <- fit$control
  do.call(drmTMBcontrol, ctl[intersect(names(ctl), names(formals(drmTMBcontrol)))])
}

# minimum of the objective of the TMB object obj over the fixed parameters from start, within the
# bounds lower and upper, by the optimiser of control; with the Hessian he (nlminb only, without
# random effects); restarted at the estimates for control$restart. The result has the elements of
# nlminb: par, objective, convergence (0 for success), message, iterations, evaluations
fitopt <- function(obj, start, lower, upper, control, he=NULL){
  opt <- runopt(obj, start, lower, upper, control, he)
  if (!isTRUE(control$restart) || !all(is.finite(opt$par))) return(opt)
  opt2 <- tryCatch(suppressWarnings(runopt(obj, opt$par, lower, upper, control, he)), error=function(e) NULL)
  better <- !is.null(opt2) && is.finite(opt2$objective) &&
    (!is.finite(opt$objective) || opt2$objective < opt$objective - 1e-8 ||
     (opt2$objective <= opt$objective + 1e-8 && opt2$convergence == 0))
  res <- if (better) opt2 else opt
  res$restart <- c(first=opt$objective, restart=if (is.null(opt2)) NA else opt2$objective)
  res
}

# one run of the optimiser of control, see fitopt
runopt <- function(obj, start, lower, upper, control, he=NULL){
  oc <- control$optcontrol
  if (control$optimiser == "nlminb"){
    ctrl <- c(list(iter.max=control$iter.max, eval.max=control$eval.max), oc)
    ctrl <- ctrl[!duplicated(names(ctrl), fromLast=TRUE)]
    return(nlminb(start=start, objective=obj$fn, gradient=obj$gr, hessian=he, lower=lower, upper=upper, control=ctrl))
  }
  method <- control$method
  bounded <- any(is.finite(lower)) || any(is.finite(upper))
  if (bounded && method != "L-BFGS-B") stop("the parameters have bounds; optim needs method=\"L-BFGS-B\" (see drmTMBcontrol).")
  # optim stops at non-finite values: a large value instead, with a zero gradient
  f0 <- obj$fn(start)
  big <- if (is.finite(f0)) abs(f0)*1e3 + 1e10 else 1e300
  fn <- function(p){
    v <- obj$fn(p)
    if (is.finite(v)) v else big
  }
  gr <- function(p){
    g <- obj$gr(p)
    if (all(is.finite(g))) g else numeric(length(p))
  }
  # tolerances close to the relative tolerance 1e-10 of nlminb, unless given
  tol <- if (method == "L-BFGS-B") list(factr=1e5) else list(reltol=1e-10)
  ctrl <- c(list(maxit=control$iter.max), tol, oc)
  ctrl <- ctrl[!duplicated(names(ctrl), fromLast=TRUE)]
  o <- if (bounded) optim(start, fn, gr, method=method, lower=lower, upper=upper, control=ctrl)
       else optim(start, fn, gr, method=method, control=ctrl)
  msg <- if (!is.null(o$message) && nzchar(o$message)) o$message
         else switch(as.character(o$convergence), "0"="converged", "1"="iteration limit maxit reached",
                     "10"="degenerate Nelder-Mead simplex", paste("convergence code", o$convergence))
  value <- if (o$value >= big) obj$fn(o$par) else o$value
  list(par=o$par, objective=value, convergence=o$convergence, message=msg,
       iterations=o$counts[[if (method == "Nelder-Mead") 1 else 2]], evaluations=o$counts)
}

# the inner control of the Laplace approximation of control$inner for a TMB object with random effects
setinner <- function(obj, control){
  inner <- control$inner
  if (is.null(obj$env$random) || length(inner) == 0) return(invisible(obj))
  # arguments of TMB's inner optimiser newton
  known <- c("trace", "maxit", "tol", "alpha", "smartsearch", "mgcmax", "super", "silent", "ustep", "power", "u0",
             "grad.tol", "step.tol", "tol10")
  bad <- setdiff(names(inner), known)
  if (length(bad) > 0) stop("unknown inner control parameters: ", paste(bad, collapse=", "), "; available: ",
                            paste(known, collapse=", "), ".")
  obj$env$inner.control[names(inner)] <- inner
  invisible(obj)
}

# name of the optimiser of a fit for summary: nlminb, or optim with its method
optlabel <- function(fit){
  ctl <- fitcontrol(fit)
  if (ctl$optimiser == "nlminb") "nlminb" else paste("optim (", ctl$method, ")", sep="")
}
