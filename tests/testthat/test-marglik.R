# marginal likelihood: Laplace approximation and importance sampling, checked against exact integrals
# (one and two dimensions) and adaptive Gauss-Hermite quadrature (nonlinear models, random effects)

library(drcData)
data(ryegrass)

# adaptive Gauss-Hermite quadrature of exp(-objective) over the optimised parameters, with m nodes per
# dimension, centred at the posterior mode and scaled by the covariance of the fit
aghq <- function(fit, m){
  J <- matrix(0, m, m)
  for (i in seq_len(m - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i/2)
  e <- eigen(J, symmetric=TRUE)
  x <- e$values
  w <- sqrt(pi)*e$vectors[1, ]^2
  par <- fit$opt$par
  k <- length(par)
  R <- chol(fit$sdrl$cov.fixed)
  idx <- as.matrix(expand.grid(rep(list(seq_len(m)), k)))
  X <- matrix(x[idx], ncol=k)
  lw <- rowSums(matrix(log(w[idx]), ncol=k)) + rowSums(X^2)
  th <- sweep(sqrt(2) * X %*% R, 2, par, "+")
  v <- apply(th, 1, function(p) -as.numeric(fit$obj$fn(p))) + lw
  mx <- max(v[is.finite(v)])
  mx + log(sum(exp(v - mx))) + k/2*log(2) + sum(log(diag(R)))
}

# log of an integral of exp(f) by integrate, shifted by a constant for numerical stability
logintegrate <- function(f, lo, hi, shift){
  log(integrate(function(z) exp(f(z) + shift), lo, hi, rel.tol=1e-10)$value) - shift
}

test_that("one parameter: Laplace and importance sampling match the exact integral", {
  set.seed(1)
  d <- data.frame(x=rep(0:4, each=3))
  d$y <- rpois(nrow(d), 4)
  f <- drmTMB(y ~ x, data=d, family="poisson", model="constant", parlink=c(b1="log"), prior=list(b1=lognormal(log(3), 0.5)))
  ex <- logintegrate(function(l) sapply(l, function(m) sum(dpois(d$y, m, log=TRUE)) + dlnorm(m, log(3), 0.5, log=TRUE)), 0, 50, 30)
  ml <- marglik(f)
  expect_s3_class(ml, "drcmarglik")
  expect_equal(attr(ml, "k"), 1)
  expect_equal(as.numeric(ml), ex, tolerance=0.01/abs(ex))
  set.seed(2)
  mi <- marglik(f, "is", n=2000)
  expect_lt(abs(mi - ex), 4*attr(mi, "se"))
  expect_gt(attr(mi, "ess"), 1000)
  expect_equal(attr(mi, "laplace"), as.numeric(ml))
  expect_output(print(ml), "Laplace approximation, 1 parameter")
  expect_output(print(mi), "Importance sampling, 2000 draws")
})

test_that("two parameters: the gaussian constant model with priors for the mean and sigma", {
  set.seed(3)
  g <- data.frame(x=rep(0:4, each=2))
  g$y <- rnorm(10, 5, 1)
  fg <- drmTMB(y ~ x, data=g, model="constant", prior=list(b1=normal(4, 2), sigma=half_normal(2)))
  inner <- function(s) sapply(s, function(si){
    exp(logintegrate(function(m) sapply(m, function(mm) sum(dnorm(g$y, mm, si, log=TRUE)) + dnorm(mm, 4, 2, log=TRUE)), -10, 20, 15)) * 2*dnorm(si, 0, 2)
  })
  ex <- log(integrate(inner, 0, 20, rel.tol=1e-10)$value)
  set.seed(4)
  mi <- marglik(fg, "is", n=4000)
  expect_lt(abs(mi - ex), 4*attr(mi, "se") + 0.002)
  # the Laplace approximation is close for 10 observations
  expect_lt(abs(marglik(fg) - ex), 0.15)
  expect_equal(aghq(fg, 15), ex, tolerance=1e-3/abs(ex))
})

test_that("nonlinear model: importance sampling matches adaptive quadrature", {
  fr <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), parlink=c(b4="log"),
               prior=list(b1=normal(0, 1), b2=normal(8, 2), b3=normal(-2, 2), b4=lognormal(log(3), 1), sigma=lognormal(log(0.5), 1)))
  # quadrature with 7 nodes is within 0.01 of 15 nodes for this fit
  ref <- aghq(fr, 7)
  set.seed(5)
  mi <- marglik(fr, "is", n=4000)
  expect_lt(abs(mi - ref), 4*attr(mi, "se") + 0.01)
  ml <- marglik(fr)
  expect_lt(abs(ml - ref), 0.5)
  # the Laplace approximation from its definition
  k <- length(fr$opt$par)
  expect_equal(as.numeric(ml), -fr$opt$objective + k/2*log(2*pi) + as.numeric(determinant(fr$sdrl$cov.fixed)$modulus)/2)
})

test_that("random effects are integrated out by the nested Laplace approximation", {
  set.seed(6)
  d <- expand.grid(dose=c(0, 0.5, 1, 2, 4, 8, 16), grp=factor(1:6))
  u <- rnorm(6, 0, 0.8)
  d$y <- 5 + u[d$grp] + rnorm(nrow(d), 0, 0.5)
  # gaussian random intercepts: the inner Laplace approximation is exact
  fc <- drmTMB(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="constant",
               prior=list(b1=normal(5, 3), sd=half_normal(1), sigma=half_normal(1)))
  ref <- aghq(fc, 11)
  set.seed(7)
  mi <- marglik(fc, "is", n=1000)
  expect_lt(abs(mi - ref), 4*attr(mi, "se") + 0.01)
  expect_lt(abs(marglik(fc) - ref), 0.2)
  # importance sampling keeps the state of the TMB object
  lpb <- fc$obj$env$last.par.best
  cf <- coef(fc)
  marglik(fc, "is", n=20)
  expect_identical(fc$obj$env$last.par.best, lpb)
  expect_identical(coef(fc), cf)

  # a nonlinear model with random effects on b2
  d$y <- 1 + (7 + u[d$grp] - 1)/(1 + (d$dose/3)^1.5) + rnorm(nrow(d), 0, 0.3)
  fl <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE),
               prior=list(b1=normal(0, 2), b2=normal(8, 3), b3=normal(-2, 2), b4=lognormal(log(3), 1),
                          sd=half_normal(1), sigma=half_normal(1)))
  set.seed(8)
  mi <- marglik(fl, "is", n=1000)
  expect_gt(attr(mi, "ess"), 300)
  expect_lt(abs(mi - attr(mi, "laplace")), 0.5)
})

test_that("missing priors, correlations and bounds", {
  f <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), prior=list(b2=normal(8, 2)))
  expect_error(marglik(f), "no prior for b1:\\(Intercept\\), b3:\\(Intercept\\), b4:\\(Intercept\\), sigma")
  f0 <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  expect_error(marglik(f0), "proper prior")
  expect_error(marglik(f, method="bridge"))

  set.seed(9)
  d <- expand.grid(dose=c(0, 0.5, 1, 2, 4, 8, 16), grp=factor(1:6))
  d$y <- 1 + (7 + rnorm(6, 0, 0.8)[d$grp] - 1)/(1 + (d$dose/3)^1.5) + rnorm(nrow(d), 0, 0.3)
  pr <- list(b1=normal(0, 2), b2=normal(8, 3), b3=normal(-2, 2), b4=lognormal(log(3), 1), sigma=half_normal(1))
  fs <- suppressWarnings(drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), prior=pr))
  expect_error(marglik(fs), "no prior for sd:b2\\|grp")
  fcor <- suppressWarnings(drmTMB(y ~ dose, rform=b1 + b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE),
                                  prior=c(pr, list(sd=half_normal(1)))))
  expect_error(marglik(fcor), "correlations of random effects have no priors")

  # a pert prior whose lower bound is close to the estimate of b4
  fb <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE),
               prior=list(b1=normal(0, 1), b2=normal(8, 2), b3=normal(-2, 2), b4=pert(2.95, 3.5, 10), sigma=half_normal(1)))
  expect_warning(marglik(fb), "within 2 standard errors of a bound")
  set.seed(10)
  mi <- marglik(fb, "is", n=500)
  expect_true(is.finite(mi))
})

test_that("Laplace and importance sampling agree with tmbstan and bridge sampling", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("DRCTMB_SLOW_TESTS"), "true"), "set DRCTMB_SLOW_TESTS=true to run")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("bridgesampling")
  # log marginal likelihood by bridge sampling from tmbstan draws of exp(-objective); with random
  # effects, the joint density of the fixed and random effects (MCMC over the random effects too)
  bridge <- function(fit){
    obj <- fit$obj
    r <- obj$env$random
    init <- if (is.null(r)) fit$opt$par else obj$env$last.par.best
    st <- suppressWarnings(tmbstan::tmbstan(obj, chains=4, iter=4000, warmup=1000, lower=fit$lower, upper=fit$upper,
                                            init=function() as.numeric(init), refresh=0, seed=1, silent=TRUE))
    X <- as.matrix(st)
    X <- X[, colnames(X) != "lp__", drop=FALSE]
    colnames(X) <- paste("p", seq_len(ncol(X)), sep="")
    lo <- hi <- NULL
    if (is.null(r)){
      lo <- fit$lower
      hi <- fit$upper
    } else {
      lo <- rep(-Inf, ncol(X))
      hi <- rep(Inf, ncol(X))
      lo[-r] <- fit$lower
      hi[-r] <- fit$upper
    }
    lpost <- function(x, data) -as.numeric(if (is.null(r)) obj$fn(x) else obj$env$f(x, order=0))
    b <- bridgesampling::bridge_sampler(X, log_posterior=lpost, data=NULL, lb=setNames(lo, colnames(X)),
                                        ub=setNames(hi, colnames(X)), silent=TRUE)
    b$logml
  }
  pr <- list(bmd=pert(0.1, 1, 5), b2=normal(8, 1), b3=normal(-2, 2), sigma=half_normal(1))
  f <- drmTMB(rootl ~ conc, data=ryegrass, model="weibull1", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
              bmdpar=list(bmr=0.05, def="relative"), prior=pr)
  ref <- bridge(f)
  set.seed(11)
  mi <- marglik(f, "is", n=10000)
  expect_lt(abs(mi - ref), 4*attr(mi, "se") + 0.02)
  expect_lt(abs(marglik(f) - ref), 0.5)

  # random effects on b2: the nested Laplace approximation integrates out the random effects
  set.seed(6)
  d <- expand.grid(dose=c(0, 0.5, 1, 2, 4, 8, 16), grp=factor(1:6))
  u <- rnorm(6, 0, 0.8)
  d$y <- 1 + (7 + u[d$grp] - 1)/(1 + (d$dose/3)^1.5) + rnorm(nrow(d), 0, 0.3)
  fl <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE),
               prior=list(b1=normal(0, 2), b2=normal(8, 3), b3=normal(-2, 2), b4=lognormal(log(3), 1),
                          sd=half_normal(1), sigma=half_normal(1)))
  ref <- bridge(fl)
  set.seed(12)
  mi <- marglik(fl, "is", n=5000)
  expect_lt(abs(mi - ref), 4*attr(mi, "se") + 0.02)
  expect_lt(abs(marglik(fl) - ref), 0.5)
})
