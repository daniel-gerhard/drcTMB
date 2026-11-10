# Gaussian models with random effects compared with nlme (maximum likelihood).
# With random effects on the asymptotes b1 and b2 (identity link), the model is linear in the
# random effects, so the Laplace approximation of drcTMB and the linearisation of nlme are both
# exact: the drcTMB objective at the nlme estimates equals the nlme log-likelihood.
# With a random effect on b4, both are approximations and only agree approximately.

skip_if_not_installed("nlme")

ll4 <- function(x, b1, b2, b3, b4) b1 + (b2 - b1)/(1 + exp(-b3*(log(x) - log(b4))))

simnlme <- function(re=c("b2", "b4"), seed=3){
  re <- match.arg(re)
  d <- simgroups(ng=10, nrep=2, seed=seed)
  u <- rnorm(10, 0, if (re == "b2") 0.4 else 0.3)
  # nlme cannot evaluate log(0); a tiny positive dose has the same mean as dose zero
  d$dose <- pmax(d$dose, 1e-8)
  b2 <- 8 + if (re == "b2") u[d$grp] else 0
  b4 <- 3 * if (re == "b4") exp(u[d$grp]) else 1
  d$y <- ll4(d$dose, 2, b2, -1.5, b4) + rnorm(nrow(d), 0, 0.3)
  d
}

# nlme estimates in the order of the drcTMB parameter vector (b1-b4, theta, log_sigma)
nlmepar <- function(fit){
  c(nlme::fixef(fit), log(as.numeric(nlme::VarCorr(fit)[1, 2])), log(fit$sigma))
}

nlmefit <- function(d, f1, random){
  nlme::nlme(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))), fixed=b1 + b2 + b3 + b4 ~ 1, random=random, groups=~grp,
             data=d, start=setNames(bcoef(f1), paste("b", 1:4, sep="")), method="ML")
}

test_that("random asymptote: drcTMB and nlme maximise the same exact likelihood", {
  d <- simnlme("b2")
  f1 <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  n1 <- nlmefit(d, f1, b2 ~ 1)

  # same likelihood function
  expect_equal(as.numeric(f1$obj$fn(nlmepar(n1))), -as.numeric(logLik(n1)), tolerance=1e-6)
  # drcTMB reaches at least the nlme optimum
  expect_lte(f1$opt$objective, -as.numeric(logLik(n1)) + 1e-6)
  expect_equal(unname(coef(f1)), unname(nlmepar(n1)), tolerance=1e-3)
})

test_that("correlated random asymptotes agree with nlme", {
  d <- simnlme("b2")
  f1 <- fit4(y ~ dose, rform=b1 + b2 ~ 1 | grp, data=d, model="loglogistic")
  n1 <- nlmefit(d, f1, nlme::pdSymm(b1 + b2 ~ 1))
  expect_lte(f1$opt$objective, -as.numeric(logLik(n1)) + 1e-6)
  expect_equal(f1$opt$objective, -as.numeric(logLik(n1)), tolerance=1e-4)
  expect_equal(unname(bcoef(f1)), unname(nlme::fixef(n1)), tolerance=1e-2)
})

test_that("random ED50 agrees approximately with nlme", {
  d <- simnlme("b4")
  f1 <- fit4(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="loglogistic")
  n1 <- nlmefit(d, f1, b4 ~ 1)
  # Laplace approximation vs linearisation: estimates within a fraction of a standard error
  se <- sqrt(diag(vcov(f1)))[1:4]
  expect_true(all(abs(bcoef(f1) - nlme::fixef(n1)) < 0.5*se))
  expect_lt(abs(f1$opt$objective + as.numeric(logLik(n1))), 0.5)
})

test_that("independent random asymptotes agree with nlme", {
  set.seed(21)
  d <- simgroups(ng=12, nrep=2, seed=21)
  u1 <- rnorm(12, 0, 0.3)
  u2 <- rnorm(12, 0, 0.6)
  d$dose <- pmax(d$dose, 1e-8)
  d$y <- ll4(d$dose, 2 + u1[d$grp], 8 + u2[d$grp], -1.5, 3) + rnorm(nrow(d), 0, 0.3)
  f1 <- fit4(y ~ dose, rform=list(b1 ~ 1 | grp, b2 ~ 1 | grp), data=d, model="loglogistic")
  n1 <- nlmefit(d, f1, nlme::pdDiag(b1 + b2 ~ 1))
  # same exact likelihood: one log sd per term, then log sigma
  p <- c(nlme::fixef(n1), log(as.numeric(nlme::VarCorr(n1)[, "StdDev"])))
  expect_equal(as.numeric(f1$obj$fn(p)), -as.numeric(logLik(n1)), tolerance=1e-6)
  expect_lte(f1$opt$objective, -as.numeric(logLik(n1)) + 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f1)))[1:4]), unname(sqrt(diag(vcov(n1)))), tolerance=2e-2)
  expect_equal(exp(unname(coef(f1)[5:7])), as.numeric(nlme::VarCorr(n1)[, "StdDev"]), tolerance=1e-2)
})

test_that("the variances start on the scale of the data (vinclozolin, response in thousands)", {
  skip_if_not_installed("drcData")
  data(vinclozolin, package="drcData")
  v <- vinclozolin
  v$dose <- pmax(v$conc, 1e-8)
  v$y <- v$effect
  # starting at standard deviations of 1, the fit stopped at the fit without random effects
  # (-logLik 390.68), with the between-experiment variation in the residual sd
  for (corr in c(TRUE, FALSE)){
    rf <- if (corr) b1 + b2 ~ 1 | exper else list(b1 ~ 1 | exper, b2 ~ 1 | exper)
    # NA/NaN evaluations on the way, e.g. at negative b4
    f1 <- suppressWarnings(fit4(y ~ dose, rform=rf, data=v, model="loglogistic"))
    n1 <- nlme::nlme(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))), fixed=b1 + b2 + b3 + b4 ~ 1,
                     random=if (corr) nlme::pdSymm(b1 + b2 ~ 1) else nlme::pdDiag(b1 + b2 ~ 1), groups=~exper, data=v,
                     start=c(b1=200, b2=2000, b3=-0.6, b4=0.05), method="ML")
    expect_lte(f1$opt$objective, -as.numeric(logLik(n1)) + 1e-4)
    expect_equal(f1$opt$objective, -as.numeric(logLik(n1)), tolerance=1e-3)
    expect_true(all(exp(coef(f1)[grepl("log_sd", names(coef(f1)))]) > 50))
  }
})

test_that("a variance stuck at zero is refitted", {
  set.seed(21)
  d <- simgroups(ng=12, nrep=2, seed=21)
  u2 <- rnorm(12, 0, 0.6)
  d$y <- ll4(pmax(d$dose, 1e-8), 2, 8 + u2[d$grp], -1.5, 3) + rnorm(nrow(d), 0, 0.3)
  f1 <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  # an optimum with the standard deviation near zero, as from a start in the flat region
  obj <- f1$obj
  p <- f1$opt$par
  p[["theta"]] <- log(1e-6)
  o <- suppressWarnings(nlminb(p, obj$fn, obj$gr, upper=replace(rep(Inf, length(p)), names(p) == "theta", log(1e-6))))
  o2 <- drcTMB:::retryvariances(obj, o, f1$reterms, rep(-Inf, length(p)), rep(Inf, length(p)), list())
  expect_lt(o2$objective, o$objective - 1)
  expect_equal(o2$objective, f1$opt$objective, tolerance=1e-6)
})
