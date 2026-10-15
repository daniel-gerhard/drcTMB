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
