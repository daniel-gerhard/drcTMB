# Some dose-response models are generalised linear (mixed) models in another parameterisation.
# Both drcTMB and lme4/glmmTMB use the Laplace approximation, so the marginal log-likelihoods
# agree exactly at the optimum, and the estimates agree after transformation.
#
#   binomial, logistic, b1 = 0, b2 = 1:    logit(p) = b3 x - b3 b4, random b4: sd = b3 sd(u)
#   poisson, exponentialdecay, log link,
#     b2 = -30 (exp(b2) ~ 0):              log(mu) = b1 - x / b4, random b1
#   beta, logistic, b1 = 0, b2 = 1:        logit(mu) as for binomial, precision phi
#   nbinom2, nbinom1: as for poisson, with size theta or dispersion phi as in glmmTMB

skip_if_not_installed("lme4")

asym01 <- list(start=c(b1=0, b2=1, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE))

simbinom <- function(seed=2){
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=seed)
  u <- rnorm(12, 0, 0.8)
  d$n <- 20
  d$y <- rbinom(nrow(d), d$n, plogis(0.9*(d$dose - (5 + u[d$grp]))))
  d
}

test_that("binomial logistic without random effects agrees with glm", {
  d <- simbinom()
  f1 <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="logistic", family="binomial",
               start=asym01$start, fix=asym01$fix)
  g1 <- glm(cbind(y, n - y) ~ dose, data=d, family=binomial)
  b <- coef(g1)
  expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-7)
  expect_equal(unname(bcoef(f1)), unname(c(b[2], -b[1]/b[2])), tolerance=1e-5)
})

test_that("binomial logistic with random b4 agrees with glmer", {
  d <- simbinom()
  f1 <- drmTMB(cbind(y, n - y) ~ dose, rform=b4 ~ 1 | grp, data=d, model="logistic", family="binomial",
               start=asym01$start, fix=asym01$fix)
  g1 <- lme4::glmer(cbind(y, n - y) ~ dose + (1 | grp), data=d, family=binomial)
  b <- lme4::fixef(g1)
  sdg <- sqrt(lme4::VarCorr(g1)$grp[1])
  expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-6)
  expect_equal(unname(bcoef(f1)), unname(c(b[2], -b[1]/b[2])), tolerance=1e-4)
  expect_equal(drcTMB:::recov(f1)[[1]]$sd, unname(sdg / b[2]), tolerance=1e-3)
})

test_that("poisson exponential decay with random b1 agrees with glmer", {
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=2)
  u <- rnorm(12, 0, 0.5)
  d$y <- rpois(nrow(d), exp(log(50) + u[d$grp] - d$dose/4))
  f1 <- drmTMB(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="exponentialdecay", family="poisson",
               link="log", start=c(b2=-30), fix=c(b2=TRUE))
  g1 <- lme4::glmer(y ~ dose + (1 | grp), data=d, family=poisson)
  b <- lme4::fixef(g1)
  expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-6)
  expect_equal(unname(bcoef(f1)), unname(c(b[1], -1/b[2])), tolerance=1e-4)
  expect_equal(drcTMB:::recov(f1)[[1]]$sd, sqrt(lme4::VarCorr(g1)$grp[1]), tolerance=1e-3)
})

test_that("beta logistic agrees with glmmTMB, with and without random effects", {
  skip_if_not_installed("glmmTMB")
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=2)
  u <- rnorm(12, 0, 0.5)
  mu <- plogis(0.5*(d$dose - 5) + u[d$grp])
  d$y <- rbeta(nrow(d), mu*30, (1 - mu)*30)

  f1 <- drmTMB(y ~ dose, data=d, model="logistic", family="beta", start=asym01$start, fix=asym01$fix)
  g1 <- glmmTMB::glmmTMB(y ~ dose, data=d, family=glmmTMB::beta_family())
  b <- glmmTMB::fixef(g1)$cond
  expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-6)
  expect_equal(unname(bcoef(f1)), unname(c(b[2], -b[1]/b[2])), tolerance=1e-4)
  expect_equal(unname(exp(coef(f1)["log_phi"])), sigma(g1), tolerance=1e-4)

  f2 <- drmTMB(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="logistic", family="beta",
               start=asym01$start, fix=asym01$fix)
  g2 <- glmmTMB::glmmTMB(y ~ dose + (1 | grp), data=d, family=glmmTMB::beta_family())
  b <- glmmTMB::fixef(g2)$cond
  expect_equal(f2$opt$objective, -as.numeric(logLik(g2)), tolerance=1e-6)
  expect_equal(unname(bcoef(f2)), unname(c(b[2], -b[1]/b[2])), tolerance=1e-4)
})

simnb <- function(seed=2){
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=seed)
  u <- rnorm(12, 0, 0.5)
  d$y <- rnbinom(nrow(d), size=4, mu=exp(log(50) + u[d$grp] - d$dose/4))
  d
}

test_that("nbinom2 exponential decay agrees with glm.nb", {
  skip_if_not_installed("MASS")
  d <- simnb()
  f1 <- drmTMB(y ~ dose, data=d, model="exponentialdecay", family="nbinom2", link="log", start=c(b2=-30), fix=c(b2=TRUE))
  g1 <- MASS::glm.nb(y ~ dose, data=d)
  b <- coef(g1)
  expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-7)
  expect_equal(unname(bcoef(f1)), unname(c(b[1], -1/b[2])), tolerance=1e-5)
  expect_equal(unname(exp(coef(f1)["log_theta"])), g1$theta, tolerance=1e-4)
})

test_that("negative binomial families with random b1 agree with glmmTMB", {
  skip_if_not_installed("glmmTMB")
  d <- simnb()
  for (fam in c("nbinom2", "nbinom1")){
    f1 <- drmTMB(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="exponentialdecay", family=fam, link="log",
                 start=c(b2=-30), fix=c(b2=TRUE))
    g1 <- glmmTMB::glmmTMB(y ~ dose + (1 | grp), data=d, family=get(fam, asNamespace("glmmTMB"))())
    b <- glmmTMB::fixef(g1)$cond
    expect_equal(f1$opt$objective, -as.numeric(logLik(g1)), tolerance=1e-6, info=fam)
    expect_equal(unname(bcoef(f1)), unname(c(b[1], -1/b[2])), tolerance=1e-4, info=fam)
    expect_equal(unname(exp(coef(f1)[length(coef(f1))])), sigma(g1), tolerance=1e-4, info=fam)
    expect_equal(drcTMB:::recov(f1)[[1]]$sd, sqrt(glmmTMB::VarCorr(g1)$cond$grp[1]), tolerance=1e-3, info=fam)
  }
})

test_that("Laplace approximation is close to the exact marginal likelihood for binomial data", {
  # the marginal likelihood of a single random effect per group, integrated numerically
  # (the log-likelihood of glmer with nAGQ > 1 has different constants, so it is not used here)
  d <- simbinom()
  f1 <- drmTMB(cbind(y, n - y) ~ dose, rform=b4 ~ 1 | grp, data=d, model="logistic", family="binomial",
               start=asym01$start, fix=asym01$fix)
  p <- f1$opt$par
  sdu <- exp(p[["theta"]])
  lik <- function(g){
    integrand <- function(u) sapply(u, function(uu) prod(dbinom(g$y, g$n, plogis(p[1]*(g$dose - p[2] - uu))))) * dnorm(u, 0, sdu)
    integrate(integrand, -Inf, Inf, rel.tol=1e-10)$value
  }
  exact <- -sum(log(sapply(split(d, d$grp), lik)))
  expect_equal(as.numeric(f1$obj$fn(p)), exact, tolerance=1e-4)
  expect_lt(abs(f1$opt$objective - exact), 0.01)
})
