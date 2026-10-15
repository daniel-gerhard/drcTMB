# Simulation-based checks: recovery of known parameters, moments of simulate(), consistency of
# the Laplace approximation, and (slow, only with DRCTMB_SLOW_TESTS=true) coverage of
# Wald confidence intervals.

# grouped loglogistic data with random effects on b1 and b4, correlation rho
simre <- function(ng=10, sd1=0.5, sd4=0, rho=0, sigma=0.3, b=c(1, 8, -1.5, 3, 0), seed=9){
  d <- simgroups(ng=ng, seed=seed)
  z1 <- rnorm(ng)
  z2 <- rho*z1 + sqrt(1 - rho^2)*rnorm(ng)
  u1 <- sd1*z1
  u4 <- sd4*z2
  d$y <- sapply(seq_len(nrow(d)), function(i){
    g <- d$grp[i]
    drcmean("loglogistic", d$dose[i], b + c(u1[g], 0, 0, u4[g], 0))
  }) + rnorm(nrow(d), 0, sigma)
  d
}

within <- function(est, truth, se, k=3) all(abs(est - truth) < k*se)

# log standard deviations of the random effects with their standard errors
logsd <- function(fit){
  ests <- fit$ssdrl[rownames(fit$ssdrl) == "theta", , drop=FALSE]
  # theta holds, for each term, the log standard deviations followed by the correlations
  np <- sapply(fit$reterms, function(tm) length(tm$pars))
  off <- cumsum(c(0, np + np*(np - 1)/2))
  id <- unlist(lapply(seq_along(np), function(i) off[i] + seq_len(np[i])))
  ests[id, , drop=FALSE]
}

test_that("fixed effects and variance components are recovered (gaussian, correlated b1 and b4)", {
  d <- simre(ng=40, sd1=0.5, sd4=0.4, rho=0.5, seed=10)
  f <- fit4(y ~ dose, rform=b1 + b4 ~ 1 | grp, data=d, model="loglogistic")
  se <- sqrt(diag(vcov(f)))[1:4]
  expect_true(within(bcoef(f), c(1, 8, -1.5, 3), se))
  ls <- logsd(f)
  expect_true(within(ls[, 1], log(c(0.5, 0.4)), ls[, 2]))
  expect_gt(drcTMB:::recov(f)[[1]]$corr[1, 2], 0)
  expect_equal(unname(exp(coef(f)["log_sigma"])), 0.3, tolerance=0.1)
})

test_that("fixed effects and variance components are recovered (binomial, random b4)", {
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=40, nrep=1, seed=11)
  u <- rnorm(40, 0, 0.6)
  d$n <- 30
  d$y <- rbinom(nrow(d), d$n, 0.05 + 0.9*plogis(1.2*(d$dose - 5 - u[d$grp])))
  # the logit link keeps the asymptotes inside (0, 1) during the optimisation
  f <- drmTMB(cbind(y, n - y) ~ dose, rform=b4 ~ 1 | grp, data=d, model="logistic", family="binomial",
              link="logit", start=c(b5=0), fix=c(b5=TRUE))
  se <- sqrt(diag(vcov(f)))[1:4]
  expect_true(within(bcoef(f), c(qlogis(0.05), qlogis(0.95), 1.2, 5), se))
  ls <- logsd(f)
  expect_true(within(ls[, 1], log(0.6), ls[, 2]))
})

test_that("simulate has the mean and variance of the fitted gaussian model with random effects", {
  d <- simre(seed=9)
  f <- fit4(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="loglogistic")
  nsim <- 2000
  sig2 <- exp(2*coef(f)[["log_sigma"]])
  # b1 enters the mean as b1 (1 - F(x)): variance (1 - F)^2 var(u) + sigma^2
  Fx <- drcmean("loglogistic", d$dose, c(0, 1, bcoef(f)[3:4], 0))
  vu <- drcTMB:::recov(f)[[1]]$sd^2

  # new random effects: mean of the population curve
  s <- simulate(f, nsim=nsim, seed=1)
  expect_equal(dim(s), c(nrow(d), nsim))
  expect_true(all(abs(rowMeans(s) - predict(f, random=FALSE)) < 4*sqrt(((1 - Fx)^2*vu + sig2)/nsim)))
  expect_equal(mean(apply(s, 1, var) / ((1 - Fx)^2*vu + sig2)), 1, tolerance=0.02)

  # fitted random effects: mean of the group curves, residual variance only
  s0 <- simulate(f, nsim=nsim, seed=1, random=FALSE)
  expect_true(all(abs(rowMeans(s0) - predict(f)) < 4*sqrt(sig2/nsim)))
  expect_equal(mean(apply(s0, 1, var)) / sig2, 1, tolerance=0.02)
})

test_that("simulate has the mean and variance of the poisson and binomial families", {
  set.seed(12)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=3), n=10)
  mu <- drcmean("loglogistic", d$dose, c(0.1, 0.9, 1.5, 2.5, 0))
  d$yb <- rbinom(nrow(d), d$n, mu)
  d$yp <- rpois(nrow(d), 40*mu)
  nsim <- 2000

  fp <- fit4(yp ~ dose, data=d, model="loglogistic", family="poisson")
  s <- simulate(fp, nsim=nsim, seed=1)
  m <- predict(fp)
  expect_true(all(abs(rowMeans(s) - m) < 4*sqrt(m/nsim)))
  expect_equal(mean(apply(s, 1, var) / m), 1, tolerance=0.02)

  fb <- fit4(cbind(yb, n - yb) ~ dose, data=d, model="loglogistic", family="binomial")
  s <- simulate(fb, nsim=nsim, seed=1)
  p <- predict(fb)
  expect_true(all(s >= 0 & s <= d$n))
  expect_true(all(abs(rowMeans(s) - d$n*p) < 4*sqrt(d$n*p*(1 - p)/nsim)))
  expect_equal(mean(apply(s, 1, var) / (d$n*p*(1 - p))), 1, tolerance=0.03)
})

test_that("the marginal score of the Laplace approximation has mean zero", {
  skip_on_cran()
  d <- simre(seed=9)
  f <- fit4(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="loglogistic")
  set.seed(1)
  cc <- summary(TMB::checkConsistency(f$obj, n=200))
  expect_gt(cc$marginal$p.value, 0.01)
  expect_true(all(abs(cc$marginal$bias) < 0.1))
})

test_that("Wald intervals have close to nominal coverage", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("DRCTMB_SLOW_TESTS"), "true"), "set DRCTMB_SLOW_TESTS=true to run")
  nsim <- 300
  res <- t(sapply(seq_len(nsim), function(i){
    d <- simre(ng=10, sd1=0.5, seed=1000 + i)
    f <- suppressWarnings(fit4(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="loglogistic"))
    se <- sqrt(diag(vcov(f)))[1:4]
    e <- ed(f, random=FALSE)
    c(abs(bcoef(f) - c(1, 8, -1.5, 3)) < qnorm(0.975)*se,
      ed=abs(e$estimates - 3) < qnorm(0.975)*sqrt(e$variance[1]))
  }))
  cover <- colMeans(res)
  message("coverage: ", paste(names(cover), round(cover, 3), collapse=", "))
  # 300 simulations: Monte Carlo standard error of the coverage about 0.013
  expect_true(all(cover > 0.90 & cover < 0.99), label=paste(round(cover, 3), collapse=", "))
})
