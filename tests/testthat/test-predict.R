# Standard errors, confidence and prediction intervals of predict

test_that("standard errors equal the delta method and agree with drc", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  nd <- data.frame(dose=c(0, 1, 3, 10))
  p <- predict(f, nd, se.fit=TRUE)
  expect_named(p, c("fit", "se.fit"))
  expect_equal(p$fit, predict(f, nd))

  skip_if_not_installed("numDeriv")
  G <- numDeriv::jacobian(function(b) drcmean("loglogistic", nd$dose, c(b, 0)), coef(f)[1:4])
  expect_equal(p$se.fit, sqrt(diag(G %*% vcov(f)[1:4, 1:4] %*% t(G))), tolerance=1e-6)

  # the unbiased residual variance, as drc
  skip_if_not_installed("drc")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4())
  se2 <- predict(f2, nd, se.fit=TRUE)[, "SE"]
  expect_equal(unname(se2), p$se.fit, tolerance=1e-4)
  pc <- predict(f, nd, interval="confidence")
  pc2 <- predict(f2, nd, interval="confidence")
  expect_equal(unname(as.matrix(pc[c("lwr", "upr")])), unname(pc2[, c("Lower", "Upper")]), tolerance=1e-4)
  pp <- predict(f, nd, interval="prediction")
  pp2 <- predict(f2, nd, interval="prediction")
  expect_equal(unname(as.matrix(pp[c("lwr", "upr")])), unname(pp2[, c("Lower", "Upper")]), tolerance=1e-4)
})

test_that("confidence and prediction intervals of gaussian models", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  # t quantiles and the unbiased residual variance
  z <- qt(0.95, nrow(d) - 4)
  ci <- predict(f, interval="confidence", level=0.9)
  expect_equal(ci$lwr, ci$fit - z*ci$se.fit)
  expect_equal(ci$upr, ci$fit + z*ci$se.fit)
  pi <- predict(f, interval="prediction", level=0.9)
  sig <- sqrt(sum(residuals(f, type="response")^2)/(nrow(d) - 4))
  expect_equal(pi$upr, pi$fit + z*sqrt(pi$se.fit^2 + sig^2))

  # maximum likelihood: normal quantiles and the ML residual standard deviation
  fml <- fit4(y ~ dose, data=d, model="loglogistic", control=list(smallsample=FALSE))
  z <- qnorm(0.95)
  ci <- predict(fml, interval="confidence", level=0.9)
  expect_equal(ci$upr, ci$fit + z*ci$se.fit)
  expect_equal(ci$se.fit, predict(f, se.fit=TRUE)$se.fit*sqrt((nrow(d) - 4)/nrow(d)))
  pi <- predict(fml, interval="prediction", level=0.9)
  sig <- exp(coef(fml)[["log_sigma"]])
  expect_equal(pi$upr, pi$fit + z*sqrt(pi$se.fit^2 + sig^2))
})

test_that("prediction intervals for groups in the data and for new groups", {
  d <- simgroups(ng=30, nrep=2, seed=3)
  u <- rnorm(30, 0, 0.8)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  nd <- data.frame(dose=c(0, 2, 16, 0, 2, 16), grp=factor(rep(c("1", "new"), each=3)))
  pc <- predict(f, nd, interval="prediction")
  pn <- predict(f, nd, random=FALSE, interval="prediction")
  expect_true(all(pc$se.fit > 0))
  # a level that is not in the data is a new group
  expect_equal(pc[4:6, ], pn[4:6, ], ignore_attr=TRUE)
  g1 <- d$grp == "1"
  expect_equal(pc$fit[1:3], predict(f)[g1][match(c(0, 2, 16), d$dose[g1])])
  # the variance of the random effects of new groups widens the intervals
  expect_true(all(pn$upr[1:3] - pn$lwr[1:3] > pc$upr[1:3] - pc$lwr[1:3]))

  # new-group intervals agree with simulated new groups, at the estimates
  set.seed(1)
  cf <- coef(f)
  sdb2 <- exp(cf[["log_sd:b2|grp"]])
  sig <- exp(cf[["log_sigma"]])
  sim <- replicate(20000, drcmean("loglogistic", nd$dose[1:3], c(cf[1], cf[2] + rnorm(1, 0, sdb2), cf[3], cf[4], 0)) + rnorm(3, 0, sig))
  q <- t(apply(sim, 1, quantile, c(0.025, 0.975)))
  z <- qnorm(0.975)
  sdp <- (pn$upr[1:3] - pn$fit[1:3])/z
  sdp0 <- sqrt(sdp^2 - pn$se.fit[1:3]^2)
  expect_equal(unname(q[, 2]), pn$fit[1:3] + z*sdp0, tolerance=0.02)
  expect_equal(unname(q[, 1]), pn$fit[1:3] - z*sdp0, tolerance=0.05)
})

test_that("confidence intervals stay in the range of the mean", {
  set.seed(2)
  dp <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16), each=3))
  dp$y <- rpois(nrow(dp), 1 + 20*exp(-dp$dose/3))
  fp <- drmTMB(y ~ dose, data=dp, model="exponentialdecay", family="poisson", link="log")
  cp <- predict(fp, interval="confidence")
  expect_true(all(cp$lwr > 0 & cp$lwr < cp$fit & cp$fit < cp$upr))

  db <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=2), n=20)
  db$s <- rbinom(nrow(db), 20, drcmean("loglogistic", db$dose, c(0.05, 0.95, -2, 2, 0)))
  fb <- fit4(cbind(s, n - s) ~ dose, data=db, model="loglogistic", family="binomial", link="logit")
  cb <- predict(fb, interval="confidence")
  expect_true(all(cb$lwr > 0 & cb$upr < 1 & cb$lwr < cb$fit & cb$fit < cb$upr))
})

test_that("bootstrap intervals agree with the delta method and the analytic prediction intervals", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic", nrep=10), model="loglogistic")
  bs <- bootstrap(f, nsim=300, seed=1, progress=FALSE)
  nd <- data.frame(dose=c(0, 1, 3, 10))
  cd <- predict(f, nd, interval="confidence")
  cb <- predict(f, nd, interval="confidence", boot=bs)
  expect_equal(cb$fit, cd$fit)
  expect_equal(cb$se.fit, cd$se.fit, tolerance=0.15)
  expect_equal(cb$upr - cb$lwr, cd$upr - cd$lwr, tolerance=0.15)
  pa <- predict(f, nd, interval="prediction")
  pb <- predict(f, nd, interval="prediction", boot=bs, seed=1)
  expect_equal(pb$upr - pb$lwr, pa$upr - pa$lwr, tolerance=0.1)
  expect_identical(pb, predict(f, nd, interval="prediction", boot=bs, seed=1))
  # a single prediction
  expect_equal(predict(f, nd[2, , drop=FALSE], interval="confidence", boot=bs), cb[2, ], ignore_attr=TRUE)
  expect_equal(nrow(predict(f, nd[2, , drop=FALSE], interval="prediction", boot=bs, seed=1)), 1)
  # too few converged fits
  b1 <- bs
  b1$estimates[-1, ] <- NA
  expect_error(predict(f, nd, interval="confidence", boot=b1), "only 1 of 300 bootstrap fits converged")
  b1$estimates[] <- NA
  expect_error(predict(f, nd, interval="prediction", boot=b1), "only 0 of 300")

  # with random effects only for population predictions
  d <- simgroups(ng=8, seed=7)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + rnorm(8, 0, 0.5)[d$grp] + rnorm(nrow(d), 0, 0.3)
  g <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  bg <- bootstrap(g, nsim=5, seed=1, progress=FALSE)
  expect_silent(predict(g, random=FALSE, interval="confidence", boot=bg))
})

test_that("bootstrap intervals with the random effects of the fitted groups", {
  d <- simgroups(ng=20, seed=8)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + rnorm(20, 0, 0.8)[d$grp] + rnorm(nrow(d), 0, 0.3)
  g <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  r <- g$obj$env$random

  # at the estimates, the drawn random effects follow their conditional distribution: mean at the
  # predicted random effects, covariance the inverse Hessian (without the uncertainty of the fixed
  # parameters, unlike sdreport)
  b0 <- list(estimates=matrix(coef(g), 4000, length(coef(g)), byrow=TRUE), nsim=4000)
  set.seed(1)
  u <- bootdraws(g, b0, random=TRUE)[r, ]
  sdr <- sdreport(g$obj, par.fixed=g$opt$par)
  expect_equal(rowMeans(u), unname(sdr$par.random), tolerance=0.02)
  g$obj$fn(g$opt$par)
  H <- as.matrix(g$obj$env$spHess(g$obj$env$last.par, random=TRUE))
  expect_equal(apply(u, 1, sd), sqrt(diag(solve(H))), tolerance=0.05)
  expect_true(all(sqrt(diag(solve(H))) < sqrt(sdr$diag.cov.random)))
  expect_true(all(bootdraws(g, b0, random=FALSE)[r, ] == 0))

  # intervals for fitted groups close to the delta method, and a new group gets the population curve
  bg <- bootstrap(g, nsim=200, seed=1, progress=FALSE)
  nd <- data.frame(dose=c(0, 2, 16, 2), grp=factor(c("1", "1", "1", "new")))
  cd <- predict(g, nd, interval="confidence")
  cb <- predict(g, nd, interval="confidence", boot=bg, seed=1)
  expect_identical(cb, predict(g, nd, interval="confidence", boot=bg, seed=1))
  expect_equal(cb$fit, cd$fit)
  expect_equal(cb$se.fit[1:3], cd$se.fit[1:3], tolerance=0.2)
  expect_equal(cb[4, ], predict(g, nd[4, ], random=FALSE, interval="confidence", boot=bg, seed=1), ignore_attr=TRUE)
  pb <- predict(g, nd, interval="prediction", boot=bg, seed=1)
  pa <- predict(g, nd, interval="prediction")
  expect_equal(pb$upr - pb$lwr, pa$upr - pa$lwr, tolerance=0.15)
})

test_that("the correlation matrix of a random effect term follows TMB", {
  d <- simgroups(ng=30, nrep=2, seed=4)
  u <- matrix(rnorm(60), 30) %*% chol(matrix(c(0.25, 0.3, 0.3, 0.64), 2))
  F <- drcmean("loglogistic", d$dose, c(0, 1, -1.5, 3, 0))
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp, 1]*(1 - F) + u[d$grp, 2]*F + rnorm(nrow(d), 0, 0.3)
  g <- fit4(y ~ dose, rform=b1 + b2 ~ 1 | grp, data=d, model="loglogistic")
  th <- g$obj$env$parList(g$opt$par)$theta
  expect_equal(corrmatrix(th[3], 2), recov(g)[[1]]$corr, ignore_attr=TRUE)
})

test_that("simulated prediction intervals of poisson and binomial models", {
  set.seed(2)
  dp <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16), each=3))
  dp$y <- rpois(nrow(dp), 1 + 20*exp(-dp$dose/3))
  fp <- drmTMB(y ~ dose, data=dp, model="exponentialdecay", family="poisson", link="log")
  pp <- predict(fp, data.frame(dose=c(0, 2, 8)), interval="prediction", seed=1)
  expect_identical(pp, predict(fp, data.frame(dose=c(0, 2, 8)), interval="prediction", seed=1))
  expect_equal(pp$lwr, round(pp$lwr))
  # wider than the poisson distribution at the estimated means
  expect_true(all(pp$lwr <= qpois(0.025, pp$fit) & pp$upr >= qpois(0.975, pp$fit)))
  expect_true(all(pp$lwr < pp$fit & pp$fit < pp$upr))

  d <- simgroups(ng=20, nrep=1, seed=6)
  d$n <- 20
  u <- rnorm(20, 0, 1)
  d$s <- rbinom(nrow(d), 20, plogis(qlogis(drcmean("loglogistic", d$dose, c(0.05, 0.9, -2, 2, 0))) + u[d$grp]))
  fb <- fit4(cbind(s, n - s) ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", family="binomial", link="logit")
  nd <- data.frame(dose=c(0, 2), grp=factor("1"), s=0, n=20)
  pc <- predict(fb, nd, interval="prediction", seed=1)
  pn <- predict(fb, nd, random=FALSE, interval="prediction", seed=1)
  expect_true(all(pc$lwr >= 0 & pc$upr <= 1 & pn$lwr >= 0 & pn$upr <= 1))
  expect_true(all(pn$upr - pn$lwr > pc$upr - pc$lwr))
  expect_error(predict(fb, data.frame(dose=1, grp="1"), interval="prediction"), "number of trials")
})
