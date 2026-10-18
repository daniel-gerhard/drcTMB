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

  # drc uses the unbiased residual variance
  skip_if_not_installed("drc")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4())
  se2 <- predict(f2, nd, se.fit=TRUE)[, "SE"]
  expect_equal(unname(se2), p$se.fit*sqrt(nrow(d)/(nrow(d) - 4)), tolerance=1e-4)
})

test_that("confidence and prediction intervals of gaussian models", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  z <- qnorm(0.95)
  ci <- predict(f, interval="confidence", level=0.9)
  expect_equal(ci$lwr, ci$fit - z*ci$se.fit)
  expect_equal(ci$upr, ci$fit + z*ci$se.fit)
  pi <- predict(f, interval="prediction", level=0.9)
  sig <- exp(coef(f)[["log_sigma"]])
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

test_that("confidence intervals stay in the range of the mean, and prediction intervals need a gaussian model", {
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
  expect_error(predict(fb, interval="prediction"), "gaussian")
})
