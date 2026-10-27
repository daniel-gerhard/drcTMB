test_that("plot gives diagnostic ggplots of the residuals", {
  d <- simgauss("loglogistic")
  f4 <- fit4(y ~ dose, data=d, model="loglogistic")
  pdf(NULL)
  on.exit(dev.off())
  p <- plot(f4, seed=1)
  expect_named(p, c("fitted", "dose", "qq"))
  expect_true(all(vapply(p, ggplot2::is_ggplot, logical(1))))
  # the residuals of the fit, quantile residuals by default
  expect_equal(p$fitted$data$res, residuals(f4, type="quantile"))
  expect_equal(p$fitted$data$xv, fitted(f4))
  # log dose axis, dose zero drawn left of the smallest positive dose
  xd <- p$dose$data$xv
  expect_true(all(xd > 0))
  expect_lt(max(xd[d$dose == 0]), min(d$dose[d$dose > 0]))
  expect_equal(xd[d$dose > 0], d$dose[d$dose > 0])
  expect_equal(ggplot2::layer_scales(p$dose)$x$get_labels()[1], "0")
  expect_null(plot(f4, which=2, logdose=FALSE, draw=FALSE)$dose$scales$get_scales("x"))
  # pearson residuals, without the fitted random effects
  pp <- plot(f4, which=c(1, 3), type="pearson", draw=FALSE)
  expect_named(pp, c("fitted", "qq"))
  expect_equal(sort(pp$qq$data$res), sort(residuals(f4)))
  expect_error(plot(f4, which=4), "no random effects")
  expect_named(plot(f4, which=3:4, draw=FALSE), "qq")
  expect_error(plot(f4, which=5), "subset of 1:4")
})

test_that("plot shows the predicted random effects of each term", {
  d <- simgroups(ng=8, nrep=1, seed=5)
  u <- rnorm(8, 0, 0.5)
  d$y <- rnorm(nrow(d), drcmean("logistic", d$dose, c(1, 8, 1, 4, 0)) + u[d$grp], 0.3)
  f <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="logistic", start=c(b5=0), fix=c(b5=TRUE))
  p <- plot(f, which=4, draw=FALSE)
  expect_equal(unique(p$ranef$data$term), "b2 | grp")
  expect_equal(p$ranef$data$res, unname(f$sdrl$par.random))
  # standard normal quantiles in the order of the random effects
  expect_equal(p$ranef$data$theo[order(p$ranef$data$res)], qnorm(ppoints(8)))
  # the linear dose scale of the logistic model
  expect_equal(plot(f, which=2, draw=FALSE)$dose$data$xv, d$dose)
})
