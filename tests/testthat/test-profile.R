# Likelihood profiles. Profile likelihood intervals are invariant to monotone reparameterisations,
# so they can be compared exactly with other software for equivalent models.

test_that("profile intervals of a gaussian model are likelihood ratio intervals", {
  d <- simgauss("loglogistic")
  d$dose[d$dose == 0] <- 1e-8
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  p <- profile(f)
  expect_s3_class(p, "profile.drmTMB")
  expect_named(p, names(bcoef(f)))
  ci <- confint(p)
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_true(all(ci[, 1] < bcoef(f) & bcoef(f) < ci[, 2]))

  # nls with b4 fixed at a limit: the deviance n log(RSS / RSSmin) equals the chi-square quantile
  rss0 <- sum((d$y - predict(f))^2)
  for (b4 in ci[4, ]){
    n1 <- nls(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))), data=d,
              start=setNames(bcoef(f)[1:3], paste("b", 1:3, sep="")))
    expect_equal(nrow(d)*log(sum(residuals(n1)^2)/rss0), qchisq(0.95, 1), tolerance=1e-3)
  }

  # close to the Wald intervals for a nearly linear parameter
  expect_equal(unname(ci[2, ]), unname(confint(f)[2, ]), tolerance=0.02)
})

test_that("profile intervals with random effects agree with lme4", {
  skip_if_not_installed("lme4")
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=2)
  u <- rnorm(12, 0, 0.5)
  d$y <- rpois(nrow(d), exp(log(50) + u[d$grp] - d$dose/4))
  f <- drmTMB(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="exponentialdecay", family="poisson",
              link="log", start=c(b2=-30), fix=c(b2=TRUE))
  ci <- confint(profile(f, parm=names(coef(f))))
  expect_equal(rownames(ci), c("b1:(Intercept)", "b4:(Intercept)", "sd:b1|grp"))

  # log(mu) = b1 - dose / b4: intercept, -1 / slope, and the same standard deviation
  g <- lme4::glmer(y ~ dose + (1 | grp), data=d, family=poisson)
  cg <- suppressMessages(suppressWarnings(confint(g, method="profile")))
  expect_equal(unname(ci[1, ]), unname(cg["(Intercept)", ]), tolerance=1e-4)
  expect_equal(unname(ci[2, ]), unname(sort(-1/cg["dose", ])), tolerance=1e-4)
  expect_equal(unname(ci[3, ]), unname(cg[".sig01", ]), tolerance=1e-4)
})

test_that("standard deviations at the boundary have a lower limit of zero", {
  d <- simgroups(dose=c(0, 1, 2, 4, 8, 16), ng=6, nrep=2, seed=14)
  d$y <- rpois(nrow(d), drcmean("loglogistic", d$dose, c(5, 45, -1.2, 3, 0)))
  f <- suppressWarnings(fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", family="poisson"))
  ci <- suppressWarnings(confint(profile(f, parm="log_sd:b2|grp")))
  expect_equal(rownames(ci), "sd:b2|grp")
  expect_equal(ci[1, 1], 0)
  expect_true(is.finite(ci[1, 2]) && ci[1, 2] > 0)
})

test_that("intervals beyond the profiled level are missing with a warning", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  p <- profile(f, parm=4, level=0.8)
  expect_warning(ci <- confint(p, level=0.999), "do not reach")
  expect_true(any(is.na(ci)))
  expect_silent(confint(p, level=0.5))
})

test_that("the profiled coefficient stays within its bounds", {
  d <- simgauss("loglogistic")
  f0 <- fit4(y ~ dose, data=d, model="loglogistic")
  up <- c(Inf, Inf, Inf, bcoef(f0)[4] + 0.1, Inf)
  f <- fit4(y ~ dose, data=d, model="loglogistic", upper=up)
  p <- profile(f, parm=4)
  expect_lte(max(p[[1]]$par), up[4] + 1e-10)
  expect_warning(ci <- confint(p), "do not reach")
  expect_true(is.na(ci[1, 2]))
})

test_that("plot draws the profiles with ggplot2", {
  skip_if_not_installed("ggplot2")
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  p <- profile(f, parm=3:4)
  g <- plot(p)
  expect_s3_class(g, "ggplot")
  expect_equal(levels(g$data$parameter), names(p))
  expect_equal(nrow(g$data), sum(sapply(p, nrow)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(print(g))
})
