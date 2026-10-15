test_that("print shows the call, model, coefficients and variance parameters", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b1 + b2 ~ 1 | grp, data=d, model="loglogistic")
  out <- capture.output(res <- print(f))
  expect_identical(res, f)
  expect_true(any(grepl("^Call:", out)))
  expect_true(any(grepl("^Model:  loglogistic", out)))
  expect_true(any(grepl("^Family: gaussian, identity link", out)))
  expect_true(any(grepl("^Fixed: b5:\\(Intercept\\) = 0", out)))
  expect_true(any(grepl("^ grp +b1", out)))
  expect_true(any(grepl("^Residual standard deviation", out)))
  expect_false(any(grepl("Std. Error", out)))

  # unused parameters are not listed as fixed; no residual line for poisson
  dp <- data.frame(dose=rep(c(0, 1, 2, 4, 8), each=3))
  dp$y <- rpois(nrow(dp), 20*exp(-dp$dose/3))
  out <- capture.output(print(drmTMB(y ~ dose, data=dp, model="exponentialdecay", family="poisson")))
  expect_false(any(grepl("^Fixed:", out)))
  expect_false(any(grepl("^Residual|^Precision", out)))
})

test_that("logLik, AIC, BIC and nobs agree with drc", {
  skip_if_not_installed("drc")
  d <- simgauss("loglogistic")
  f1 <- drmTMB(y ~ dose, data=d, model="loglogistic")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::LL.5())
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance=1e-6)
  expect_equal(attr(logLik(f1), "df"), attr(logLik(f2), "df"))
  expect_equal(AIC(f1), AIC(f2), tolerance=1e-6)
  expect_equal(BIC(f1), BIC(f2), tolerance=1e-6)
  expect_equal(nobs(f1), nrow(d))
})

test_that("logLik counts the variance parameters of mixed models", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b1 + b2 ~ 1 | grp, data=d, model="loglogistic")
  # 4 coefficients, 2 standard deviations, 1 correlation, residual sd
  expect_equal(attr(logLik(f), "df"), 8)
  expect_equal(as.numeric(logLik(f)), -f$opt$objective)
})

test_that("summary contains the coefficient table, sample sizes and variance components", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=list(b2 ~ 1 | grp, b4 ~ 1 | grp), data=d, model="loglogistic")
  s <- summary(f)
  expect_s3_class(s, "summary.drmTMB")
  expect_equal(colnames(coef(s)), c("Estimate", "Std. Error"))
  expect_equal(coef(s)[, 1], bcoef(f))
  expect_equal(unname(coef(s)[, 2]), unname(sqrt(diag(vcov(f)))[1:4]))
  expect_equal(s$nobs, nrow(d))
  expect_equal(s$ngroups, c(grp=8))
  expect_equal(s$direction, "decreasing")
  expect_equal(s$residpar, unname(exp(coef(f)["log_sigma"])))
  expect_length(s$varcor, 2)
  expect_true(s$pdHess)

  out <- capture.output(print(s))
  expect_true(any(grepl("^Log-likelihood", out)))
  expect_true(any(grepl("groups: grp 8", out)))
  expect_true(any(grepl("^Optimiser", out)))
  expect_false(any(grepl("z value|Pr\\(", out)))
})

test_that("summary reports a singular fit and the Pearson dispersion", {
  # no variation between groups: the standard deviation goes to zero
  set.seed(14)
  d <- simgroups(dose=c(0, 1, 2, 4, 8, 16), ng=6, nrep=2, seed=14)
  d$y <- rpois(nrow(d), drcmean("loglogistic", d$dose, c(5, 45, -1.2, 3, 0)))
  f <- suppressWarnings(fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", family="poisson"))
  s <- summary(f)
  expect_equal(s$singular, "grp b2")
  expect_equal(unname(s$dispersion["df"]), nrow(d) - attr(logLik(f), "df"))
  expect_true(abs(s$dispersion["ratio"] - 1) < 0.5)
  expect_true(any(grepl("singular fit", capture.output(print(s)))))
})

test_that("summary warns when the Hessian is not positive definite", {
  d <- simgauss("loglogistic")
  s <- summary(fit4(y ~ dose, data=d, model="loglogistic"))
  s$pdHess <- FALSE
  expect_true(any(grepl("not positive definite", capture.output(print(s)))))
})

test_that("fitted values and response residuals", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  expect_equal(fitted(f), predict(f))
  expect_equal(fitted(f, random=FALSE), predict(f, random=FALSE))
  expect_equal(residuals(f, type="response"), d$y - predict(f))
  expect_equal(residuals(f, type="response", random=FALSE), d$y - predict(f, random=FALSE))
  sig <- exp(coef(f)[["log_sigma"]])
  expect_equal(residuals(f), (d$y - predict(f))/sig)
  expect_equal(residuals(f, type="quantile"), residuals(f))
})

test_that("pearson residuals agree with glm", {
  set.seed(15)
  d <- data.frame(dose=rep(c(0, 1, 2, 4, 6, 8, 10), each=3), n=20)
  d$y <- rbinom(nrow(d), d$n, plogis(0.9*(d$dose - 5)))
  d$cnt <- rpois(nrow(d), exp(log(50) - d$dose/4))
  w <- rep(1:3, length.out=nrow(d))

  fb <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="logistic", family="binomial",
               start=c(b1=0, b2=1, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE))
  gb <- glm(cbind(y, n - y) ~ dose, data=d, family=binomial)
  expect_equal(residuals(fb), unname(residuals(gb, type="pearson")), tolerance=1e-5)
  expect_equal(residuals(fb, type="response"), unname(residuals(gb, type="response")), tolerance=1e-5)

  fp <- drmTMB(cnt ~ dose, data=d, model="exponentialdecay", family="poisson", link="log",
               start=c(b2=-30), fix=c(b2=TRUE), weights=w)
  gp <- glm(cnt ~ dose, data=d, family=poisson, weights=w)
  expect_equal(residuals(fp), unname(residuals(gp, type="pearson")), tolerance=1e-5)

  # sum of squares: the Pearson statistic of summary
  expect_equal(sum(residuals(fp)^2), unname(summary(fp)$dispersion["chisq"]))
})

test_that("quantile residuals are standard normal under the model", {
  set.seed(16)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=60), n=10)
  mu <- drcmean("loglogistic", d$dose, c(0.1, 0.9, 1.5, 2.5, 0))
  d$yb <- rbinom(nrow(d), d$n, mu)
  d$yp <- rpois(nrow(d), 10*mu)
  d$yr <- rbeta(nrow(d), 20*mu, 20*(1 - mu))
  fits <- list(
    binomial=fit4(cbind(yb, n - yb) ~ dose, data=d, model="loglogistic", family="binomial"),
    poisson=fit4(yp ~ dose, data=d, model="loglogistic", family="poisson"),
    beta=fit4(yr ~ dose, data=d, model="loglogistic", family="beta")
  )
  for (nm in names(fits)){
    set.seed(17)
    r <- residuals(fits[[nm]], type="quantile")
    expect_true(all(is.finite(r)), label=nm)
    expect_lt(abs(mean(r)), 0.1)
    expect_lt(abs(sd(r) - 1), 0.1)
    expect_gt(suppressWarnings(ks.test(r, "pnorm"))$p.value, 0.01)
  }
  # discrete families: randomised, reproducible with a seed
  set.seed(1); r1 <- residuals(fits$poisson, type="quantile")
  set.seed(1); r2 <- residuals(fits$poisson, type="quantile")
  expect_identical(r1, r2)
})
