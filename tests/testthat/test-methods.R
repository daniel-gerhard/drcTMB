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

test_that("coef names the variance parameters", {
  d <- simgroups(ng=30, seed=18)
  u <- matrix(rnorm(120, 0, 0.4), 30)
  d$y <- 1 + u[d$grp, 1] + (7 + u[d$grp, 2] - u[d$grp, 1]) /
    (1 + exp((1.5 + u[d$grp, 4])*(log(d$dose) - log(3*exp(u[d$grp, 3]))))) + rnorm(nrow(d), 0, 0.2)
  f <- fit4(y ~ dose, rform=list(b1 + b2 + b4 ~ 1 | grp, b3 ~ 1 | grp), data=d, model="loglogistic")
  cf <- coef(f)
  expect_equal(names(cf)[-(1:4)],
               c("log_sd:b1|grp", "log_sd:b2|grp", "log_sd:b4|grp",
                 "corr_par:b2,b1|grp", "corr_par:b4,b1|grp", "corr_par:b4,b2|grp",
                 "log_sd:b3|grp", "log_sigma"))
  expect_equal(rownames(vcov(f)), names(cf))

  # the correlation matrix from the named parameters, as in TMB: lower triangle of L by row
  L <- diag(3)
  L[2, 1] <- cf[["corr_par:b2,b1|grp"]]
  L[3, 1] <- cf[["corr_par:b4,b1|grp"]]
  L[3, 2] <- cf[["corr_par:b4,b2|grp"]]
  expect_equal(cov2cor(L %*% t(L)), drcTMB:::recov(f)[[1]]$corr, tolerance=1e-6)
  expect_equal(exp(unname(cf[c("log_sd:b1|grp", "log_sd:b2|grp", "log_sd:b4|grp")])), unname(drcTMB:::recov(f)[[1]]$sd))

  # beta: log precision
  set.seed(19)
  db <- data.frame(dose=rep(c(0, 1, 2, 4, 8), each=6))
  mu <- drcmean("loglogistic", db$dose, c(0.1, 0.9, 1.5, 2.5, 0))
  db$y <- rbeta(nrow(db), 30*mu, 30*(1 - mu))
  expect_equal(tail(names(coef(fit4(y ~ dose, data=db, model="loglogistic", family="beta"))), 1), "log_phi")
})

test_that("predict uses the factor levels of the fit for new data", {
  d <- simgauss("loglogistic")
  d$trt <- factor(rep(c("a", "b"), length.out=nrow(d)))
  d$y <- d$y + 3*(d$trt == "b")
  f <- fit4(y ~ dose, fform=b2 ~ trt, data=d, model="loglogistic")
  p <- predict(f)
  i <- match(c("a", "b"), d$trt)
  # levels in a different order, a single level, and a character variable
  nd <- data.frame(dose=d$dose[i], trt=factor(c("a", "b"), levels=c("b", "a")))
  expect_equal(predict(f, newdata=nd), p[i])
  expect_equal(predict(f, newdata=data.frame(dose=d$dose[i[2]], trt=factor("b"))), p[i[2]])
  expect_equal(predict(f, newdata=data.frame(dose=d$dose[i], trt=c("a", "b"))), p[i])
})

test_that("observations with missing values are left out", {
  d <- simgauss("loglogistic")
  d$w <- rep(1:2, length.out=nrow(d))
  f <- fit4(y ~ dose, data=d, weights=d$w, model="loglogistic")
  dn <- rbind(d, data.frame(dose=c(NA, 1, 2), y=c(5, NA, 5), w=c(1, 1, NA)))
  fn <- fit4(y ~ dose, data=dn, weights=dn$w, model="loglogistic")
  expect_equal(coef(fn), coef(f), tolerance=1e-6)
  expect_equal(nobs(fn), nrow(d))
  expect_length(predict(fn), nrow(d))
  expect_length(residuals(fn), nrow(d))
  expect_error(fit4(y ~ dose, data=d, weights=1:3, model="loglogistic"), "weights")
})

test_that("ed includes the residual variance only for the gaussian family", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  e0 <- ed(f)
  e1 <- ed(f, residual=TRUE)
  expect_equal(coef(e1), coef(e0))
  expect_gt(vcov(e1)[1, 1], vcov(e0)[1, 1])
  set.seed(3)
  dp <- data.frame(dose=rep(c(0, 1, 2, 4, 8), each=3))
  dp$y <- rpois(nrow(dp), 1 + 20*exp(-dp$dose/3))
  fp <- drmTMB(y ~ dose, data=dp, model="exponentialdecay", family="poisson", link="log")
  expect_error(ed(fp, residual=TRUE), "gaussian")
})

test_that("unknown families and links, and linfct of the wrong size, give clear errors", {
  d <- simgauss("loglogistic")
  expect_error(drmTMB(y ~ dose, data=d, family="gausian"), "Unknown family")
  expect_error(drmTMB(y ~ dose, data=d, link="logg"), "Unknown link")
  expect_error(drmTMBstart(y ~ dose, data=d, family="gausian"), "Unknown family")
  d$trt <- factor(rep(c("a", "b"), length.out=nrow(d)))
  f <- fit4(y ~ dose, fform=b2 ~ trt, data=d, model="loglogistic")
  expect_error(ed(f), "6 model coefficients")
  expect_error(ed(f, linfct=list(diag(6)[-3, ], diag(5))), "5 rows")
  expect_length(coef(ed(f, linfct=diag(6)[-3, ])), 1)
})

test_that("print returns the effective doses invisibly", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  e <- ed(f)
  out <- capture.output(res <- withVisible(print(e)))
  expect_false(res$visible)
  expect_identical(res$value, e)
  expect_true(any(grepl("^Effective dose estimates", out)))
})

test_that("methods for the negative binomial families", {
  set.seed(3)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16, 32), each=5))
  d$y <- rnbinom(nrow(d), size=5, mu=drcmean("loglogistic", d$dose, c(2, 60, -1.5, 3, 0)))
  f2 <- fit4(y ~ dose, data=d, model="loglogistic", family="nbinom2")
  f1 <- fit4(y ~ dose, data=d, model="loglogistic", family="nbinom1")
  expect_equal(names(coef(f2))[5], "log_theta")
  expect_equal(names(coef(f1))[5], "log_phi")
  expect_true(any(grepl("^Negative binomial size \\(theta\\)", capture.output(print(f2)))))
  expect_true(any(grepl("^Negative binomial dispersion \\(phi\\)", capture.output(print(summary(f1))))))
  expect_null(summary(f2)$dispersion)

  # Pearson residuals with the negative binomial variance; quantile residuals of a correct model
  mu <- fitted(f2)
  th <- exp(coef(f2)[["log_theta"]])
  expect_equal(residuals(f2), (d$y - mu)/sqrt(mu*(1 + mu/th)))
  ph <- exp(coef(f1)[["log_phi"]])
  expect_equal(residuals(f1), (d$y - fitted(f1))/sqrt(fitted(f1)*(1 + ph)))
  set.seed(1)
  expect_true(abs(mean(residuals(f2, type="quantile"))) < 0.5)

  # intervals: log scale for the mean, integer quantiles for new counts, size on the natural scale
  pc <- predict(f2, data.frame(dose=c(0, 3)), interval="confidence")
  expect_true(all(pc$lwr > 0 & pc$lwr < pc$fit & pc$fit < pc$upr))
  pp <- predict(f2, data.frame(dose=c(0, 3)), interval="prediction", seed=1)
  expect_true(all(pp$lwr == round(pp$lwr) & pp$upr == round(pp$upr)))
  expect_true(all(pp$lwr < pc$lwr & pp$upr > pc$upr))
  expect_equal(rownames(confint(f2, parm="log_theta", method="wald")), "theta")
  expect_equal(dim(simulate(f2, 2, seed=1)), c(nrow(d), 2))
  # the overdispersed counts are much better described than by the poisson model
  fp <- fit4(y ~ dose, data=d, model="loglogistic", family="poisson")
  expect_lt(AIC(f2), AIC(fp) - 50)
})

test_that("methods for the beta-binomial family", {
  set.seed(4)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16, 32), each=4), n=25)
  p <- drcmean("loglogistic", d$dose, c(0.05, 0.9, -1.5, 3, 0))
  d$y <- rbinom(nrow(d), d$n, rbeta(nrow(d), p*6, (1 - p)*6))
  f <- fit4(cbind(y, n - y) ~ dose, data=d, model="loglogistic", family="betabinomial", link="logit", b4link="log")
  expect_equal(names(coef(f))[5], "log_phi")
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("^Beta-binomial precision \\(phi\\)", out)))
  expect_true(any(grepl("number of trials: 800", out)))
  expect_null(summary(f)$dispersion)

  # Pearson residuals on the scale of proportions with the beta-binomial variance
  mu <- fitted(f)
  ph <- exp(coef(f)[["log_phi"]])
  expect_equal(residuals(f), (d$y/d$n - mu)/sqrt(mu*(1 - mu)/d$n*(1 + (d$n - 1)/(ph + 1))))
  # distribution function: complete at n, binomial for a large precision
  expect_equal(drcTMB:::pbetabinom(20, 20, 0.3, 5), 1)
  expect_equal(drcTMB:::pbetabinom(c(-1, 3, 8), 20, 0.3, 1e8), c(0, pbinom(c(3, 8), 20, 0.3)), tolerance=1e-6)
  set.seed(1)
  expect_true(abs(mean(residuals(f, type="quantile"))) < 0.5)

  # prediction intervals of proportions need the number of trials
  pp <- predict(f, data.frame(dose=c(1, 3), y=0, n=25), interval="prediction", seed=1)
  expect_true(all(pp$lwr*25 == round(pp$lwr*25)))
  expect_error(predict(f, data.frame(dose=1), interval="prediction"), "number of trials")
  expect_equal(dim(simulate(f, 2, seed=1)), c(nrow(d), 2))
  # overdispersed proportions are much better described than by the binomial model
  fb <- fit4(cbind(y, n - y) ~ dose, data=d, model="loglogistic", family="binomial", link="logit", b4link="log")
  expect_lt(AIC(f), AIC(fb) - 20)
  # binary responses
  d$bin <- as.integer(d$y > 12)
  expect_error(fit4(bin ~ dose, data=d, model="loglogistic", family="betabinomial"), "not identifiable for binary")
})
