test_that("anova gives the likelihood ratio test of drc", {
  skip_if_not_installed("drc")
  d <- simgauss("loglogistic")
  f4 <- fit4(y ~ dose, data=d, model="loglogistic")
  f5 <- drmTMB(y ~ dose, data=d, model="loglogistic")
  a <- anova(f5, f4)
  expect_s3_class(a, "anova")
  # ordered by the number of parameters
  expect_equal(rownames(a), c("f4", "f5"))
  expect_equal(a$npar, c(5, 6))
  d4 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4())
  d5 <- drc::drm(y ~ dose, data=d, fct=drc::LL.5())
  ad <- drc:::anova.drc(d4, d5, test="Chisq", details=FALSE)
  expect_equal(a$Chisq[2], ad[2, "LR value"], tolerance=1e-5)
  expect_equal(a[["Pr(>Chisq)"]][2], ad[2, "p value"], tolerance=1e-5)
  expect_equal(a$AIC, c(AIC(f4), AIC(f5)))
  expect_equal(a$deviance, -2*a$logLik)
  expect_true(is.na(a$Chisq[1]) && is.na(a[["Pr(>Chisq)"]][1]))
})

test_that("anova of models with and without random effects agrees with lme4", {
  skip_if_not_installed("lme4")
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=2)
  u <- rnorm(12, 0, 0.8)
  d$n <- 20
  d$y <- rbinom(nrow(d), d$n, plogis(0.9*(d$dose - (5 + u[d$grp]))))
  a01 <- list(start=c(b1=0, b2=1, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE))
  m0 <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="logistic", family="binomial", start=a01$start, fix=a01$fix)
  m1 <- drmTMB(cbind(y, n - y) ~ dose, rform=b4 ~ 1 | grp, data=d, model="logistic", family="binomial",
               start=a01$start, fix=a01$fix)
  g0 <- glm(cbind(y, n - y) ~ dose, data=d, family=binomial)
  g1 <- lme4::glmer(cbind(y, n - y) ~ dose + (1 | grp), data=d, family=binomial)
  a <- anova(m0, m1)
  ag <- anova(g1, g0)
  expect_equal(a$logLik, ag$logLik, tolerance=1e-6)
  expect_equal(a$Chisq[2], ag$Chisq[2], tolerance=1e-5)
  expect_equal(a[["Pr(>Chisq)"]][2], ag[["Pr(>Chisq)"]][2], tolerance=1e-4)
})

test_that("anova checks its arguments", {
  d <- simgauss("loglogistic")
  f4 <- fit4(y ~ dose, data=d, model="loglogistic")
  f5 <- drmTMB(y ~ dose, data=d, model="loglogistic")
  expect_error(anova(f4), "at least two models")
  expect_error(anova(f4, lm(y ~ dose, data=d)), "drmTMB objects")
  d2 <- d
  d2$y[1] <- d2$y[1] + 1
  expect_error(anova(f4, fit4(y ~ dose, data=d2, model="loglogistic")), "same data")
  # a larger model with a lower log-likelihood
  f5$opt$objective <- f4$opt$objective + 1
  expect_warning(a <- anova(f4, f5), "check the convergence")
  expect_equal(a[["Pr(>Chisq)"]][2], 1)
})
