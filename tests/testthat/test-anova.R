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

test_that("F-tests of gaussian models agree with drc", {
  skip_if_not_installed("drc")
  d <- simgauss("loglogistic")
  # lower asymptote fixed at 0 instead of its true value 1
  f3 <- drmTMB(y ~ dose, data=d, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))
  f4 <- fit4(y ~ dose, data=d, model="loglogistic")
  a <- anova(f3, f4, test="F")
  expect_equal(rownames(a), c("f3", "f4"))
  d3 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4(fixed=c(NA, 0, NA, NA)))
  d4 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4())
  ad <- drc:::anova.drc(d3, d4, details=FALSE)
  expect_equal(a$Res.Df, ad[, "ModelDf"])
  expect_equal(a$RSS, ad[, "RSS"], tolerance=1e-5)
  expect_gt(a$F[2], 10)
  expect_equal(a$F[2], ad[2, "F value"], tolerance=1e-4)
  expect_equal(a[["Pr(>F)"]][2], ad[2, "p value"], tolerance=1e-4)

  # several models: the residual mean square of the largest model, as for anova.lm
  f5 <- drmTMB(y ~ dose, data=d, model="loglogistic")
  a3 <- anova(f3, f4, f5, test="F")
  s2 <- a3$RSS[3]/a3$Res.Df[3]
  expect_equal(a3$F[2:3], -diff(a3$RSS)/s2)
  expect_equal(a3[["Pr(>F)"]][3], pf(a3$F[3], 1, a3$Res.Df[3], lower.tail=FALSE))
})

test_that("F-tests are only for gaussian models without random effects", {
  set.seed(3)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16, 32), each=5))
  d$y <- rpois(nrow(d), drcmean("loglogistic", d$dose, c(2, 60, -1.5, 3, 0)))
  p4 <- fit4(y ~ dose, data=d, model="loglogistic", family="poisson")
  p5 <- drmTMB(y ~ dose, data=d, model="loglogistic", family="poisson")
  expect_error(anova(p4, p5, test="F"), "only available for gaussian")
})
