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

test_that("the constant model tests for a dose effect, as noEffect in drc", {
  skip_if_not_installed("drc")
  d <- simgauss("loglogistic")
  f0 <- drmTMB(y ~ dose, data=d, model="constant")
  f4 <- fit4(y ~ dose, data=d, model="loglogistic")
  expect_equal(as.numeric(logLik(f0)), as.numeric(logLik(lm(y ~ 1, data=d))), tolerance=1e-8)
  expect_equal(unname(coef(f0)[1]), mean(d$y), tolerance=1e-6)
  expect_equal(sum(!f0$fix), 1)
  a <- anova(f0, f4)
  ne <- drc::noEffect(drc::drm(y ~ dose, data=d, fct=drc::LL.4()))
  expect_equal(a$Chisq[2], unname(ne[1]), tolerance=1e-5)
  expect_equal(a$Df[2], unname(ne[2]))
  expect_equal(a[["Pr(>Chisq)"]][2], unname(ne[3]), tolerance=1e-5)
  # the same mean at all doses, also at dose zero
  p <- predict(f0, newdata=data.frame(dose=c(0, 1, 100)))
  expect_equal(p, rep(mean(d$y), 3), tolerance=1e-6)
  expect_error(ed(f0), "no effective doses")
  expect_error(drmTMB(y ~ dose, rform=b4 ~ 1 | dose, data=d, model="constant"), "does not use b4")
})

test_that("the constant model agrees with glm and lme4", {
  skip_if_not_installed("lme4")
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=8, nrep=1, seed=4)
  u <- rnorm(8, 0, 0.6)
  d$n <- 20
  d$y <- rbinom(nrow(d), d$n, plogis(0.3 + u[d$grp]))
  b0 <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="constant", family="binomial", link="logit")
  g0 <- glm(cbind(y, n - y) ~ 1, data=d, family=binomial)
  expect_equal(as.numeric(logLik(b0)), as.numeric(logLik(g0)), tolerance=1e-8)
  expect_equal(unname(coef(b0)[1]), unname(coef(g0)), tolerance=1e-6)
  r0 <- drmTMB(cbind(y, n - y) ~ dose, rform=b1 ~ 1 | grp, data=d, model="constant", family="binomial", link="logit")
  l0 <- lme4::glmer(cbind(y, n - y) ~ 1 + (1 | grp), data=d, family=binomial, nAGQ=1)
  expect_equal(as.numeric(logLik(r0)), as.numeric(logLik(l0)), tolerance=1e-5)
  expect_equal(unname(coef(r0)[1]), unname(lme4::fixef(l0)), tolerance=1e-4)
})

test_that("parametric bootstrap likelihood ratio tests", {
  d <- simgauss("loglogistic", nrep=2)
  f0 <- drmTMB(y ~ dose, data=d, model="constant")
  f3 <- drmTMB(y ~ dose, data=d, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))
  f4 <- fit4(y ~ dose, data=d, model="loglogistic")
  a <- anova(f3, f4, f0, test="boot", nsim=19, seed=1, progress=FALSE)
  ac <- anova(f3, f4, f0)
  expect_s3_class(a, "anova")
  expect_equal(rownames(a), c("f0", "f3", "f4"))
  expect_equal(a[names(ac)], ac, ignore_attr=TRUE)
  expect_true(is.na(a[["Pr(boot)"]][1]))
  lr <- attr(a, "boot")
  expect_named(lr, c("f0 vs f3", "f3 vs f4"))
  expect_true(all(vapply(lr, length, integer(1)) == 19))
  expect_true(all(unlist(lr) >= 0))
  # p-value from the simulated statistics; a strong dose effect exceeds all of them
  expect_equal(a[["Pr(boot)"]][2:3], vapply(2:3, function(i) (1 + sum(lr[[i - 1]] >= a$Chisq[i]))/20, numeric(1)))
  expect_equal(a[["Pr(boot)"]][2], 1/20)
  # reproducible with a seed
  expect_equal(anova(f3, f4, test="boot", nsim=5, seed=2, progress=FALSE)[["Pr(boot)"]],
               anova(f3, f4, test="boot", nsim=5, seed=2, progress=FALSE)[["Pr(boot)"]])
})

test_that("bootstrap tests of a random effect on the boundary", {
  d <- simgroups(ng=6, nrep=1, seed=4)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + rnorm(nrow(d), 0, 0.4)
  r0 <- fit4(y ~ dose, data=d, model="loglogistic")
  r1 <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  a <- anova(r0, r1, test="boot", nsim=19, seed=1, progress=FALSE)
  lr <- attr(a, "boot")[[1]]
  # a point mass at zero: the standard deviation is estimated at zero for many simulations
  expect_gt(sum(lr < 1e-6), 0)
  expect_true(a[["Pr(boot)"]][2] >= 1/20 && a[["Pr(boot)"]][2] <= 1)
})
