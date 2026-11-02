# Threshold model nec, F(x) = 1 - exp(-(b3 (x - b4)_+)^exp(b5)). With b5 = 0, drc's NEC.4 maps to
# b1 = d, b2 = c, b3 = b, b4 = e.

skip_if_not_installed("drc")
skip_if_not_installed("drcData")

# drc coefficients named b, c, d, e
drccoef <- function(fit){
  p <- coef(fit)
  names(p) <- sub(":.*", "", names(p))
  p
}

# free and fixed coefficients b1-b5 of a fit
bfull <- function(fit){
  b <- fit$start
  b[!fit$fix] <- coef(fit)[seq_len(sum(!fit$fix))]
  b
}

necmean <- function(x, b) drcTMB:::drcmeanR(15, x, b[1], b[2], b[3], b[4], b[5])

simnec <- function(b5=0, seed=1){
  set.seed(seed)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 16), each=3))
  d$y <- necmean(d$dose, c(8, 1, 0.4, 1.8, b5)) + rnorm(nrow(d), 0, 0.3)
  d
}

data(ryegrass, package="drcData")

test_that("gaussian nec model agrees with drc NEC.4", {
  f1 <- drmTMB(rootl ~ conc, data=ryegrass, model="nec", start=c(b5=0), fix=c(b5=TRUE))
  # drc's default starting values end in a local optimum with a step at e = 3.4; from e = 1.6 it
  # finds the optimum of the starting values of drcTMB
  f2 <- drc::drm(rootl ~ conc, data=ryegrass, fct=drc::NEC.4(), start=c(0.5, 0.7, 7.4, 1.6), control=drc::drmc(relTol=1e-12))
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(p[c("d", "c", "b", "e")]), tolerance=1e-4)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-6)
  # small-sample standard errors, as drc
  se2 <- sqrt(diag(vcov(f2)))
  expect_equal(unname(sqrt(diag(vcov(f1)))[1:4]), unname(se2[c(3, 2, 1, 4)]), tolerance=1e-3)
  f3 <- drc::drm(rootl ~ conc, data=ryegrass, fct=drc::NEC.4())
  expect_lt(f1$opt$objective, -as.numeric(logLik(f3)) - 1)
})

test_that("nec models with fixed asymptotes agree with drc NEC.3 and binomial NEC.3", {
  f1 <- drmTMB(rootl ~ conc, data=ryegrass, model="nec", start=c(b2=0, b5=0), fix=c(b2=TRUE, b5=TRUE))
  f2 <- drc::drm(rootl ~ conc, data=ryegrass, fct=drc::NEC.3(), start=c(0.5, 7.4, 1.6), control=drc::drmc(relTol=1e-12))
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(p[c("d", "b", "e")]), tolerance=1e-4)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-6)

  set.seed(1)
  q <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16), each=3), n=20)
  q$y <- rbinom(nrow(q), q$n, necmean(q$dose, c(0.9, 0, 0.3, 2.5, 0)))
  f3 <- drmTMB(cbind(y, n - y) ~ dose, data=q, family="binomial", model="nec", start=c(b2=0, b5=0), fix=c(b2=TRUE, b5=TRUE))
  f4 <- drc::drm(y/n ~ dose, weights=n, data=q, type="binomial", fct=drc::NEC.3(), start=c(0.3, 0.9, 2.5))
  p <- drccoef(f4)
  expect_equal(unname(bcoef(f3)), unname(p[c("d", "b", "e")]), tolerance=1e-4)
  expect_equal(f3$opt$objective, -as.numeric(logLik(f4)), tolerance=1e-6)
  se4 <- sqrt(diag(vcov(f4)))
  expect_equal(unname(sqrt(diag(vcov(f3)))), unname(se4[c(2, 1, 3)]), tolerance=1e-3)
})

test_that("sigmoidal nec model maximises the likelihood and anova tests the shape", {
  d <- simnec(b5=0.7)
  f1 <- drmTMB(y ~ dose, data=d, model="nec")
  f0 <- drmTMB(y ~ dose, data=d, model="nec", start=c(b5=0), fix=c(b5=TRUE))
  # least squares from the drcTMB estimates by a different optimiser
  rss <- function(b) sum((d$y - necmean(d$dose, b))^2)
  o <- optim(bcoef(f1), rss, method="BFGS", control=list(reltol=1e-14, maxit=1000))
  expect_equal(unname(bcoef(f1)), unname(o$par), tolerance=1e-4)
  expect_lte(rss(bcoef(f1)), o$value + 1e-8)
  expect_equal(unname(bcoef(f1)), c(8, 1, 0.4, 1.8, 0.7), tolerance=0.15)
  a <- anova(f0, f1)
  expect_lt(a[2, "Pr(>Chisq)"], 0.05)
  expect_output(print(summary(f1)), "The shape parameter b5 is on the log scale.")
  # b4link="log" gives the same fit without random effects on b4
  f2 <- drmTMB(y ~ dose, data=d, model="nec", b4link="log")
  expect_equal(f2$opt$objective, f1$opt$objective, tolerance=1e-6)
  expect_equal(exp(unname(bcoef(f2)[4])), unname(bcoef(f1)[4]), tolerance=1e-4)
})

test_that("effective doses of the nec model have the delta method standard errors", {
  skip_if_not_installed("numDeriv")
  f <- drmTMB(y ~ dose, data=simnec(b5=0.7), model="nec")
  b <- bfull(f)
  e <- ed(f, c(10, 50, 90))
  edR <- function(b, p) b[4] + (-log(1 - p/100))^exp(-b[5]) / b[3]
  expect_equal(unname(coef(e)), sapply(c(10, 50, 90), function(p) edR(b, p)))
  # the curve reaches the response of the effective doses at their doses
  expect_equal(unname(necmean(coef(e), b)), b[1] + (b[2] - b[1]) * c(0.1, 0.5, 0.9))
  g <- t(sapply(c(10, 50, 90), function(p) numDeriv::grad(edR, b, p=p)))
  expect_equal(unname(vcov(e)), unname(g %*% vcov(f)[1:5, 1:5] %*% t(g)), tolerance=1e-6)
  # profile intervals contain the estimates
  ci <- confint(e, parm=1:2)
  expect_true(all(ci[, 1] < coef(e)[1:2] & ci[, 2] > coef(e)[1:2]))
})

test_that("nec model with a random asymptote maximises the exact likelihood of nlme", {
  skip_if_not_installed("nlme")
  set.seed(4)
  d <- expand.grid(rep=1:2, dose=c(0, 0.5, 1, 2, 4, 8, 16), grp=factor(1:10))
  u <- rnorm(10, 0, 0.5)
  d$y <- necmean(d$dose, c(8, 1, 0.4, 1.5, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f1 <- drmTMB(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="nec", start=c(b5=0), fix=c(b5=TRUE))
  n1 <- nlme::nlme(y ~ b1 + (b2 - b1)*(1 - exp(-b3*pmax(dose - b4, 0))), fixed=b1 + b2 + b3 + b4 ~ 1, random=b1 ~ 1,
                   groups=~grp, data=d, start=setNames(bcoef(f1), paste("b", 1:4, sep="")), method="ML")
  np <- c(nlme::fixef(n1), log(as.numeric(nlme::VarCorr(n1)[1, 2])), log(n1$sigma))
  expect_equal(as.numeric(f1$obj$fn(np)), -as.numeric(logLik(n1)), tolerance=1e-6)
  expect_lte(f1$opt$objective, -as.numeric(logLik(n1)) + 1e-6)
  expect_equal(unname(bcoef(f1)), unname(nlme::fixef(n1)), tolerance=1e-3)
})

test_that("nec model warns for thresholds at a dose or above the highest dose", {
  f <- drmTMB(rootl ~ conc, data=ryegrass, model="nec", start=c(b5=0), fix=c(b5=TRUE))
  expect_silent(drcTMB:::necwarn(f))
  f$opt$par[4] <- 1.88
  expect_warning(drcTMB:::necwarn(f), "at the dose 1.88")
  f$opt$par[4] <- 40
  expect_warning(drcTMB:::necwarn(f), "above the highest dose")
})

test_that("starting values of the nec model keep given values", {
  sv <- drmTMBstart(rootl ~ conc, data=ryegrass, model="nec", start=c(NA, NA, NA, 2, 0), fix=c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(sv[[4]], 2)
  expect_true(sv[[3]] > 0)
  # b4link="log": positive thresholds only
  sv <- drmTMBstart(rootl ~ conc, data=ryegrass, model="nec", b4link="log")
  expect_true(is.finite(sv[[4]]))
})
