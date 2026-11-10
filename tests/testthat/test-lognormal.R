# lognormal responses around a decreasing log-logistic median curve
simlognormal <- function(seed=4, sd=0.2){
  set.seed(seed)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=4))
  d$y <- (2 + 8/(1 + exp(1.5*(log(pmax(d$dose, 1e-12)) - log(3))))) * exp(rnorm(nrow(d), 0, sd))
  d
}
ll4 <- function(x, p) p[1] + (p[2] - p[1])/(1 + exp(-p[3]*(log(x) - log(p[4]))))

test_that("the lognormal family is a nonlinear regression of log(y)", {
  d <- simlognormal()
  f <- fit4(y ~ dose, data=d, family="lognormal", model="loglogistic")
  n <- nls(log(y) ~ log(b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4))))), data=d,
           start=list(b1=2, b2=10, b3=-1.5, b4=3))
  expect_equal(unname(bcoef(f)[1:4]), unname(coef(n)), tolerance=1e-5)
  expect_equal(residpar(f), summary(n)$sigma, tolerance=1e-6)
  expect_equal(df.residual(f), df.residual(n))
  # density of y: the log-likelihood of log(y) minus the Jacobian
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(n)) - sum(log(d$y)), tolerance=1e-8)
  # standard errors from the observed information, with the unbiased residual variance
  rss <- function(p) sum((log(d$y) - log(ll4(d$dose, p)))^2)
  H <- numDeriv::hessian(rss, coef(n)) / (2*summary(n)$sigma^2)
  expect_equal(unname(sqrt(diag(vcov(f)))[1:4]), sqrt(diag(solve(H))), tolerance=1e-5)
  # residuals of log(y)
  expect_equal(unname(residuals(f)), as.vector(residuals(n))/summary(n)$sigma, tolerance=1e-5)
  expect_error(drmTMB(y ~ dose, data=transform(d, y=y - 3), family="lognormal", model="loglogistic"), "positive")
})

test_that("predictions of the lognormal median", {
  d <- simlognormal()
  f <- fit4(y ~ dose, data=d, family="lognormal", model="loglogistic")
  nd <- data.frame(dose=c(0.1, 3, 20))
  pc <- predict(f, newdata=nd, interval="confidence")
  pp <- predict(f, newdata=nd, interval="prediction")
  expect_equal(pc$fit, ll4(nd$dose, bcoef(f)[1:4]), tolerance=1e-8)
  # intervals symmetric on the log scale, with t quantiles and the residual sd of log(y)
  z <- qt(0.975, df.residual(f))
  expect_equal(log(pc$upr/pc$fit), z*pc$se.fit/pc$fit, tolerance=1e-8)
  expect_equal(log(pp$upr/pp$fit), z*sqrt((pp$se.fit/pp$fit)^2 + residpar(f)^2), tolerance=1e-8)
  expect_equal(log(pp$upr/pp$fit), log(pp$fit/pp$lwr), tolerance=1e-8)
  s <- simulate(f, 3, seed=1)
  expect_true(all(unlist(s) > 0))
})

test_that("F-tests of nested lognormal models use the residual sums of squares of log(y)", {
  d <- simlognormal()
  f4 <- fit4(y ~ dose, data=d, family="lognormal", model="loglogistic")
  f3 <- drmTMB(y ~ dose, data=d, family="lognormal", model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))
  a <- anova(f3, f4, test="F")
  n4 <- nls(log(y) ~ log(b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4))))), data=d, start=list(b1=2, b2=10, b3=-1.5, b4=3))
  n3 <- nls(log(y) ~ log(b2/(1 + exp(-b3*(log(dose) - log(b4))))), data=d, start=list(b2=10, b3=-1.5, b4=3))
  an <- anova(n3, n4)
  expect_equal(a$F[2], an$`F value`[2], tolerance=1e-5)
  expect_equal(a$`Pr(>F)`[2], an$`Pr(>F)`[2], tolerance=1e-5)
  expect_error(suppressWarnings(anova(f4, fit4(y ~ dose, data=d, model="loglogistic"), test="F")), "lognormal")
})

test_that("BMDs of lognormal responses", {
  d <- simlognormal()
  f <- fit4(y ~ dose, data=d, family="lognormal", model="loglogistic")
  mu0 <- bcoef(f)[[2]]
  # relative deviation of the median, which is also that of the mean
  b <- bmd(f, 0.1, def="relative")
  expect_equal(unname(coef(b)), unname(coef(ed(f, 0.9*mu0, type="absolute"))), tolerance=1e-8)
  # sd: standard deviations of log(y)
  bs <- bmd(f, 1, def="sd")
  expect_equal(unname(coef(bs)), unname(coef(ed(f, mu0*exp(-residpar(f)), type="absolute"))), tolerance=1e-8)
  pf <- f$sdrl$par.fixed
  val <- function(p){
    pv <- f$start
    pv[!f$fix] <- p[seq_len(sum(!f$fix))]
    bmdvalue(f, bs$bmd[[1]], pv, numeric(0), p[["log_sigma"]])
  }
  g <- numDeriv::grad(val, pf)
  expect_equal(sqrt(vcov(bs)[1, 1]), sqrt(drop(t(g) %*% fixedcov(f) %*% g)), tolerance=1e-6)
  ci <- confint(bs, level=0.9, progress=FALSE)
  expect_true(ci[1, 1] < coef(bs) && ci[1, 2] > coef(bs))
})

test_that("lognormal models with random effects", {
  set.seed(5)
  d <- simgroups(ng=10)
  u2 <- rnorm(10, 0, 0.3)
  u4 <- rnorm(10, 0, 0.4)
  d$y <- with(d, (1 + (8*exp(u2[grp]) - 1)/(1 + exp(1.5*(log(pmax(dose, 1e-12)) - log(3) - u4[grp]))))) * exp(rnorm(nrow(d), 0, 0.15))
  f <- fit4(y ~ dose, rform=b4 ~ 1 | grp, data=d, family="lognormal", model="loglogistic", link=c(b4="log"))
  expect_null(df.residual(f))
  expect_equal(f$opt$convergence, 0)
  bm <- bmd(f, 0.1, def="relative")
  bc <- bmd(f, 0.1, def="relative", scale="conditional")
  expect_true(bm$integrated[[1]])
  expect_true(coef(bm) < coef(bc))
  expect_error(bmd(f, 1, def="sd"), "random effects")
})
