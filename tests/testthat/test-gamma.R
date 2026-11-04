# Gamma model, F(x) = pgamma(x/b4, shape b3). drc's gammadr maps to b1 = c, b2 = d, b3 = e,
# b4 = 1/b (drc estimates the rate b); drc has no effective doses for gammadr.

skip_if_not_installed("drc")

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

simgamma <- function(seed=1){
  set.seed(seed)
  d <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16, 32), each=4))
  d$y <- drcTMB:::drcmeanR(13, d$dose, 1, 8, 2.5, 3, 0) + rnorm(nrow(d), 0, 0.4)
  d
}

test_that("gaussian gamma model agrees with drc", {
  d <- simgamma()
  f1 <- drmTMB(y ~ dose, data=d, model="gamma")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::gammadr(), control=drc::drmc(relTol=1e-12))
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(c(p["c"], p["d"], p["e"], 1/p["b"])), tolerance=1e-4)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-6)
  # drc's vcov in the order b, c, d, e; the standard error of b4 = 1/b by the delta method
  se2 <- sqrt(diag(vcov(f2)))
  expect_equal(unname(sqrt(diag(vcov(f1)))[1:4]), unname(c(se2[2:4], se2[1]/p[["b"]]^2)), tolerance=1e-3)
})

test_that("quantal gamma model with background risk agrees with drc", {
  set.seed(2)
  q <- data.frame(dose=c(0, 0.5, 1, 2, 4, 8, 16), n=50)
  q$dead <- rbinom(nrow(q), q$n, drcTMB:::drcmeanR(13, q$dose, 0.05, 1, 2, 2, 0))
  f1 <- drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="gamma", start=c(b2=1), fix=c(b2=TRUE))
  f2 <- drc::drm(dead/n ~ dose, weights=n, data=q, type="binomial", fct=drc::gammadr(fixed=c(NA, NA, 1, NA)))
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(c(p["c"], p["e"], 1/p["b"])), tolerance=1e-4)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-6)
  se2 <- sqrt(diag(vcov(f2)))
  expect_equal(unname(sqrt(diag(vcov(f1)))[1:3]), unname(c(se2[2:3], se2[1]/p[["b"]]^2)), tolerance=1e-3)
})

test_that("effective doses of the gamma model have the delta method standard errors", {
  skip_if_not_installed("numDeriv")
  f <- drmTMB(y ~ dose, data=simgamma(), model="gamma")
  b <- bfull(f)
  e <- ed(f, c(10, 50, 90))
  edR <- function(b, p) b[4]*qgamma(p/100, b[3])
  expect_equal(unname(coef(e)), sapply(c(10, 50, 90), function(p) edR(b, p)))
  g <- t(sapply(c(10, 50, 90), function(p) numDeriv::grad(edR, b[1:4], p=p)))
  expect_equal(unname(vcov(e)), unname(g %*% vcov(f)[1:4, 1:4] %*% t(g)), tolerance=1e-6)
  # mean of the template and of R
  nd <- data.frame(dose=c(0, 0.5, 3, 20))
  expect_equal(predict(f, newdata=nd), drcTMB:::drcmeanR(13, nd$dose, b[1], b[2], b[3], b[4], 0))
  # profile intervals of the coefficients and effective doses contain the estimates
  ci <- suppressMessages(confint(e, parm=2, progress=FALSE))
  expect_true(ci[1, 1] < coef(e)[2] && coef(e)[2] < ci[1, 2])
})

test_that("the gamma model with shape 1 is the asymptoticregression model", {
  d <- simgamma()
  f1 <- drmTMB(y ~ dose, data=d, model="gamma", start=c(b3=1), fix=c(b3=TRUE))
  f2 <- drmTMB(y ~ dose, data=d, model="asymptoticregression")
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-6)
  expect_equal(unname(bcoef(f1)), unname(bcoef(f2)), tolerance=1e-4)
  f3 <- drmTMB(y ~ dose, data=d, model="gamma")
  expect_equal(anova(f2, f3)$Df[2], 1)
})

test_that("gamma mixed model recovers the simulated curve", {
  set.seed(4)
  d <- simgroups(dose=c(0, 1, 2, 4, 8, 16, 32), ng=10)
  u4 <- rnorm(10, 0, 0.3)
  d$y <- drcTMB:::drcmeanR(13, d$dose, 1, 8, 2.5, 3*exp(u4[d$grp]), 0) + rnorm(nrow(d), 0, 0.3)
  f <- drmTMB(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="gamma", link=c(b4="log"))
  expect_equal(f$opt$convergence, 0)
  expect_equal(unname(bcoef(f)), c(1, 8, 2.5, log(3)), tolerance=0.1)
  # effective dose of a fitted group: its predicted curve is halfway between the asymptotes there
  e <- ed(f, 50, newdata=data.frame(grp="1"))
  b <- bfull(f)
  expect_equal(predict(f, newdata=data.frame(dose=coef(e), grp="1")), (b[[1]] + b[[2]])/2, tolerance=1e-8)
})
