# substrate inhibition model haldane: f(x) = b1 + (b2 - b1) x / (b4 + x + b5 x^2)

simhaldane <- function(Ki=20, seed=1){
  set.seed(seed)
  d <- data.frame(S=rep(c(0, 0.5, 1, 2, 4, 8, 16, 32, 64), each=3))
  d$v <- 10 * d$S / (2 + d$S + d$S^2/Ki) + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("haldane agrees with nls, with b5 = 1/Ki", {
  d <- simhaldane()
  f <- drmTMB(v ~ S, data=d, model="haldane", start=c(b1=0), fix=c(b1=TRUE))
  n <- nls(v ~ Vm*S/(Km + S + S^2/Ki), data=d, start=list(Vm=10, Km=2, Ki=20))
  cn <- coef(n)
  expect_equal(unname(coef(f)[1:3]), unname(c(cn[["Vm"]], cn[["Km"]], 1/cn[["Ki"]])), tolerance=1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(n)), tolerance=1e-6)
  expect_equal(f$lower[3], 0)
  # mean of the template and of R
  nd <- data.frame(S=c(0, 1, 6, 50))
  expect_equal(predict(f, nd, se.fit=TRUE)$fit, predict(f, nd))
  expect_equal(predict(f, nd), unname(predict(n, nd)), tolerance=1e-5)
  expect_error(ed(f), "no effective doses")
})

test_that("the peak of haldane is sqrt(Km Ki), with delta method standard errors", {
  d <- simhaldane()
  f <- drmTMB(v ~ S, data=d, model="haldane", start=c(b1=0), fix=c(b1=TRUE))
  p <- peak(f)
  b <- coef(f)[1:3]
  pk <- sqrt(b[[2]]/b[[3]])
  expect_equal(unname(coef(p)), c(pk, b[[1]]*pk/(b[[2]] + pk + b[[3]]*pk^2)))
  g <- numDeriv::jacobian(function(b) sqrt(b[2]/b[3]), b)
  expect_equal(vcov(p)[1, 1], (g %*% vcov(f)[1:3, 1:3] %*% t(g))[1, 1], tolerance=1e-6)
  # profile intervals of the dose of the peak, close to the Wald intervals
  ci <- suppressMessages(confint(p))
  cw <- confint(p, method="wald")
  expect_true(ci[1, 1] < pk && ci[1, 2] > pk)
  expect_equal(ci[1, ], cw[1, ], tolerance=0.01)
})

test_that("haldane without inhibition is michaelismenten", {
  d <- simhaldane(Ki=Inf)
  fm <- drmTMB(v ~ S, data=d, model="michaelismenten")
  f0 <- drmTMB(v ~ S, data=d, model="haldane", start=c(b5=0), fix=c(b5=TRUE))
  expect_equal(coef(f0), coef(fm), tolerance=1e-5)
  expect_equal(as.numeric(logLik(f0)), as.numeric(logLik(fm)), tolerance=1e-8)
  # anova tests for inhibition
  f <- drmTMB(v ~ S, data=simhaldane(), model="haldane")
  fm <- drmTMB(v ~ S, data=simhaldane(), model="michaelismenten")
  a <- anova(fm, f)
  expect_equal(a$Df[2], 1)
  expect_lt(a[2, "Pr(>Chisq)"], 1e-10)
  # without inhibition, there is no peak
  expect_warning(p <- peak(f0), "no peak")
  expect_true(all(is.na(coef(p))))
  expect_error(peak(fm), "hormesis and substrate inhibition models")
})

test_that("haldane with random effects on the half-saturation constant", {
  set.seed(3)
  d <- expand.grid(rep=1:2, S=c(0, 0.5, 1, 2, 4, 8, 16, 32, 64), grp=factor(1:8))
  u <- rnorm(8, 0, 0.3)
  d$v <- 10 * d$S / (2*exp(u[d$grp]) + d$S + d$S^2/20) + rnorm(nrow(d), 0, 0.3)
  f <- drmTMB(v ~ S, rform=b4 ~ 1 | grp, data=d, model="haldane", link=c(b4="log"), start=c(b1=0), fix=c(b1=TRUE))
  expect_equal(f$opt$convergence, 0)
  expect_equal(unname(exp(coef(f)[2])), 2, tolerance=0.2)
  expect_equal(unname(1/coef(f)[3]), 20, tolerance=0.2)
  # peaks of fitted groups: sqrt(b4/b5) with the group's b4
  nd <- data.frame(grp=c("1", "2"))
  p <- peak(f, newdata=nd)
  b <- predict(f, nd, type="parameters")
  expect_equal(unname(coef(p)[1:2]), unname(sqrt(b[, "b4"]/b[, "b5"])))
})
