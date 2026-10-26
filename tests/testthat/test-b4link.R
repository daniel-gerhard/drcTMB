# grouped loglogistic data with b4 varying by a factor between groups
simlogb4 <- function(ng=12, sd=0.5, seed=3){
  d <- simgroups(dose=c(0.5, 1, 2, 4, 8, 16, 32), ng=ng, seed=seed)
  u <- rnorm(ng, 0, sd)
  d$y <- drcTMB:::drcmeanR(2, d$dose, 1, 8, -1.5, 3*exp(u[d$grp]), 0) + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("random effects on log b4 give the logistic model of log dose", {
  d <- simlogb4()
  d$ldose <- log(d$dose)
  fl <- fit4(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="loglogistic", b4link="log")
  fg <- fit4(y ~ ldose, rform=b4 ~ 1 | grp, data=d, model="logistic")
  expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(fg)), tolerance=1e-6)
  expect_equal(coef(fl), coef(fg), tolerance=1e-4)
  expect_equal(predict(fl), predict(fg), tolerance=1e-6)
  expect_equal(coef(ed(fl, random=FALSE)), exp(coef(ed(fg, random=FALSE))), tolerance=1e-4)
  nd <- data.frame(dose=c(1, 3, 10), grp=factor(c(1, 2, 99)))
  nd$ldose <- log(nd$dose)
  expect_equal(predict(fl, nd, interval="prediction"), predict(fg, nd, interval="prediction"), tolerance=1e-4)
})

test_that("predict has the same standard errors and intervals with both links of b4", {
  d <- simgauss("loglogistic")
  fi <- fit4(y ~ dose, data=d, model="loglogistic")
  fl <- fit4(y ~ dose, data=d, model="loglogistic", b4link="log")
  nd <- data.frame(dose=c(0, 1, 3, 10))
  expect_equal(predict(fl, nd, interval="confidence"), predict(fi, nd, interval="confidence"), tolerance=1e-4)
  expect_equal(predict(fl, nd, interval="prediction"), predict(fi, nd, interval="prediction"), tolerance=1e-4)
})

test_that("without random effects both links of b4 give the same fit and effective doses", {
  d <- simgauss("loglogistic")
  for (m in c("loglogistic", "weibull1", "lognormal", "michaelismenten", "asymptoticregression")){
    fi <- drmTMB(y ~ dose, data=d, model=m)
    fl <- drmTMB(y ~ dose, data=d, model=m, b4link="log")
    expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(fi)), tolerance=1e-8, info=m)
    expect_equal(exp(bcoef(fl)[4]), bcoef(fi)[4], tolerance=1e-5, info=m)
    expect_equal(predict(fl), predict(fi), tolerance=1e-5, info=m)
    expect_equal(coef(ed(fl, c(10, 50))), coef(ed(fi, c(10, 50))), tolerance=1e-5, info=m)
    expect_equal(vcov(ed(fl, c(10, 50))), vcov(ed(fi, c(10, 50))), tolerance=1e-4, info=m)
  }
  # grid search over the default range of log b4
  fg <- drmTMB(y ~ dose, data=d, model="loglogistic", b4link="log", gsn=5)
  expect_equal(as.numeric(logLik(fg)), as.numeric(logLik(drmTMB(y ~ dose, data=d, model="loglogistic"))), tolerance=1e-8)
})

test_that("effective doses with log b4 have the same profile and bootstrap intervals", {
  d <- simgauss("loglogistic")
  ei <- ed(fit4(y ~ dose, data=d, model="loglogistic"))
  el <- ed(fit4(y ~ dose, data=d, model="loglogistic", b4link="log"))
  expect_equal(confint(profile(el, progress=FALSE)), confint(profile(ei, progress=FALSE)), tolerance=1e-3)
  ci <- suppressMessages(confint(ei, method="boot", nsim=20, seed=1, progress=FALSE))
  cl <- suppressMessages(confint(el, method="boot", nsim=20, seed=1, progress=FALSE))
  expect_equal(cl, ci, tolerance=1e-4)
})

test_that("print and summary show the log link of b4", {
  d <- simlogb4(ng=8)
  f <- fit4(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="loglogistic", b4link="log")
  expect_true(any(grepl("log link for b4", capture.output(print(f)))))
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("log link for b4", out)))
  expect_true(any(grepl("coefficients of b4 are on the log scale", out)))
})

test_that("b4link is checked", {
  d <- simgauss("logistic")
  expect_error(drmTMB(y ~ dose, data=d, model="logistic", b4link="log"), "positive b4")
  expect_error(drmTMB(y ~ dose, data=d, model="gompertz", b4link="log"), "positive b4")
  expect_error(drmTMB(y ~ dose, data=d, model="loglogistic", b4link="logit"), "Unknown b4link")
})
