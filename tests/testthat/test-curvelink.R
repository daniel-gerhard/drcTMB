# link of the whole curve (curvelink): the mean, median or probability is the inverse link of the curve

library(drcData)
data(ryegrass)
data(daphnids)
ry <- ryegrass
ry$ly <- log(ry$rootl)
d24 <- daphnids[daphnids$time == "24h", ]

test_that("lognormal with the log curve link is the gaussian model of log(y)", {
  fl <- drmTMB(rootl ~ conc, data=ry, model="weibull2", family="lognormal", curvelink="log", start=c(b5=0), fix=c(b5=TRUE))
  fg <- drmTMB(ly ~ conc, data=ry, model="weibull2", family="gaussian", start=c(b5=0), fix=c(b5=TRUE), control=list(smallsample=FALSE))
  expect_equal(coef(fl), coef(fg), tolerance=1e-6)
  # density of y: the Jacobian of log(y)
  expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(fg)) - sum(ry$ly), tolerance=1e-8)
  expect_equal(predict(fl), exp(predict(fg)), tolerance=1e-6)
  # effective doses refer to the curve
  expect_equal(coef(ed(fl, c(10, 50))), coef(ed(fg, c(10, 50))), tolerance=1e-6)
  # absolute effective doses take responses on the scale of the response
  med <- function(x) predict(fl, newdata=data.frame(conc=x))
  expect_equal(unname(coef(ed(fl, 4, type="absolute"))), uniroot(function(x) med(x) - 4, c(1e-6, 30), tol=1e-10)$root, tolerance=1e-6)
  expect_equal(fl$curvelink, "log")
  expect_output(print(fl), "log link of the curve, the mean \\(median\\) is exp\\(curve\\)")
  expect_output(print(summary(fl)), "log link of the curve")
})

test_that("probit and logit curve links of the power model, and the gaussian and poisson log link, are glms", {
  for (lk in c("probit", "logit")){
    fb <- drmTMB(cbind(no, total - no) ~ dose, data=d24, family="binomial", model="power", curvelink=lk, start=c(b3=1), fix=c(b3=TRUE))
    gb <- glm(cbind(no, total - no) ~ dose, data=d24, family=binomial(lk))
    expect_equal(unname(coef(fb)), unname(coef(gb)), tolerance=1e-5)
    expect_equal(as.numeric(logLik(fb)), as.numeric(logLik(gb)), tolerance=1e-8)
    expect_equal(predict(fb, newdata=data.frame(dose=c(0, 1000))), unname(predict(gb, newdata=data.frame(dose=c(0, 1000)), type="response")), tolerance=1e-5)
  }
  fq <- drmTMB(rootl ~ conc, data=ry, model="power", curvelink="log", start=c(b3=1), fix=c(b3=TRUE), control=list(smallsample=FALSE))
  gq <- glm(rootl ~ conc, data=ry, family=gaussian("log"))
  expect_equal(unname(coef(fq)[1:2]), unname(coef(gq)), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fq)), as.numeric(logLik(gq)), tolerance=1e-8)
  set.seed(1)
  dp <- data.frame(x=rep(0:5, each=4))
  dp$y <- rpois(nrow(dp), exp(2 - 0.3*dp$x))
  fpo <- drmTMB(y ~ x, data=dp, family="poisson", model="power", curvelink="log", start=c(b3=1), fix=c(b3=TRUE))
  gpo <- glm(y ~ x, data=dp, family=poisson)
  expect_equal(unname(coef(fpo)), unname(coef(gpo)), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fpo)), as.numeric(logLik(gpo)), tolerance=1e-8)
})

test_that("benchmark doses on the scale of the response, and the BMD parameterisation", {
  fl <- drmTMB(rootl ~ conc, data=ry, model="weibull2", family="lognormal", curvelink="log", start=c(b5=0), fix=c(b5=TRUE))
  med <- function(x) predict(fl, newdata=data.frame(conc=x))
  b <- bmd(fl, 0.05, "relative")
  expect_equal(unname(coef(b)), uniroot(function(x) med(x) - 0.95*med(0), c(1e-6, 30), tol=1e-10)$root, tolerance=1e-5)
  fb <- drmTMB(rootl ~ conc, data=ry, model="weibull2", family="lognormal", curvelink="log", start=c(b5=0), fix=c(b5=TRUE),
               bmdpar=list(bmr=0.05, def="relative"))
  expect_equal(exp(coef(fb)[["log_bmd"]]), unname(coef(b)), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fb)), as.numeric(logLik(fl)), tolerance=1e-7)

  fp <- drmTMB(cbind(no, total - no) ~ dose, data=d24, family="binomial", model="power", curvelink="probit")
  pr <- function(x) predict(fp, newdata=data.frame(dose=x))
  p0 <- pr(0)
  bp <- bmd(fp, 0.1, "excess")
  expect_equal(unname(coef(bp)), uniroot(function(x) pr(x) - (p0 + 0.1*(1 - p0)), c(1e-6, 5000), tol=1e-10)$root, tolerance=1e-5)
  fpb <- drmTMB(cbind(no, total - no) ~ dose, data=d24, family="binomial", model="power", curvelink="probit",
                bmdpar=list(bmr=0.1, def="excess"))
  expect_equal(exp(coef(fpb)[["log_bmd"]]), unname(coef(bp)), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fpb)), as.numeric(logLik(fp)), tolerance=1e-7)
})

test_that("marginal BMD with random effects on b2 and the log curve link", {
  set.seed(3)
  dd <- expand.grid(dose=c(0, 0.5, 1, 2, 4, 8, 16), grp=factor(1:8), rep=1:2)
  u <- rnorm(8, 0, 0.3)
  dd$y <- exp(0.2 + (2 + u[dd$grp] - 0.2)/(1 + (dd$dose/3)^1.5) + rnorm(nrow(dd), 0, 0.15))
  fr <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=dd, model="loglogistic", family="lognormal", curvelink="log",
               start=c(b5=0), fix=c(b5=TRUE))
  bm <- bmd(fr, 0.1, "relative")
  bc <- bmd(fr, 0.1, "relative", scale="conditional")
  # with a curve link, random effects on an asymptote do not enter the median linearly
  expect_true(abs(coef(bm) - coef(bc)) > 1e-3)
  cf <- coef(fr)
  s2 <- exp(cf[["log_sd:b2|grp"]])
  z <- qnorm(ppoints(20000))
  mm <- function(x) mean(exp(drcmeanR(modelindex("loglogistic"), x, cf[[1]], cf[[2]] + s2*z, cf[[3]], cf[[4]], 1)))
  expect_equal(unname(coef(bm)), uniroot(function(x) mm(x) - 0.9*mm(0), c(1e-4, 20), tol=1e-10)$root, tolerance=1e-4)
})

test_that("simulation on the scale of the response; checks of the curve link", {
  fl <- drmTMB(rootl ~ conc, data=ry, model="weibull2", family="lognormal", curvelink="log", start=c(b5=0), fix=c(b5=TRUE))
  s <- simulate(fl, nsim=2000, seed=1)
  expect_equal(apply(log(as.matrix(s)), 1, mean), log(predict(fl)), tolerance=0.02)

  expect_error(drmTMB(rootl ~ conc, data=ry, model="weibull2", curvelink="probit"), "not available for the gaussian family")
  expect_error(drmTMB(rootl ~ conc, data=ry, model="weibull2", curvelink="log", parlink=c(b1="log")), "identity parlink")
  expect_error(drmTMB(rootl ~ conc, data=ry, model="weibull2", curvelink=c("log", "log")), "single link")
  expect_error(ed(fl, 50, residual=TRUE), "not available for fits with a curvelink")
  data(lettuce)
  fh <- drmTMB(weight ~ conc, data=lettuce, model="braincousens", curvelink="log")
  expect_error(peak(fh), "not available for fits with a curvelink")
  # beta-binomial with the probit link of the curve, for overdispersed litters (the daphnids are not
  # overdispersed): the parameters are recovered
  set.seed(5)
  dl <- data.frame(dose=rep(c(0, 5, 10, 20, 40), each=30), n=10)
  p <- pnorm(-1.5 + 0.05*dl$dose)
  pl <- rbeta(nrow(dl), p*10, (1 - p)*10)
  dl$y <- rbinom(nrow(dl), dl$n, pl)
  fbb <- drmTMB(cbind(y, n - y) ~ dose, data=dl, family="betabinomial", model="power", curvelink="probit", start=c(b3=1), fix=c(b3=TRUE))
  expect_equal(fbb$opt$convergence, 0)
  expect_equal(unname(coef(fbb)[1:2]), c(-1.5, 0.05), tolerance=0.15)
})
