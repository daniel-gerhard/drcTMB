# quantal data: background 0.1, loglogistic increase
simmaquantal <- function(seed=2){
  set.seed(seed)
  d <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16), each=3), n=20)
  d$y <- rbinom(nrow(d), d$n, 0.1 + 0.9*plogis(2*(log(pmax(d$dose, 1e-9)) - log(4))))
  d
}

quantalfits <- function(d){
  q <- function(m) drmTMB(cbind(y, n - y) ~ dose, data=d, family="binomial", model=m, start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  list(LL=q("loglogistic"), W1=q("weibull1"), W2=q("weibull2"))
}

test_that("the averaged curve BMD is the root of the AIC-weighted mean curve", {
  fs <- quantalfits(simmaquantal())
  ma <- bmdma(LL=fs$LL, W1=fs$W1, W2=fs$W2, bmr=c(0.05, 0.1), nsim=0)
  aic <- vapply(fs, AIC, numeric(1))
  w <- exp(-aic/2)/sum(exp(-aic/2))
  expect_equal(ma$models$weight, unname(w), tolerance=1e-10)
  m <- function(d) sum(w * vapply(fs, function(f) unname(predict(f, newdata=data.frame(dose=d))), numeric(1)))
  for (i in 1:2){
    bmr <- c(0.05, 0.1)[i]
    target <- m(0) + bmr*(1 - m(0))
    r <- uniroot(function(d) m(d) - target, c(0.01, 16), tol=1e-12)$root
    expect_equal(unname(coef(ma)[i]), r, tolerance=1e-6)
  }
  # the BMDs of the models are those of bmd
  expect_equal(unname(coef(ma$bmds$W1)), unname(coef(bmd(fs$W1, c(0.05, 0.1)))))
})

test_that("averaged BMDs and Buckland intervals", {
  fs <- quantalfits(simmaquantal())
  ma <- bmdma(LL=fs$LL, W1=fs$W1, W2=fs$W2, bmr=0.1, average="bmd", nsim=0)
  w <- ma$models$weight
  b <- lapply(fs, bmd, bmr=0.1)
  d <- vapply(b, function(x) unname(coef(x)), numeric(1))
  s <- vapply(b, function(x) sqrt(vcov(x)[1, 1]), numeric(1))
  est <- sum(w*d)
  expect_equal(unname(coef(ma)), est)
  se <- sum(w*sqrt(s^2 + (d - est)^2))
  expect_equal(as.vector(confint(ma)), est + c(-1, 1)*qnorm(0.95)*se)
  expect_equal(as.vector(confint(ma, level=0.95, method="buckland")), est + c(-1, 1)*qnorm(0.975)*se)
  expect_error(confint(bmdma(fs$LL, fs$W1, bmr=0.1, nsim=0)), "average=\"bmd\"")
})

test_that("averages of a single model are its BMD", {
  data(ryegrass, package="drc")
  f <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  b <- coef(bmd(f, 0.1, def="relative"))
  for (a in c("curve", "bmd")){
    ma <- bmdma(f, f, bmr=0.1, def="relative", average=a, nsim=0)
    expect_equal(ma$models$weight, c(0.5, 0.5))
    expect_equal(unname(coef(ma)), unname(b), tolerance=1e-8)
  }
  expect_equal(rownames(ma$models), c("f", "f.1"))
  # sd: weighted mean of the standard deviations of gaussian fits
  ma <- bmdma(a=f, b=f, bmr=1, def="sd", nsim=0)
  expect_equal(unname(coef(ma)), unname(coef(bmd(f, 1, def="sd"))), tolerance=1e-8)
})

test_that("marginal BMDs of models with random effects are averaged", {
  d <- simbmdgroups()
  f <- fit4(y ~ dose, rform=b2 + b4 ~ 1 | grp, data=d, model="loglogistic", link=c(b4="log"))
  bm <- bmd(f, 0.1, def="relative")
  ma <- bmdma(f, f, bmr=0.1, def="relative", nsim=0)
  expect_equal(unname(coef(ma)), unname(coef(bm)), tolerance=1e-6)
  ma <- bmdma(f, f, bmr=0.1, def="relative", scale="conditional", nsim=0)
  expect_equal(unname(coef(ma)), unname(coef(bmd(f, 0.1, def="relative", scale="conditional"))), tolerance=1e-6)
  expect_error(bmdma(f, f, bmr=0.1, def="relative", newdata=data.frame(grp=d$grp[1]), nsim=0), "fitted groups")
})

test_that("bootstrap intervals are reproducible and cover the estimate", {
  fs <- quantalfits(simmaquantal())
  ma <- bmdma(LL=fs$LL, W1=fs$W1, bmr=0.1, nsim=40, seed=3, progress=FALSE)
  expect_equal(dim(ma$boot), c(40, 1))
  expect_true(all(ma$source %in% 1:2))
  expect_equal(rowSums(ma$bootweights), rep(1, 40), tolerance=1e-10)
  ci <- confint(ma)
  expect_true(ci[1] < coef(ma) && coef(ma) < ci[2])
  expect_equal(colnames(ci), c("5 %", "95 %"))
  mb <- bmdma(LL=fs$LL, W1=fs$W1, bmr=0.1, nsim=40, seed=3, progress=FALSE)
  expect_equal(mb$boot, ma$boot)
  expect_output(print(ma), "BMDL")
})

test_that("lognormal models get weight 0 for non-positive bootstrap responses", {
  data(ryegrass, package="drc")
  fg <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  fl <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", family="lognormal", start=c(b5=0), fix=c(b5=TRUE))
  ma <- suppressMessages(bmdma(fg, fl, bmr=0.1, def="relative", nsim=30, seed=1, progress=FALSE))
  expect_true(any(ma$bootweights[, 2] == 0, na.rm=TRUE))
  expect_true(all(ma$bootweights[ma$source == 2, 2] > 0, na.rm=TRUE))
})

test_that("bmdma checks its input", {
  fs <- quantalfits(simmaquantal())
  fo <- drmTMB(cbind(y, n - y) ~ dose, data=simmaquantal(seed=5), family="binomial", model="loglogistic", start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  expect_error(bmdma(fs$LL), "at least two")
  expect_error(bmdma(fs$LL, fo, nsim=0), "same doses")
  expect_error(bmdma(fs$LL, lm(y ~ dose, data=simmaquantal()), nsim=0), "drmTMB")
  data(ryegrass, package="drc")
  fg <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  fl <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", family="lognormal", start=c(b5=0), fix=c(b5=TRUE))
  expect_error(bmdma(fg, fl, def="sd", nsim=0), "gaussian families")
  expect_s3_class(bmdma(fg, fl, def="sd", average="bmd", nsim=0), "bmdma")
  expect_error(bmdma(fg, fl, nsim=0), "model fg: the excess risk")
})
