# loglogistic data with b4 depending on a treatment factor and a numeric covariate
simcov <- function(seed=5){
  set.seed(seed)
  d <- expand.grid(rep=1:3, dose=c(0, 0.5, 1, 2, 4, 8, 16, 32), trt=c("A", "B"), temp=c(10, 20))
  d$rep <- NULL
  b4 <- 2*ifelse(d$trt == "B", 2, 1)*exp(0.03*(d$temp - 15))
  d$y <- drcTMB:::drcmeanR(2, d$dose, 1, 8, -1.5, b4, 0) + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("ed with newdata gives the effective doses of the linear functions of each row", {
  d <- simcov()
  f <- fit4(y ~ dose, fform=b4 ~ trt + temp, data=d, model="loglogistic", b4link="log")
  # coefficients b1, b2, b3, b4:(Intercept), b4:trtB, b4:temp, b5
  lf <- function(trt, temp) rbind(diag(7)[1:3, ], c(0, 0, 0, 1, trt, temp, 0), diag(7)[7, ])
  el <- ed(f, rep(c(10, 50), each=4), linfct=list(lf(0, 10), lf(1, 10), lf(0, 20), lf(1, 20), lf(0, 10), lf(1, 10), lf(0, 20), lf(1, 20)))
  nd <- expand.grid(trt=c("A", "B"), temp=c(10, 20))
  en <- ed(f, c(10, 50), newdata=nd)
  expect_equal(unname(coef(en)), unname(coef(el)))
  expect_equal(unname(vcov(en)), unname(vcov(el)))
  # every level of respLev for every row, named by level and covariate values
  expect_equal(names(coef(en))[1:4], c("10:trt=A,temp=10", "10:trt=B,temp=10", "10:trt=A,temp=20", "10:trt=B,temp=20"))
  expect_equal(names(coef(en))[5], "50:trt=A,temp=10")
  # row names of newdata, other variables and the dose are ignored
  nr <- data.frame(trt=c("A", "B"), temp=15, dose=99, row.names=c("ctrl", "treated"))
  expect_equal(names(coef(ed(f, newdata=nr))), c("50:ctrl", "50:treated"))
  # profile and bootstrap intervals use the linear functions of the rows
  e1 <- ed(f, newdata=nd[2, ])
  expect_equal(confint(profile(e1, progress=FALSE)), confint(profile(ed(f, linfct=lf(1, 10)), progress=FALSE)),
               tolerance=1e-6, ignore_attr=TRUE)
})

test_that("ed with newdata works without covariates", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  e <- ed(f, newdata=data.frame(dose=1))
  expect_equal(unname(coef(e)), unname(coef(ed(f))))
  expect_equal(names(coef(e)), "50:1")
})

test_that("ed checks newdata", {
  d <- simcov()
  f <- fit4(y ~ dose, fform=b4 ~ trt, data=d, model="loglogistic")
  expect_error(ed(f), "newdata")
  expect_error(ed(f, newdata=data.frame(trt="A"), linfct=diag(6)), "either linfct or newdata")
  expect_error(ed(f, newdata=data.frame(dose=1)), "covariates of fform: trt")
  expect_error(ed(f, newdata=data.frame(trt=NA_character_)), "missing values")
  expect_error(ed(f, newdata=data.frame(trt=character(0))), "at least one row")
})

# grouped loglogistic data with random effects on b1 and log b4
simre <- function(ng=10, seed=3){
  d <- simgroups(dose=c(0, 0.5, 1, 2, 4, 8, 16, 32), ng=ng, seed=seed)
  u4 <- rnorm(ng, 0, 0.4)
  u1 <- rnorm(ng, 0, 0.3)
  d$y <- drcTMB:::drcmeanR(2, d$dose, 1 + u1[d$grp], 8, -1.5, 3*exp(u4[d$grp]), 0) + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("ed with fitted groups in newdata gives the effective doses of the group curves", {
  d <- simre()
  f <- fit4(y ~ dose, rform=b1 + b4 ~ 1 | grp, data=d, model="loglogistic", b4link="log")
  nd <- data.frame(grp=c("2", "7", "new"))
  e <- ed(f, c(20, 50), newdata=nd)
  expect_equal(names(coef(e))[1:3], c("20:grp=2", "20:grp=7", "20:grp=new"))
  expect_equal(unname(e$conditional), rep(c(TRUE, TRUE, FALSE), 2))

  # the dose at which the predicted curve of the group is halfway between its asymptotes
  pr <- function(dose, g) predict(f, newdata=data.frame(dose=dose, grp=g))
  for (g in c("2", "7")){
    lo <- pr(1e-8, g)
    hi <- pr(1e8, g)
    ed50 <- uniroot(function(x) pr(x, g) - (lo + hi)/2, c(0.01, 100), tol=1e-10)$root
    expect_equal(coef(e)[[paste("50:grp=", g, sep="")]], ed50, tolerance=1e-6)
  }

  # standard errors from the joint normal approximation of the coefficients and random effects
  set.seed(1)
  dr <- drcTMB:::jointdraws(f, 4000)
  r <- f$obj$env$random
  sim <- sapply(4:5, function(i) apply(dr, 2, function(par){
    pv <- f$start
    pv[!f$fix] <- par[1:4]
    lp <- as.vector(e$linfct[[i]] %*% pv + e$uweights[[i]] %*% par[r])
    drcTMB:::edfun(2, e$q[[i]], lp[3], exp(lp[4]), lp[5])
  }))
  expect_equal(unname(sqrt(diag(vcov(e)))[4:5]), apply(sim, 2, sd), tolerance=0.05)

  # a group not in the fit, and random=FALSE: the population effective dose
  ep <- ed(f, c(20, 50))
  expect_equal(unname(coef(e)[c(3, 6)]), unname(coef(ep)))
  expect_equal(unname(diag(vcov(e))[c(3, 6)]), unname(diag(vcov(ep))))
  e0 <- ed(f, 50, newdata=nd, random=FALSE)
  expect_equal(unname(coef(e0)), rep(unname(coef(ed(f, 50, random=FALSE))), 3))
  expect_false(any(e0$conditional))
})

test_that("intervals of effective doses of fitted groups", {
  d <- simre()
  f <- fit4(y ~ dose, rform=b1 + b4 ~ 1 | grp, data=d, model="loglogistic", b4link="log")
  e <- ed(f, newdata=data.frame(grp=c("2", "7")))
  expect_message(ci <- confint(e), "not available for effective doses of fitted groups")
  expect_equal(ci, confint(e, method="wald"))
  expect_error(profile(e), "fitted groups")
  cb <- suppressMessages(confint(e, method="boot", nsim=20, seed=2, progress=FALSE))
  expect_true(all(cb[, 1] < coef(e) & coef(e) < cb[, 2]))
  # with a new group, the variance of the random effects is included: Wald intervals
  en <- ed(f, newdata=data.frame(grp="new"))
  expect_message(confint(en, method="boot"), "only available for population")
})

test_that("predict uses the population curve for rows with a missing group", {
  d <- simre(ng=6)
  f <- fit4(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="loglogistic")
  p <- predict(f, newdata=data.frame(dose=2, grp=c("2", NA, "3")))
  expect_equal(p[2], predict(f, newdata=data.frame(dose=2, grp="new")))
  expect_equal(p[c(1, 3)], predict(f, newdata=data.frame(dose=2, grp=c("2", "3"))))
})

test_that("absolute effective doses: doses at which the curve reaches a response", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  e <- ed(f, c(2, 4, 6), type="absolute")
  expect_equal(names(coef(e)), c("y=2:lf1", "y=4:lf2", "y=6:lf3"))
  for (i in 1:3){
    y <- c(2, 4, 6)[i]
    expect_equal(coef(e)[[i]], uniroot(function(x) predict(f, data.frame(dose=x)) - y, c(0.01, 100), tol=1e-10)$root, tolerance=1e-6)
  }
  # delta method with the dependence on the asymptotes
  b <- bcoef(f)[1:4]
  edabs <- function(b, y) b[4] * (((y - b[1])/(b[2] - b[1]))^(-1) - 1)^(-1/b[3])
  g <- numDeriv::grad(function(b) edabs(b, 4), b)
  expect_equal(sqrt(vcov(e)[2, 2]), sqrt(drop(g %*% vcov(f)[1:4, 1:4] %*% g)), tolerance=1e-6)
  # same dose as the relative effective dose at the corresponding response
  y50 <- (b[[1]] + b[[2]])/2
  expect_equal(unname(coef(ed(f, y50, type="absolute"))), unname(coef(ed(f, 50))))
  # responses outside the curve
  expect_warning(ea <- ed(f, c(4, 20), type="absolute"), "do not reach the response for y=20:lf2")
  expect_true(is.na(coef(ea)[2]) && !is.na(coef(ea)[1]))
  expect_error(profile(ea, parm=2), "NA")

  # profile and bootstrap intervals
  ci <- confint(profile(e, parm=2, progress=FALSE))
  expect_true(ci[1] < coef(e)[2] && coef(e)[2] < ci[2])
  expect_equal(ci, confint(e, parm=2, method="wald"), tolerance=0.05)
  cb <- suppressMessages(confint(e, method="boot", nsim=20, seed=1, progress=FALSE))
  expect_true(all(cb[, 1] < coef(e) & coef(e) < cb[, 2]))
})

test_that("absolute effective doses agree with drc", {
  skip_if_not_installed("drc")
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4())
  expect_equal(unname(coef(ed(f, c(2, 4, 6), type="absolute"))),
               unname(drc::ED(f2, c(2, 4, 6), type="absolute", display=FALSE)[, 1]), tolerance=1e-4)
})

test_that("absolute effective doses of fitted groups depend on their random asymptotes", {
  d <- simre()
  f <- fit4(y ~ dose, rform=b1 + b4 ~ 1 | grp, data=d, model="loglogistic", b4link="log")
  e <- ed(f, 4, newdata=data.frame(grp=c("2", "new")), type="absolute")
  root <- uniroot(function(x) predict(f, data.frame(dose=x, grp="2")) - 4, c(0.01, 100), tol=1e-10)$root
  expect_equal(coef(e)[[1]], root, tolerance=1e-6)
  # the random effects of b1 change the absolute effective dose of a new group, not the relative one
  er <- ed(f, newdata=data.frame(grp="new"))
  er0 <- ed(f, newdata=data.frame(grp="new"), random=FALSE)
  ea0 <- ed(f, 4, newdata=data.frame(grp="new"), random=FALSE, type="absolute")
  expect_gt(vcov(e)[2, 2] - vcov(ea0)[1, 1], 0)
  expect_equal(c(vcov(er) - vcov(er0)), unname(coef(er)^2 * recov(f)[[1]]$cov[2, 2]), tolerance=1e-8)
})
