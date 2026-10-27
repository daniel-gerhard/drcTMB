# Random slopes: a random coefficient for each parameter and each column of the term, e.g.
# b2 ~ temp | grp gives correlated random intercepts and slopes in temp of b2.

# loglogistic data with temp varying within groups; random intercept and slope of b2 or of log b4
simslopes <- function(par=c("b2", "b4"), ng=10, nrep=1, seed=11){
  par <- match.arg(par)
  set.seed(seed)
  d <- expand.grid(rep=seq_len(nrep), dose=c(0, 0.5, 1, 2, 4, 8, 16, 32), temp=c(-1, 0, 1), grp=factor(seq_len(ng)))
  d$rep <- NULL
  u0 <- rnorm(ng, 0, 0.4)
  u1 <- 0.5*u0 + rnorm(ng, 0, 0.2)
  ug <- u0[d$grp] + u1[d$grp]*d$temp
  b2 <- 8 + d$temp + if (par == "b2") ug else 0
  b4 <- 3*exp(0.4*d$temp + if (par == "b4") 0.7*ug else 0)
  d$y <- drcTMB:::drcmeanR(2, d$dose, 2, b2, -1.5, b4, 0) + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("random slopes of an asymptote maximise the same likelihood as nlme", {
  skip_if_not_installed("nlme")
  d <- simslopes("b2")
  # nlme cannot evaluate log(0); a tiny positive dose has the same mean as dose zero
  d$dose <- pmax(d$dose, 1e-8)
  f <- fit4(y ~ dose, fform=b2 ~ temp, rform=b2 ~ temp | grp, data=d, model="loglogistic")
  n1 <- nlme::nlme(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))),
                   fixed=list(b1 ~ 1, b2 ~ temp, b3 ~ 1, b4 ~ 1), random=nlme::pdSymm(b2 ~ temp), groups=~grp,
                   data=d, start=unname(bcoef(f)), method="ML")
  # linear in the random effects: the Laplace approximation and the linearisation of nlme are exact
  expect_lte(f$opt$objective, -as.numeric(logLik(n1)) + 1e-6)
  expect_equal(f$opt$objective, -as.numeric(logLik(n1)), tolerance=1e-3)
  expect_equal(unname(bcoef(f)), unname(nlme::fixef(n1)), tolerance=1e-2)
  vc <- nlme::VarCorr(n1)
  expect_equal(recov(f)[[1]]$sd, as.numeric(vc[1:2, "StdDev"]), tolerance=1e-2, ignore_attr=TRUE)
  expect_equal(recov(f)[[1]]$corr[2, 1], as.numeric(vc[2, 3]), tolerance=1e-2)
  expect_equal(names(coef(f))[6:8], c("log_sd:b2|grp", "log_sd:b2:temp|grp", "corr_par:b2:temp,b2|grp"))
  out <- capture.output(print(f))
  expect_true(any(grepl("^ +b2:temp ", out)))
})

test_that("predict with random slopes", {
  d <- simslopes("b2")
  f <- fit4(y ~ dose, fform=b2 ~ temp, rform=b2 ~ temp | grp, data=d, model="loglogistic")
  cf <- coef(f)
  u <- matrix(f$sdrl$par.random, ncol=2)
  # fitted group 3 at temp 1, by hand
  p <- predict(f, data.frame(dose=2, temp=1, grp="3"))
  expect_equal(p, drcTMB:::drcmeanR(2, 2, cf[[1]], cf[[2]] + cf[[3]] + u[3, 1] + u[3, 2], cf[[4]], cf[[5]], 0))
  # a new group: the variance of the random intercept and slope at temp 1 in the prediction interval
  nd <- data.frame(dose=2, temp=1, grp="new")
  pc <- predict(f, nd, interval="confidence")
  pp <- predict(f, nd, interval="prediction")
  Fx <- (pc$fit - cf[[1]])/(cf[[2]] + cf[[3]] - cf[[1]])
  z <- c(1, 1)
  vnew <- Fx^2 * drop(z %*% recov(f)[[1]]$cov %*% z)
  expect_equal(((pp$upr - pp$fit)/qnorm(0.975))^2, pc$se.fit^2 + residpar(f)^2 + vnew, tolerance=1e-4)
  # simulation keeps working
  expect_equal(dim(simulate(f, 2, seed=1)), c(nrow(d), 2))
})

test_that("effective doses with random slopes of log b4", {
  d <- simslopes("b4", ng=12, nrep=2, seed=12)
  f <- fit4(y ~ dose, fform=b4 ~ temp, rform=b4 ~ temp | grp, data=d, model="loglogistic", b4link="log")
  nd <- data.frame(temp=c(1, 1, -1), grp=c("4", "new", "new"))
  e <- ed(f, 50, newdata=nd)
  expect_equal(names(coef(e)), c("50:temp=1,grp=4", "50:temp=1,grp=new", "50:temp=-1,grp=new"))
  # fitted group: the dose halfway between the asymptotes of its predicted curve
  pr <- function(dose) predict(f, data.frame(dose=dose, temp=1, grp="4"))
  ed50 <- uniroot(function(x) pr(x) - (pr(1e-8) + pr(1e8))/2, c(0.01, 100), tol=1e-10)$root
  expect_equal(coef(e)[[1]], ed50, tolerance=1e-6)
  # new group: variance of the random intercept and slope at temp, with dED/du = ED for log b4
  e0 <- ed(f, 50, newdata=nd, random=FALSE)
  S <- recov(f)[[1]]$cov
  for (i in 2:3){
    z <- c(1, nd$temp[i])
    expect_equal(diag(vcov(e))[[i]] - diag(vcov(e0))[[i]], coef(e)[[i]]^2 * drop(z %*% S %*% z), tolerance=1e-8)
  }
})

test_that("random slopes need their covariates", {
  d <- simslopes("b4", ng=12, nrep=2, seed=12)
  # temp only in the random effects
  f <- fit4(y ~ dose, rform=b4 ~ temp | grp, data=d, model="loglogistic", b4link="log")
  expect_error(ed(f), "covariates of the random slopes")
  expect_length(coef(ed(f, random=FALSE)), 1)
  expect_error(ed(f, newdata=data.frame(grp="4")), "covariates of the random slopes: temp")
  expect_error(ed(f, newdata=data.frame(temp=NA_real_, grp="4")), "missing values")
  expect_error(predict(f, data.frame(dose=2, grp="4")), "covariates of the random slopes: temp")
  expect_error(predict(f, data.frame(dose=2), random=FALSE, interval="prediction"), "covariates of the random slopes")
  expect_length(predict(f, data.frame(dose=2), random=FALSE), 1)
  expect_true(is.na(predict(f, data.frame(dose=2, temp=NA, grp="4"))))
})

test_that("uncorrelated random intercepts and slopes with ||", {
  d <- simslopes("b4", ng=12, nrep=2, seed=12)
  f <- fit4(y ~ dose, fform=b4 ~ temp, rform=b4 ~ temp || grp, data=d, model="loglogistic", b4link="log")
  expect_length(f$reterms, 2)
  expect_equal(names(coef(f))[6:7], c("log_sd:b4|grp", "log_sd:b4:temp|grp"))
  expect_equal(summary(f)$ngroups, c(grp=12))
})

test_that("random slopes of a factor use the levels and contrasts of the fit for any newdata", {
  set.seed(5)
  ng <- 10
  d <- expand.grid(rep=1:2, dose=c(0, 0.5, 1, 2, 4, 8, 16, 32), time=factor(c("24h", "48h", "72h")), grp=factor(seq_len(ng)))
  d$rep <- NULL
  U <- matrix(rnorm(ng*3, 0, 0.25), ng)
  b4 <- 3*exp(c(0, -0.3, -0.6)[d$time] + U[cbind(as.integer(d$grp), as.integer(d$time))])
  d$y <- drcTMB:::drcmeanR(2, d$dose, 1, 8, -1.5, b4, 0) + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, fform=b4 ~ time, rform=b4 ~ time | grp, data=d, model="loglogistic", b4link="log")
  expect_equal(names(coef(f))[7:9], c("log_sd:b4|grp", "log_sd:b4:time48h|grp", "log_sd:b4:time72h|grp"))

  lev <- c("24h", "48h", "72h")
  full <- predict(f, data.frame(dose=2, time=lev, grp="3"), interval="prediction")
  # a single level, a subset, and the levels in another order of the factor
  expect_equal(predict(f, data.frame(dose=2, time="72h", grp="3"), interval="prediction"), full[3, ], ignore_attr=TRUE)
  expect_equal(predict(f, data.frame(dose=2, time=c("48h", "72h"), grp="3"), interval="prediction"), full[2:3, ], ignore_attr=TRUE)
  nd <- data.frame(dose=2, time=factor(lev, levels=rev(lev)), grp="3")
  expect_equal(predict(f, nd, interval="prediction"), full, ignore_attr=TRUE)

  # effective doses of a fitted and a new group
  ef <- ed(f, newdata=expand.grid(time=lev, grp=c("3", "new")))
  e1 <- ed(f, newdata=data.frame(time="72h", grp=c("3", "new")))
  expect_equal(unname(coef(e1)), unname(coef(ef)[c(3, 6)]))
  expect_equal(unname(vcov(e1)), unname(vcov(ef)[c(3, 6), c(3, 6)]))
})
