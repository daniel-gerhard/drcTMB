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
