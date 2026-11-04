# links of the parameters b1-b5

# quantal data with a background response b1 and b2 at 1: a loglogistic curve on the probability scale
simquantal <- function(seed=1){
  set.seed(seed)
  q <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=4), n=20)
  q$dead <- rbinom(nrow(q), q$n, drcTMB:::drcmeanR(2, q$dose, 0.1, 1, 1.5, 3, 1))
  q
}

test_that("the links of a fit are the defaults of the model, by parameter", {
  d <- simgauss("loglogistic")
  f <- drmTMB(y ~ dose, data=d, model="loglogistic")
  expect_equal(f$link, c(b1="identity", b2="identity", b3="identity", b4="identity", b5="log"))
  f <- drmTMB(y ~ dose, data=simgauss("weibull1"), model="weibull1", link="log")
  expect_equal(f$link, c(b1="log", b2="log", b3="identity", b4="identity", b5="identity"))
  # the constant model does not use b2
  f <- drmTMB(y ~ dose, data=d, model="constant", link="log")
  expect_equal(f$link[1:2], c(b1="log", b2="identity"))
  expect_equal(unname(exp(coef(f)[1])), mean(d$y), tolerance=1e-6)
})

test_that("the asymptotes can have different links", {
  q <- simquantal()
  # background on the logit scale, b2 fixed at 1 on the probability scale
  f <- drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="loglogistic",
              link=c(b1="logit"), start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  expect_equal(f$link[1:2], c(b1="logit", b2="identity"))
  cf <- coef(f)
  nd <- data.frame(dose=c(0, 1, 3, 10))
  expect_equal(predict(f, newdata=nd), drcTMB:::drcmeanR(2, nd$dose, plogis(cf[[1]]), 1, cf[[2]], cf[[3]], 1))
  # the same fit as with the identity link for both asymptotes
  f0 <- drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="loglogistic",
               start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f0)), tolerance=1e-6)
  expect_equal(plogis(cf[[1]]), coef(f0)[[1]], tolerance=1e-4)
  # effective doses with their standard errors do not depend on the parameterisation
  expect_equal(coef(ed(f, c(10, 50))), coef(ed(f0, c(10, 50))), tolerance=1e-4)
  expect_equal(vcov(ed(f, c(10, 50))), vcov(ed(f0, c(10, 50))), tolerance=1e-3)
  expect_output(print(f), "logit link for b1, identity link for b2")
})

test_that("link is checked", {
  d <- simgauss("loglogistic")
  expect_error(drmTMB(y ~ dose, data=d, model="loglogistic", link="probit"), "Unknown link")
  expect_error(drmTMB(y ~ dose, data=d, model="loglogistic", link=c("log", "log")), "single link")
  expect_error(drmTMB(y ~ dose, data=d, model="loglogistic", link=c(b6="log")), "parameter names b1-b5")
  expect_error(drmTMB(y ~ dose, data=d, model="loglogistic", link=c(b3="log")), "not available for b3")
  expect_error(drmTMB(y ~ dose, data=d, model="loglogistic", link=c(b5="identity")), "not available for b5")
  expect_error(drmTMB(y ~ dose, data=d, model="logistic", link=c(b4="log")), "positive b4")
})
