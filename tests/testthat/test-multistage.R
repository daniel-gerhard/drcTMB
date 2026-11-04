# Multistage model, F(x) = 1 - exp(-(b3 x + b4 x^2 + b5 x^3)), with b3-b5 bounded below by 0 as in
# BMDS. drc has no multistage model; the reference is the binomial likelihood maximised by nlminb
# with the same bounds, on rescaled coefficients, as the coefficients of the powers of the dose
# differ by orders of magnitude.

simms <- function(seed=5){
  set.seed(seed)
  q <- data.frame(dose=c(0, 1, 2, 4, 8, 16), n=50)
  q$dead <- rbinom(nrow(q), q$n, drcTMB:::drcmeanR(14, q$dose, 0.05, 1, 0.01, 0.004, 0))
  q
}

# multistage fit of degree deg to quantal data q, with b2 fixed at 1
fitms <- function(q, deg){
  start <- c(b2=1, b4=0, b5=0)[c(TRUE, deg < 2, deg < 3)]
  fix <- c(b2=TRUE, b4=deg < 2, b5=deg < 3)
  drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="multistage", start=start, fix=fix)
}

# reference fit: background and coefficients of degree deg
refms <- function(q, deg){
  nll <- function(p){
    beta <- c(p[-1], 0, 0, 0)[1:3]
    P <- p[1] + (1 - p[1])*(1 - exp(-(beta[1]*q$dose + beta[2]*q$dose^2 + beta[3]*q$dose^3)))
    -sum(dbinom(q$dead, q$n, pmin(P, 1 - 1e-12), log=TRUE))
  }
  sc <- c(1, max(q$dose)^(1:deg))
  o <- nlminb(c(0.05, rep(0.5, deg)), function(p) nll(p/sc), lower=c(1e-8, rep(0, deg)), upper=c(0.99, rep(Inf, deg)),
              control=list(rel.tol=1e-14, eval.max=1e4, iter.max=1e4))
  list(par=o$par/sc, nll=o$objective)
}

# free and fixed coefficients b1-b5 of a fit
bfull <- function(fit){
  b <- fit$start
  b[!fit$fix] <- coef(fit)[seq_len(sum(!fit$fix))]
  b
}

test_that("multistage fits agree with the maximum likelihood under the bounds", {
  q <- simms()
  for (deg in 1:3){
    f <- fitms(q, deg)
    r <- refms(q, deg)
    expect_equal(f$opt$objective, r$nll, tolerance=1e-8)
    expect_equal(unname(bcoef(f)), r$par, tolerance=1e-5)
  }
  # the cubic fit has b4 at its bound 0: listed by summary, with the projected gradient
  s <- summary(fitms(q, 3))
  expect_equal(names(s$atbound), "b4:(Intercept)")
  expect_lt(s$maxgrad, 1e-4)
  expect_output(print(s), "At a bound: b4")
})

test_that("the multistage model of degree 1 is the asymptoticregression model", {
  q <- simms()
  f1 <- fitms(q, 1)
  f2 <- drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="asymptoticregression", start=c(b2=1), fix=c(b2=TRUE))
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-8)
  expect_equal(coef(f1)[[2]], 1/coef(f2)[[2]], tolerance=1e-5)
  expect_equal(unname(coef(ed(f1, c(5, 10)))), unname(coef(ed(f2, c(5, 10)))), tolerance=1e-5)
})

test_that("effective doses of the multistage model have the delta method standard errors", {
  skip_if_not_installed("numDeriv")
  q <- simms()
  for (deg in 2:3){
    f <- fitms(q, deg)
    b <- bfull(f)
    free <- !f$fix
    e <- ed(f, c(5, 10))
    edR <- function(bf, p){
      bb <- b
      bb[free] <- bf
      drcTMB:::msroot(p/100, bb[3:5])
    }
    expect_equal(unname(coef(e)), sapply(c(5, 10), function(p) edR(b[free], p)))
    g <- t(sapply(c(5, 10), function(p) numDeriv::grad(edR, b[free], p=p)))
    expect_equal(unname(vcov(e)), unname(g %*% vcov(f)[seq_len(sum(free)), seq_len(sum(free))] %*% t(g)), tolerance=1e-6)
    # the extra risk at the effective dose
    P <- function(x) drcTMB:::drcmeanR(14, x, b[1], b[2], b[3], b[4], b[5])
    expect_equal((P(coef(e)[[2]]) - b[[1]])/(1 - b[[1]]), 0.1)
    # each coefficient of the polynomial from the effective dose, as in the profiles
    for (ps in which(b[3:5] > 0) + 2) expect_equal(drcTMB:::msbeta(0.1, coef(e)[[2]], b, ps), b[[ps]])
  }
})

test_that("profile intervals of multistage effective doses respect the bounds", {
  q <- simms()
  for (deg in 2:3){
    f <- fitms(q, deg)
    e <- ed(f, 10)
    pr <- profile(e, level=0.95, progress=FALSE)
    ci <- confint(pr, level=0.90)
    expect_true(all(is.finite(ci)))
    expect_true(ci[1, 1] < coef(e) && coef(e) < ci[1, 2])
    # the profile deviance at the limits equals the cutoff, with the coefficients within their bounds
    expect_true(all(pr[[1]]$deviance >= -1e-6))
  }
})

test_that("multistage bounds and options", {
  q <- simms()
  f <- fitms(q, 3)
  expect_equal(f$lower[1:3], c(-Inf, 0, 0))
  expect_equal(f$lower[4], 0)
  # bounds given by the user are used instead
  f2 <- drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="multistage", start=c(b2=1, b5=0),
               fix=c(b2=TRUE, b5=TRUE), lower=c(0, 1, 0, -1, 0))
  expect_equal(f2$lower[1:3], c(0, 0, -1))
  expect_error(drmTMB(cbind(dead, n - dead) ~ dose, data=q, family="binomial", model="multistage", link=c(b4="log")),
               "not available for the multistage")
  # mean of the template and of R
  b <- bfull(f)
  nd <- data.frame(dose=c(0, 3, 10))
  expect_equal(predict(f, newdata=nd), drcTMB:::drcmeanR(14, nd$dose, b[1], b[2], b[3], b[4], b[5]))
})

test_that("standard errors of coefficients on a small scale are finite", {
  # selenium, type 1: b4 of the quadratic model is about 2.5e-6, smaller than the absolute steps of
  # the finite-difference Hessian of sdreport
  skip_if_not_installed("drcData")
  data(selenium, package="drcData")
  se1 <- subset(selenium, type == 1)
  f <- drmTMB(cbind(dead, total - dead) ~ conc, data=se1, family="binomial", model="multistage",
              start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  expect_true(all(is.finite(sqrt(diag(vcov(f))))))
  expect_equal(unname(vcov(f)), unname(solve(f$obj$he(f$opt$par))), tolerance=1e-6)
  e <- ed(f, 10)
  expect_true(all(is.finite(confint(e, method="wald"))))
  # b5 of the cubic model is about 5e-9 but not at its bound: only b4 is
  f3 <- drmTMB(cbind(dead, total - dead) ~ conc, data=se1, family="binomial", model="multistage",
               start=c(b2=1), fix=c(b2=TRUE))
  expect_equal(names(summary(f3)$atbound), "b4:(Intercept)")
})
