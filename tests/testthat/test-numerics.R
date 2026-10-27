# Internal numerical checks: automatic derivatives against finite differences, effective doses
# against a numerical inverse of the mean function, and invariances of the fitted models.

skip_if_not_installed("numDeriv")

# a point away from the optimum, where the gradient is not zero
offopt <- function(fit) fit$opt$par * 1.05 + 0.01

test_that("gradients and hessians agree with finite differences for all models", {
  for (m in models){
    fit <- fit4(y ~ dose, data=simgauss(m), model=m)
    p <- offopt(fit)
    expect_equal(as.vector(fit$obj$gr(p)), numDeriv::grad(fit$obj$fn, p), tolerance=1e-6, label=m)
    expect_equal(unname(fit$obj$he(p)), numDeriv::jacobian(fit$obj$gr, p), tolerance=1e-6, label=m)
  }
})

test_that("gradients agree with finite differences for the families and links", {
  set.seed(5)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=4), n=20)
  mu <- drcmean("loglogistic", d$dose, c(0.1, 0.9, 1.5, 2.5, 0))
  d$yb <- rbinom(nrow(d), d$n, mu)
  d$yp <- rpois(nrow(d), 40*mu)
  d$yr <- rbeta(nrow(d), 30*mu, 30*(1 - mu))
  fits <- list(
    binomial=fit4(cbind(yb, n - yb) ~ dose, data=d, model="loglogistic", family="binomial"),
    binomiallogit=fit4(cbind(yb, n - yb) ~ dose, data=d, model="loglogistic", family="binomial", link="logit"),
    poissonlog=fit4(yp ~ dose, data=d, model="loglogistic", family="poisson", link="log"),
    beta=fit4(yr ~ dose, data=d, model="loglogistic", family="beta", link="logit"),
    asymmetricll=drmTMB(y ~ dose, data=simgauss("loglogistic"), model="loglogistic"),
    asymmetricl=drmTMB(y ~ dose, data=simgauss("logistic"), model="logistic")
  )
  for (nm in names(fits)){
    fit <- fits[[nm]]
    p <- offopt(fit)
    expect_equal(as.vector(fit$obj$gr(p)), numDeriv::grad(fit$obj$fn, p), tolerance=1e-6, label=nm)
  }
})

test_that("gradients of the Laplace approximation agree with finite differences", {
  d <- simgroups(ng=8, seed=6)
  u <- cbind(rnorm(8, 0, 0.5), rnorm(8, 0, 0.3))
  d$y <- 1 + u[d$grp, 1] + (7 - u[d$grp, 1])/(1 + exp(1.5*(log(d$dose) - log(3*exp(u[d$grp, 2]))))) + rnorm(nrow(d), 0, 0.3)
  fits <- list(
    single=fit4(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="loglogistic"),
    correlated=fit4(y ~ dose, rform=b1 + b4 ~ 1 | grp, data=d, model="loglogistic"),
    independent=fit4(y ~ dose, rform=list(b1 ~ 1 | grp, b4 ~ 1 | grp), data=d, model="loglogistic")
  )
  for (nm in names(fits)){
    fit <- fits[[nm]]
    p <- offopt(fit)
    expect_equal(as.vector(fit$obj$gr(p)), numDeriv::grad(fit$obj$fn, p), tolerance=1e-5, label=nm)
  }
})

# effective dose by numerical inversion of the mean function: the dose at which the response
# has changed by p percent from the control asymptote towards the asymptote at high doses
numed <- function(m, b, p){
  mod <- match(m, models)
  control <- if (mod %in% 6:8 || b[3] > 0) b[1] else b[2]
  other <- if (control == b[1]) b[2] else b[1]
  target <- control + p/100 * (other - control)
  # search on the log-dose scale for all models with positive doses
  logdose <- mod %in% 2:8
  g <- function(z) drcmean(m, if (logdose) exp(z) else z, b) - target
  z <- uniroot(g, c(-5, 5), extendInt="yes", tol=1e-12)$root
  if (logdose) exp(z) else z
}

test_that("effective doses and their standard errors agree with a numerical inverse", {
  for (m in models){
    fit <- fit4(y ~ dose, data=simgauss(m), model=m)
    b <- numeric(5)
    b[!fit$fix] <- bcoef(fit)
    b[fit$fix] <- fit$start[fit$fix]
    e <- ed(fit, respLev=c(10, 50, 90))
    expect_equal(unname(e$estimates), sapply(c(10, 50, 90), function(p) numed(m, b, p)), tolerance=1e-6, label=m)

    # delta method with numerical derivatives of the inverse
    free <- which(!fit$fix)
    V <- vcov(fit)[seq_along(free), seq_along(free)]
    se <- sapply(c(10, 50, 90), function(p){
      g <- numDeriv::grad(function(bf){ bb <- b; bb[free] <- bf; numed(m, bb, p) }, b[free])
      sqrt(drop(t(g) %*% V %*% g))
    })
    expect_equal(unname(sqrt(diag(e$variance))), se, tolerance=1e-4, label=m)
  }
})

test_that("effective doses are defined relative to the control for increasing and decreasing curves", {
  d <- simgauss("loglogistic")
  dn <- d
  dn$y <- 9 - d$y
  fd <- fit4(y ~ dose, data=d, model="loglogistic")
  fi <- fit4(y ~ dose, data=dn, model="loglogistic")
  # mirrored response: same ED, and ED 10 below ED 90 for both
  expect_equal(unname(ed(fd, respLev=c(10, 90))$estimates), unname(ed(fi, respLev=c(10, 90))$estimates), tolerance=1e-5)
  expect_true(diff(ed(fd, respLev=c(10, 90))$estimates) > 0)
})

test_that("predict equals the mean function at the estimates", {
  for (m in models){
    d <- simgauss(m)
    fit <- fit4(y ~ dose, data=d, model=m)
    b <- numeric(5)
    b[!fit$fix] <- bcoef(fit)
    b[fit$fix] <- fit$start[fit$fix]
    expect_equal(predict(fit), drcmean(m, d$dose, b), tolerance=1e-10, label=m)
  }
})

test_that("weights are equivalent to replicated observations", {
  # for the likelihood; small-sample inference counts the observations, as lm
  d <- simgauss("loglogistic")
  w <- rep(1:2, length.out=nrow(d))
  ml <- list(smallsample=FALSE)
  f1 <- fit4(y ~ dose, data=d, model="loglogistic", weights=w, control=ml)
  f2 <- fit4(y ~ dose, data=d[rep(seq_len(nrow(d)), w), ], model="loglogistic", control=ml)
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-8)
  expect_equal(coef(f1), coef(f2), tolerance=1e-5)
  expect_equal(vcov(f1), vcov(f2), tolerance=1e-4)
})

test_that("the order of the observations does not matter", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  set.seed(8)
  dp <- d[sample(nrow(d)), ]
  f1 <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  f2 <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=dp, model="loglogistic")
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-8)
  expect_equal(coef(f1), coef(f2), tolerance=1e-5)
})

test_that("contrasts of the fixed effects do not change the fit", {
  d <- rbind(cbind(simgauss("loglogistic", seed=1), trt="a"),
             cbind(simgauss("loglogistic", seed=2), trt="b"),
             cbind(simgauss("loglogistic", seed=3), trt="c"))
  d$trt <- factor(d$trt)
  d$y[d$trt == "b"] <- d$y[d$trt == "b"] + 1
  f1 <- fit4(y ~ dose, fform=b2 + b4 ~ trt, data=d, model="loglogistic")
  d2 <- d
  contrasts(d2$trt) <- contr.sum(3)
  f2 <- fit4(y ~ dose, fform=b2 + b4 ~ trt, data=d2, model="loglogistic")
  f3 <- fit4(y ~ dose, fform=b2 + b4 ~ 0 + trt, data=d, model="loglogistic")
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-7)
  expect_equal(f1$opt$objective, f3$opt$objective, tolerance=1e-7)
  expect_equal(predict(f1), predict(f2), tolerance=1e-5)
  expect_equal(predict(f1), predict(f3), tolerance=1e-5)
})

test_that("rescaling the dose rescales b4 only", {
  d <- simgauss("loglogistic")
  d2 <- d
  d2$dose <- d$dose * 1000
  f1 <- fit4(y ~ dose, data=d, model="loglogistic")
  f2 <- fit4(y ~ dose, data=d2, model="loglogistic")
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-7)
  expect_equal(unname(bcoef(f2)), unname(bcoef(f1) * c(1, 1, 1, 1000)), tolerance=1e-4)
  expect_equal(unname(ed(f2)$estimates), unname(ed(f1)$estimates * 1000), tolerance=1e-4)
})

test_that("a non-finite objective at the estimates gives a warning", {
  # responses at dose zero are impossible with the lower asymptote fixed at 0
  set.seed(13)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=3), n=20)
  d$y <- rbinom(nrow(d), d$n, drcmean("loglogistic", d$dose, c(0.1, 0.9, 2, 2.5, 0)))
  expect_warning(drmTMB(cbind(y, n - y) ~ dose, data=d, model="loglogistic", family="binomial",
                        start=c(b1=0, b2=1, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE)),
                 "not finite")
})
