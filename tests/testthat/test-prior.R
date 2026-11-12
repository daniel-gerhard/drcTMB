# priors: posterior modes, with the log Jacobian of the transformations and truncated priors

library(drcData)
data(ryegrass)

fitry <- function(...) drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0),
                              fix=c(b1=TRUE, b5=TRUE), ...)

# the estimates of fit with priors are a stationary point of the negative log-likelihood of fit0 (without
# priors) minus the log prior lp of the optimised parameters, computed in R: its gradient by numerical
# derivatives, relative to the gradient of the log-likelihood alone, and the objective of fit
expect_posterior_mode <- function(fit, fit0, lp, tol=1e-4){
  p <- fit$opt$par
  nlpost <- function(q) as.numeric(fit0$obj$fn(q) - lp(q))
  expect_equal(fit$opt$objective, nlpost(p), tolerance=1e-8)
  g <- numDeriv::grad(nlpost, p)
  gl <- numDeriv::grad(fit0$obj$fn, p)
  expect_lt(max(abs(g)/pmax(abs(gl), 1)), tol)
}

test_that("without priors the fit is unchanged", {
  f0 <- fitry()
  f1 <- fitry(prior=list())
  expect_identical(coef(f0), coef(f1))
  expect_null(f1$prior)
  expect_null(f1$nlprior)
})

test_that("posterior modes with priors on the natural scale, and the log-likelihood without the priors", {
  f0 <- fitry()
  pr <- list(b2=normal(8, 1), b4=pert(0.5, 3, 10), sigma=half_normal(1))
  f <- fitry(prior=pr)
  # the log prior, with the log Jacobian of log_sigma; the half normal is not truncated
  lp <- function(p) dnorm(p[1], 8, 1, log=TRUE) + lpdfprior(pert(0.5, 3, 10), p[3]) +
    log(2) + dnorm(exp(p[4]), 0, 1, log=TRUE) + p[4]
  expect_equal(f$nlprior, -unname(lp(f$opt$par)), tolerance=1e-10)
  expect_equal(as.numeric(logLik(f)), -f0$obj$fn(f$opt$par), tolerance=1e-10)
  expect_posterior_mode(f, f0, lp)
  # the pert prior bounds b4; small-sample inference is not used
  expect_equal(unname(f$lower[3]), 0.5)
  expect_equal(unname(f$upper[3]), 10)
  expect_null(f$rdf)
  expect_equal(f$prior$parameter, c("b2", "b4", "sigma"))
  expect_output(print(summary(f)), "Priors:")
  expect_output(print(f), "posterior modes")
  # a very informative prior holds the parameter at its mean
  fs <- fitry(prior=list(b2=normal(7, 1e-4)))
  expect_equal(coef(fs)[["b2:(Intercept)"]], 7, tolerance=1e-6)
})

test_that("priors after the log and logit links include the log Jacobian", {
  # log link for b4: lognormal prior for b4, the log Jacobian is log(b4)
  f0 <- fitry(link=c(b4="log"))
  f <- fitry(link=c(b4="log"), prior=list(b4=lognormal(log(2), 0.3)))
  lp <- function(p) dlnorm(exp(p[3]), log(2), 0.3, log=TRUE) + p[3]
  expect_posterior_mode(f, f0, lp)
  # logit link of the background b1 of a binomial model: beta-like pert prior on the probability
  data(daphnids)
  d <- daphnids[daphnids$time == "24h", ]
  g0 <- drmTMB(cbind(no, total - no) ~ dose, data=d, family="binomial", model="loglogistic", link=c(b1="logit"),
               start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  g <- drmTMB(cbind(no, total - no) ~ dose, data=d, family="binomial", model="loglogistic", link=c(b1="logit"),
              start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE), prior=list(b1=pert(0, 0.05, 0.3)))
  lpg <- function(p){
    z <- plogis(p[1])
    if (z >= 0.3) return(-Inf)
    lpdfprior(pert(0, 0.05, 0.3), z) + log(z) + log(1 - z)
  }
  expect_posterior_mode(g, g0, lpg)
  expect_equal(unname(g$upper[1]), qlogis(0.3))
})

test_that("priors on coefficients on the scale of the link, and truncated priors", {
  data(S.alba)
  f0 <- drmTMB(DryMatter ~ Dose, fform=b4 ~ Herbicide, data=S.alba, model="loglogistic",
               start=c(b5=0), fix=c(b5=TRUE), link=c(b4="log"))
  cn <- f0$coefnames[!f0$fix]
  pr <- list(normal(0, 0.5))
  names(pr) <- cn[5]
  f <- drmTMB(DryMatter ~ Dose, fform=b4 ~ Herbicide, data=S.alba, model="loglogistic",
              start=c(b5=0), fix=c(b5=TRUE), link=c(b4="log"), prior=pr)
  lp <- function(p) dnorm(p[5], 0, 0.5, log=TRUE)
  expect_posterior_mode(f, f0, lp)
  expect_equal(f$prior$scale, "link")
  expect_error(drmTMB(DryMatter ~ Dose, fform=b4 ~ Herbicide, data=S.alba, model="loglogistic",
                      start=c(b5=0), fix=c(b5=TRUE), prior=list(b4=normal(0, 1))), "single coefficient")

  # the power b3 of the exponential model is bounded below by 0: a normal prior is truncated at 0 and
  # renormalised, which does not change the estimates
  set.seed(1)
  d <- data.frame(dose=rep(c(0, 1, 2, 5, 10, 20), each=5))
  d$y <- 10*exp(0.02*d$dose^1.3) + rnorm(nrow(d), 0, 0.5)
  e1 <- drmTMB(y ~ dose, data=d, model="exponential", prior=list(b3=normal(1, 1)))
  expect_equal(e1$nlprior, -(dnorm(coef(e1)[["b3:(Intercept)"]], 1, 1, log=TRUE) - log(pnorm(1))), tolerance=1e-10)
})

test_that("priors on the standard deviations of random effects", {
  set.seed(4)
  d <- expand.grid(dose=c(0, 0.5, 1, 2, 4, 8), rep=1:2, grp=factor(1:5))
  u <- rnorm(5, 0, 0.5)
  d$y <- 2 + (8 + u[d$grp] - 2)/(1 + exp(2*(log(d$dose) - log(2)))) + rnorm(nrow(d), 0, 0.3)
  f0 <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  pr <- list(sd=half_student_t(3, 1), sigma=lognormal(log(0.3), 0.5))
  f <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), prior=pr)
  expect_equal(f$prior$parameter, c("sd:b2|grp", "sigma"))
  lp <- function(p) log(2) + dt(exp(p[5]), 3, log=TRUE) + p[5] + dlnorm(exp(p[6]), log(0.3), 0.5, log=TRUE) + p[6]
  expect_equal(f$nlprior, -unname(lp(f$opt$par)), tolerance=1e-8)
  expect_posterior_mode(f, f0, lp)
  # a specific name takes precedence over sd; a pert prior bounds the standard deviation
  f2 <- drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE),
               prior=list(sd=half_normal(1), "sd:b2|grp"=pert(0.1, 0.5, 2)))
  expect_equal(f2$prior$prior, "pert(0.1, 0.5, 2)")
  expect_equal(unname(f2$lower[5]), log(0.1))
  expect_true(all(is.finite(sqrt(diag(vcov(f2))))))
  # methods with priors
  expect_message(suppressWarnings(anova(f0, f)), "approximate")
  expect_true(all(is.finite(coef(ed(f, 50)))))
})

test_that("invalid priors", {
  expect_error(fitry(prior=normal(0, 1)), "named")
  expect_error(fitry(prior=list(b1=normal(0, 1))), "fixed")
  expect_error(fitry(prior=list(b9=normal(0, 1))), "unknown name")
  expect_error(fitry(prior=list(sd=half_normal(1))), "random effects")
  expect_error(fitry(prior=list(b2="normal")), "priors")
  expect_error(fitry(prior=list(phi=half_normal(1))), "unknown name")
  expect_error(fitry(link=c(b4="log"), prior=list(b4=pert(-2, -1, 0))), "no support")
  expect_error(pert(1, 0, 2), "min < mode < max")
  expect_error(normal(0, -1), "positive")
})
