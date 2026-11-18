# Bayesian model-averaged BMD: weights from the marginal likelihoods, mixture quantiles, checks

data(ryegrass)

prb <- list(bmd=pert(0.1, 1, 5), b2=normal(8, 1), b3=normal(-2, 2), sigma=half_normal(1))
fitb <- function(m, prior=prb, bmdpar=list(bmr=0.05, def="relative"), ...){
  drmTMB(rootl ~ conc, data=ryegrass, model=m, start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE), bmdpar=bmdpar, prior=prior, ...)
}
LL <- fitb("loglogistic")
W1 <- fitb("weibull1")
W2 <- fitb("weibull2")

test_that("a single model gives the normal posterior of its log BMD", {
  b <- bmdbayes(LL)
  i <- which(names(coef(LL)) == "log_bmd")
  mu <- coef(LL)[[i]]
  s <- sqrt(vcov(LL)[i, i])
  expect_equal(b$models$weight, 1)
  expect_equal(unname(coef(b)), exp(mu), tolerance=1e-8)
  expect_equal(unname(confint(b)[1, ]), exp(mu + qnorm(c(0.05, 0.95))*s), tolerance=1e-8)
  expect_equal(unname(confint(b, level=0.5)[1, ]), exp(mu + qnorm(c(0.25, 0.75))*s), tolerance=1e-8)
  expect_equal(rownames(b$models), "LL")
})

test_that("posterior model probabilities and the mixture of the models", {
  b <- bmdbayes(LL, W1, weib2=W2)
  expect_equal(rownames(b$models), c("LL", "W1", "weib2"))
  ml <- c(marglik(LL), marglik(W1), marglik(W2))
  expect_equal(b$models$logML, ml)
  expect_equal(b$models$weight, exp(ml - max(ml))/sum(exp(ml - max(ml))))
  # quantiles of the mixture against a large sample of it
  set.seed(1)
  N <- 400000
  k <- sample(3, N, replace=TRUE, prob=b$models$weight)
  mu <- vapply(b$posterior, function(p) p$mean, numeric(1))
  s <- vapply(b$posterior, function(p) p$sd, numeric(1))
  z <- rnorm(N, mu[k], s[k])
  expect_equal(log(c(coef(b), confint(b))), unname(quantile(z, c(0.5, 0.05, 0.95))), tolerance=0.005, ignore_attr=TRUE)
  # the BMDs of the models are their posterior modes, with normal credible intervals
  bm <- bmdbayesmodels(b)
  expect_equal(bm$BMD, unname(exp(mu)), tolerance=1e-8)
  expect_equal(bm$BMDL, unname(exp(mu + qnorm(0.05)*s)), tolerance=1e-8)
  expect_output(print(b), "Bayesian model-averaged benchmark dose, relative deviation 0.05, prior pert\\(0.1, 1, 5\\)")
  expect_output(print(b), "model average")

  # prior probabilities: a model with prior probability 0 is left out
  b0 <- bmdbayes(LL, W1, W2, prior.weights=c(1, 0, 1))
  b2 <- bmdbayes(LL, W2)
  expect_equal(b0$models$weight[2], 0)
  expect_equal(coef(b0), coef(b2), tolerance=1e-8)
  expect_equal(confint(b0), confint(b2), tolerance=1e-8)
  bp <- bmdbayes(LL, W2, prior.weights=c(3, 1))
  r <- exp(b2$models$logML[1] - b2$models$logML[2])
  expect_equal(bp$models$weight[1], 3*r/(3*r + 1))
})

test_that("importance sampling of the marginal likelihoods and the posterior of the BMD", {
  set.seed(2)
  bi <- bmdbayes(LL, W2, method="is", n=4000)
  set.seed(2)
  m1 <- marglik(LL, "is", n=4000)
  expect_equal(bi$models$logML[1], as.numeric(m1))
  expect_equal(bi$models$ess[1], attr(m1, "ess"))
  # close to the Laplace approximation for these well-identified fits
  b <- bmdbayes(LL, W2)
  expect_equal(coef(bi), coef(b), tolerance=0.05)
  expect_equal(confint(bi), confint(b), tolerance=0.05)
  expect_output(print(bi), "Importance sampling")
  # weighted quantiles
  expect_equal(wquantile(c(3, 1, 2), c(1, 1, 2), c(0.2, 0.25, 0.5, 0.75, 1)), c(1, 1, 2, 2, 3))
})

test_that("random effects on b2", {
  set.seed(6)
  d <- expand.grid(dose=c(0, 0.5, 1, 2, 4, 8, 16), grp=factor(1:6))
  u <- rnorm(6, 0, 0.8)
  d$y <- 1 + (7 + u[d$grp] - 1)/(1 + (d$dose/3)^1.5) + rnorm(nrow(d), 0, 0.3)
  pr <- list(bmd=pert(0.05, 0.5, 3), b1=normal(0, 2), b2=normal(8, 3), b3=normal(-2, 2), sd=half_normal(1), sigma=half_normal(1))
  fr <- function(m) drmTMB(y ~ dose, rform=b2 ~ 1 | grp, data=d, model=m, start=c(b5=0), fix=c(b5=TRUE),
                           bmdpar=list(bmr=0.1, def="relative"), prior=pr)
  b <- bmdbayes(fr("loglogistic"), fr("lognormal"))
  expect_equal(sum(b$models$weight), 1)
  ci <- confint(b)
  expect_true(ci[1] < coef(b) && coef(b) < ci[2])
})

test_that("checks of the models, failed fits and notes", {
  f0 <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE), prior=prb[-1])
  expect_error(bmdbayes(LL, f0), "f0 is not fitted in the BMD parameterisation")
  expect_error(bmdbayes(LL, fitb("weibull1", bmdpar=list(bmr=0.1, def="relative"))), "same benchmark response")
  expect_error(bmdbayes(LL, fitb("weibull1", prior=c(list(bmd=pert(0.1, 2, 5)), prb[-1]))), "same prior for the BMD")
  expect_error(bmdbayes(LL, fitb("weibull1", prior=prb[-1])), "need a prior for the BMD")
  expect_error(bmdbayes(LL, fitb("weibull1", prior=prb[-4])), "model fitb\\(\"weibull1\", prior = prb\\[-4\\]\\): .*no prior for sigma")
  expect_error(bmdbayes(LL, drmTMB(rootl ~ conc, data=ryegrass[-1, ], model="weibull1", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
                                   bmdpar=list(bmr=0.05, def="relative"), prior=prb)), "same doses")
  expect_error(bmdbayes(LL, W1, prior.weights=c(1, -1)), "prior.weights")
  expect_error(bmdbayes(LL, W1, prior.weights=1), "prior.weights")
  expect_error(bmdbayes(LL, W1, method="bridge"))

  # a model whose fit did not converge gets the weight 0
  Wbad <- W1
  Wbad$opt$convergence <- 1
  expect_warning(bb <- bmdbayes(LL, Wbad, W2), "model Wbad gets the weight 0: the fit did not converge")
  expect_equal(bb$models$weight[2], 0)
  expect_true(is.na(bmdbayesmodels(bb)$BMD[2]))
  expect_equal(coef(bb), coef(bmdbayes(LL, W2)), tolerance=1e-8)
  expect_error(suppressWarnings(bmdbayes(Wbad)), "none of the models")

  # notes for wide credible intervals
  b <- bmdbayes(LL)
  expect_no_match(paste(capture.output(print(b)), collapse="\n"), "Note")
  b$posterior[[1]]$sd <- 3
  b$estimates <- c(BMD=exp(bmdbayesq(b, 0.5)))
  expect_output(print(b), "Note \\(EFSA 2022\\): BMD/BMDL above 20, BMDU/BMDL above 50, BMDL below a tenth of the lowest dose")
})

test_that("importance sampling agrees with the posterior of tmbstan for a curved posterior", {
  # weibull1: the log BMD is strongly correlated with the skewed slope b3, so that the normal
  # approximation misses the upper tail; adaptive importance sampling follows it
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("DRCTMB_SLOW_TESTS"), "true"), "set DRCTMB_SLOW_TESTS=true to run")
  skip_if_not_installed("tmbstan")
  st <- suppressWarnings(tmbstan::tmbstan(W1$obj, chains=4, iter=6000, warmup=1000, lower=W1$lower, upper=W1$upper,
                                          init=function() as.numeric(W1$opt$par), refresh=0, seed=1, silent=TRUE))
  X <- as.matrix(st)
  lb <- X[, which(parnames(W1) == "log_bmd")]
  set.seed(3)
  bi <- bmdbayes(W1, method="is", n=20000)
  # absolute differences on the log scale: the Monte Carlo error of the 5% and 95% quantiles is about 0.02
  expect_lt(max(abs(log(c(coef(bi), confint(bi))) - quantile(lb, c(0.5, 0.05, 0.95)))), 0.05)
})
