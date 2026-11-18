# EFSA Bayesian BMD defaults: negate and exponential priors, the quadratic exponential BMD
# parameterisation, efsapriors, bmdefsa (against EFSA EN-9902) and fitpert

data(ryegrass)
data(TCDD)
data(liver.tumor)
fq <- cbind(incidence, total - incidence) ~ conc

test_that("priors of the negative and of the exponential of a parameter", {
  f0 <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))
  f <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
              prior=list(b3=negate(lognormal(1, 1))))
  nlpost <- function(q) unname(f0$obj$fn(q) - dlnorm(-q[2], 1, 1, log=TRUE))
  expect_equal(f$opt$objective, nlpost(f$opt$par), tolerance=1e-8)
  expect_lt(max(abs(numDeriv::grad(nlpost, f$opt$par))), 1e-4)
  expect_equal(unname(f$upper[2]), 0)
  expect_equal(f$prior$prior, "negate(lognormal(1, 1))")
  expect_error(negate(1), "needs a prior")
  expect_error(negate(negate(normal())), "parameter itself")

  # a PERT prior of the median exp(b2) with the log curve link: density of exp(b2) with the Jacobian
  fl0 <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", family="lognormal", curvelink="log", start=c(b5=0), fix=c(b5=TRUE))
  fl <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", family="lognormal", curvelink="log", start=c(b5=0), fix=c(b5=TRUE),
               prior=list(b2=onexp(pert(4, 7.5, 15))))
  a <- 1 + 4*3.5/11
  b <- 1 + 4*7.5/11
  nl2 <- function(q) unname(fl0$obj$fn(q) - (dbeta((exp(q[2]) - 4)/11, a, b, log=TRUE) - log(11) + q[2]))
  expect_equal(fl$opt$objective, nl2(fl$opt$par), tolerance=1e-8)
  expect_lt(max(abs(numDeriv::grad(nl2, fl$opt$par))), 1e-4)
  expect_equal(unname(c(fl$lower[2], fl$upper[2])), log(c(4, 15)))
  expect_equal(lpdfprior(onexp(pert(4, 7.5, 15)), log(6)), lpdfprior(pert(4, 7.5, 15), 6) + log(6))
})

test_that("priors of the background probability pnorm(b1) or plogis(b1) of the probit and logit models", {
  pb <- pert(0.001, 0.03, 0.1)
  for (cl in c("probit", "logit")){
    ifun <- if (cl == "probit") pnorm else plogis
    p <- oninv(pb, if (cl == "probit") "pnorm" else "plogis")
    f0 <- drmTMB(fq, data=TCDD, model="power", family="binomial", curvelink=cl, start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
    f <- drmTMB(fq, data=TCDD, model="power", family="binomial", curvelink=cl, start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE),
                prior=list(b1=p))
    # density of the background probability with the Jacobian of the inverse link
    nl <- function(q) unname(f0$obj$fn(q) - (lpdfprior(pb, ifun(q[1])) + log(numDeriv::grad(ifun, q[1]))))
    expect_equal(f$opt$objective, nl(f$opt$par), tolerance=1e-8)
    expect_lt(max(abs(numDeriv::grad(nl, f$opt$par))), 1e-4)
    expect_equal(unname(ifun(c(f$lower[1], f$upper[1]))), c(0.001, 0.1))
  }
  expect_error(oninv(normal(0, 1), "pnorm"), "support within")
  expect_error(oninv(pert(0, 0.5, 2), "plogis"), "support within")
})

test_that("BMD parameterisation of the quantal probit and logit models", {
  for (cl in c("probit", "logit")){
    f0 <- drmTMB(fq, data=TCDD, model="power", family="binomial", curvelink=cl, start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
    f1 <- drmTMB(fq, data=TCDD, model="power", family="binomial", curvelink=cl, start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE),
                 bmdpar=list(bmr=0.1, def="excess"))
    expect_equal(c(logLik(f1)), c(logLik(f0)), tolerance=1e-6)
    expect_equal(unname(exp(coef(f1)["log_bmd"])), unname(coef(bmd(f0, 0.1, "excess"))), tolerance=1e-5)
    # closed form: Phi^-1(p0 + q (1 - p0)) = b1 + b4 BMD^b3
    cf <- coef(f0)
    ifun <- if (cl == "probit") pnorm else plogis
    qfun <- if (cl == "probit") qnorm else qlogis
    p0 <- ifun(cf[[1]])
    expect_equal(unname(exp(coef(f1)["log_bmd"])), ((qfun(p0 + 0.1*(1 - p0)) - cf[[1]])/cf[[3]])^(1/cf[[2]]), tolerance=1e-5)
  }
  e <- efsapriors(fq, TCDD, "probit", "binomial")
  expect_equal(c(e$name, e$model, e$curvelink), c("probit", "power", "probit"))
  expect_match(e$prior$b1$name, "of pnorm$")
  expect_equal(e$prior$b3$name, "lognormal(1, 1)")
  expect_equal(efsapriors(fq, TCDD, "logit", "binomial")$curvelink, "logit")
  expect_error(efsapriors(rootl ~ conc, ryegrass, "probit"), "quantal responses")
})

test_that("the quadratic exponential BMD parameterisation of the multistage model of degree 2", {
  set.seed(2)
  d <- data.frame(dose=rep(c(0, 5, 10, 20, 40), each=4), n=25)
  d$y <- rbinom(nrow(d), d$n, 0.05 + 0.95*(1 - exp(-(0.01*d$dose + 0.0008*d$dose^2))))
  f0 <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="multistage", family="binomial", start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  fq <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="multistage", family="binomial", start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE),
               bmdpar=list(bmr=0.1, def="excess"))
  expect_equal(names(coef(fq)), c("b1:(Intercept)", "log_d", "log_bmd"))
  expect_equal(as.numeric(logLik(fq)), as.numeric(logLik(f0)), tolerance=1e-8)
  expect_equal(exp(coef(fq)[["log_bmd"]]), unname(coef(bmd(f0, 0.1, "excess"))), tolerance=1e-6)
  # both coefficients of the polynomial from the BMD and log d
  r <- fq$obj$report()
  expect_equal(c(r$f3bmd, r$f4bmd), unname(coef(f0)[2:3]), tolerance=1e-5)
  expect_output(print(fq), "log_d is the technical parameter of the quadratic exponential model")
})

test_that("default priors of efsapriors", {
  # ryegrass: a decreasing curve; lognormal: priors of the medians, the background on b2
  e <- efsapriors(rootl ~ conc, ryegrass, "loglogistic", "lognormal")
  y0 <- exp(mean(log(ryegrass$rootl[ryegrass$conc == 0])))
  sl <- mean(tapply(log(ryegrass$rootl), ryegrass$conc, sd))
  expect_false(e$increasing)
  expect_equal(e$curvelink, "log")
  expect_equal(names(e$prior), c("b1", "b2", "b3", "bmd", "sigma"))
  expect_equal(e$prior$b2$hyp[1:2], c(0.5*y0, 2*y0))
  expect_equal(e$prior$b2$post, 3L)
  expect_equal(e$prior$b3$name, "negate(lognormal(1, 1))")
  expect_equal(e$prior$bmd$hyp[1:2], c(0, 3*30))
  expect_equal(e$prior$sigma$hyp[1:2], c(log(1.5*sl), 0.5))
  expect_equal(e$bmdpar, list(bmr=0.05, def="relative"))
  expect_output(print(e), "decreasing curve")
  # gaussian: the maximum response on b1 (the lower asymptote), on the scale of the response
  eg <- efsapriors(rootl ~ conc, ryegrass, "weibull2", "gaussian", extended=FALSE, dprior="alternative")
  yN <- mean(ryegrass$rootl[ryegrass$conc == 30])
  y0g <- mean(ryegrass$rootl[ryegrass$conc == 0])
  expect_equal(eg$prior$b1$hyp[2], if (eg$flat) y0g*(1 - 0.05 - 0.01) else y0g*(1 - 0.05 - 0.01))
  expect_equal(eg$prior$bmd$hyp[2], 30)
  expect_equal(eg$prior$b3$name, "negate(lognormal(0.4, 0.707107))")
  # quantal background from the controls
  eq <- efsapriors(fq, TCDD, "gamma", "binomial")
  ci <- prop.test(2, 86)$conf.int
  expect_equal(eq$prior$b1$hyp[1:2], c(max(ci[1]/2, 1/860), min(3*ci[2]/2, 1 - 1/860)))
  expect_equal(eq$prior$b3$name, "lognormal(1, 1)")
  expect_equal(eq$start, c(b5=0, b2=1))
  expect_equal(efsapriors(fq, TCDD, "multistage", "binomial")$prior$b3$name, "normal(0, 1)")
  # informative priors
  ei <- efsapriors(fq, TCDD, "gamma", "binomial", background=c(0.01, NA, 0.05), bmd=pert(1, 5, 50))
  expect_equal(ei$prior$b1$hyp[1:2], c(0.01, 0.05))
  expect_equal(ei$prior$bmd$name, "pert(1, 5, 50)")
  # checks
  expect_error(efsapriors(rootl ~ conc, ryegrass[ryegrass$conc > 0, ], "gamma"), "control group")
  expect_error(efsapriors(fq, TCDD, "logistic", "binomial"), "model must be one of")
  dec <- TCDD
  dec$incidence <- rev(dec$incidence)
  expect_error(efsapriors(fq, dec, "gamma", "binomial"), "must increase")
  expect_error(efsapriors(fq, TCDD, "gamma", "binomial", maxresponse=c(0.5, 0.8, 1)), "no maximum response")
  # the curve levels off: last slope below half the largest
  expect_true(fpflat(c(0, 0.1, 0.3, 1), c(1, 5, 9, 9.5), TRUE))
  expect_false(fpflat(c(0, 0.1, 0.3, 1), c(1, 1.2, 2, 9), TRUE))
})

test_that("bmdefsa reproduces the TCDD example of EFSA EN-9902 (Table 5.1)", {
  # EFSA: 8 models (with probit and logistic), bridge sampling: historical (Kociba) BMD 10.63 (4.52,
  # 22.67), current (NTP) 24.67 (20.45, 28.33), with the informative MLE PERT prior 23.87 (19.75, 27.83)
  set.seed(1)
  h <- suppressWarnings(bmdefsa(fq, data=TCDD, family="binomial"))
  expect_s3_class(h, "bmdbayes")
  expect_equal(nrow(h$models), 8)
  expect_length(h$fits, 8)
  expect_equal(h$models["probit_Q", "model"], "probit")
  expect_equal(unname(c(coef(h), confint(h))), c(10.63, 4.52, 22.67), tolerance=0.06)
  # the quadratic exponential model gets the largest weight, as in EFSA
  expect_equal(rownames(h$models)[which.max(h$models$weight)], "multistage_Q")
  set.seed(1)
  cur <- suppressWarnings(bmdefsa(fq, data=liver.tumor, family="binomial"))
  expect_equal(unname(c(coef(cur), confint(cur))), c(24.67, 20.45, 28.33), tolerance=0.01)
  # informative BMD prior from the historical study: EFSA pert(0.980, 9.373, 76.236), shape 10.04
  p <- fitpert(h)
  expect_equal(p$hyp[1:2], c(0.980, 76.236), tolerance=0.1)
  set.seed(1)
  ci <- suppressWarnings(bmdefsa(fq, data=liver.tumor, family="binomial", bmd=p))
  expect_equal(unname(c(coef(ci), confint(ci))), c(23.87, 19.75, 27.83), tolerance=0.01)
  expect_equal(ci$prior, p$name)
})

test_that("fitpert recovers a PERT distribution, discounts and mixes posteriors", {
  set.seed(4)
  # PERT(1, 3, 10) with shape 4: a scaled beta
  x <- 1 + 9*rbeta(1e5, 1 + 4*2/9, 1 + 4*7/9)
  p <- fitpert(x)
  # the bounds are the extreme quantile midpoints (0.05% and 99.95%), as in EN-9902; the fitted PERT has
  # the central quantiles of the distribution
  expect_equal(p$hyp[1:2], unname(quantile(x, c(0.0005, 0.9995))), tolerance=0.01)
  qp <- function(q, pr) q$hyp[1] + (q$hyp[2] - q$hyp[1])*qbeta(pr, q$hyp[3], q$hyp[4])
  expect_equal(qp(p, c(0.05, 0.5, 0.95)), unname(quantile(x, c(0.05, 0.5, 0.95))), tolerance=0.06)
  pd <- fitpert(x, discount=TRUE)
  sh <- function(q) (q$hyp[3] - 1)*(q$hyp[2] - q$hyp[1])/(q$center - q$hyp[1])
  expect_lt(sh(pd), sh(p))
  # weighted draws, and a mixture of two posteriors with equal weights
  pw <- fitpert(x, weights=rep(1, length(x)))
  expect_equal(pw$hyp, p$hyp, tolerance=1e-3)
  pm <- fitpert(x, x + 5)
  expect_equal(pm$hyp[1:2], unname(quantile(c(x, x + 5), c(0.0005, 0.9995))), tolerance=0.01)
  expect_error(fitpert("a"), "bmdbayes objects or numeric")
})
