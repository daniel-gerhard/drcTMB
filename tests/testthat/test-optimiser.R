# optimiser choice (nlminb, optim), control parameters of the optimiser and of the inner optimisation,
# and the restart at the estimates

data(ryegrass)

fitry <- function(control=list(), ...) drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", control=control, ...)

set.seed(1)
doses <- c(0.5, 1, 2, 4, 8, 16, 32)
dr <- expand.grid(k=seq_along(doses), subject=factor(1:12))
dr$dose <- doses[dr$k]
dr$y <- 1 + 7/(1 + exp(1.5*(log(dr$dose) - log(4)))) + rnorm(12, sd=0.5)[dr$subject] + rnorm(84, sd=0.4)
fitre <- function(control=list()) drmTMB(y ~ dose, rform=b2 ~ 1 | subject, data=dr, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), control=control)

test_that("optim reaches the optimum of nlminb", {
  f0 <- fitry()
  for (ctl in list(list(optimiser="optim"), list(optimiser="optim", method="BFGS"), list(optimiser="optim", method="Nelder-Mead", iter.max=5000))){
    f <- fitry(ctl)
    expect_equal(f$opt$convergence, 0)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f0)), tolerance=1e-8)
    expect_equal(coef(f), coef(f0), tolerance=1e-3)
  }
  # bounds (default lower bound 0 of nec) need L-BFGS-B
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="nec", control=list(optimiser="optim", method="BFGS")), "L-BFGS-B")
  fn <- drmTMB(rootl ~ conc, data=ryegrass, model="nec", control=list(optimiser="optim"))
  expect_equal(as.numeric(logLik(fn)), as.numeric(logLik(drmTMB(rootl ~ conc, data=ryegrass, model="nec"))), tolerance=1e-6)
  # BMD parameterisation: an infinite negative log-likelihood for b4 < 0
  bp <- list(start=c(b5=0), fix=c(b5=TRUE), bmdpar=list(bmr=0.1, def="relative"))
  fb <- do.call(fitry, c(list(control=list(optimiser="optim")), bp))
  expect_equal(coef(fb)[["log_bmd"]], coef(do.call(fitry, bp))[["log_bmd"]], tolerance=1e-4)
  # random effects
  fr0 <- fitre()
  fr <- fitre(list(optimiser="optim"))
  expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fr0)), tolerance=1e-8)
  expect_equal(coef(fr), coef(fr0), tolerance=1e-4)
  expect_match(capture.output(print(summary(fr))), "Optimiser: optim \\(L-BFGS-B\\)", all=FALSE)
})

test_that("control parameters of the optimiser and of the inner optimisation", {
  # optcontrol replaces iter.max
  expect_warning(fitry(list(optcontrol=list(iter.max=1))), "did not converge")
  expect_warning(fitry(list(optimiser="optim", optcontrol=list(maxit=1))), "did not converge")
  # a looser tolerance stops earlier
  f <- fitry(list(optcontrol=list(rel.tol=1e-3)))
  expect_lt(f$opt$iterations, fitry()$opt$iterations)
  expect_error(drmTMBcontrol(optcontrol=list(reltol=1e-8)), "unknown control parameters of nlminb")
  expect_error(drmTMBcontrol(optimiser="optim", optcontrol=list(rel.tol=1e-8)), "unknown control parameters of optim")
  expect_error(drmTMBcontrol(optimiser="optim", optcontrol=list(fnscale=-1)), "fnscale")
  expect_error(drmTMBcontrol(optcontrol=list(1)), "named list")
  expect_error(drmTMBcontrol(restart=NA), "TRUE or FALSE")
  # inner optimisation of the random effects
  fr <- fitre(list(inner=list(maxit=200, tol10=1e-6)))
  expect_equal(fr$obj$env$inner.control$maxit, 200)
  expect_equal(fr$obj$env$inner.control$tol10, 1e-6)
  expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fitre())), tolerance=1e-8)
  expect_error(fitre(list(inner=list(foo=1))), "unknown inner control parameters: foo")
  # fits from older versions without the new control parameters
  f0 <- fitry()
  f0$control <- f0$control[c("optimiser", "iter.max", "eval.max", "smallsample")]
  expect_equal(fitcontrol(f0), drmTMBcontrol())
})

test_that("restart at the estimates", {
  f0 <- fitry()
  f <- fitry(list(restart=TRUE))
  expect_equal(coef(f), coef(f0), tolerance=1e-6)
  expect_named(f$opt$restart, c("first", "restart"))
  expect_match(capture.output(print(summary(f))), "Restarted at the estimates: no better optimum", all=FALSE)
  # a fit stopped by the iteration limit is completed by the restart, which keeps the limit: two
  # iterations reach further than one
  f1 <- suppressWarnings(fitry(list(iter.max=1)))
  f2 <- suppressWarnings(fitry(list(iter.max=1, restart=TRUE)))
  expect_lt(f2$opt$objective, f1$opt$objective)
  expect_match(capture.output(print(summary(f2))), "decreased by", all=FALSE)
  # random effects, and the bootstrap fits with the settings of the fit
  fr <- fitre(list(optimiser="optim", restart=TRUE))
  expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fitre())), tolerance=1e-8)
  b <- bootstrap(fitry(list(optimiser="optim", restart=TRUE)), nsim=3, seed=1, progress=FALSE)
  expect_true(all(!is.na(b$estimates)))
})
