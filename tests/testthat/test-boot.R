# Parametric bootstrap intervals

test_that("bootstrap refits equal drmTMB fits to the simulated responses", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  bs <- bootstrap(f, nsim=3, seed=1)
  expect_s3_class(bs, "bootstrap.drmTMB")
  b <- bs$estimates
  expect_equal(dim(b), c(3, length(coef(f))))
  expect_equal(colnames(b), names(coef(f)))

  sims <- simulate(f, nsim=3, seed=1)
  for (i in 1:3){
    d$y <- sims[[i]]
    fi <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
    expect_equal(b[i, ], coef(fi), tolerance=1e-4)
  }
})

test_that("bootstrap intervals are reproducible and close to Wald intervals for large samples", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic", nrep=20), model="loglogistic")
  expect_message(cb <- confint(f, method="boot", nsim=200, seed=1), "bootstrap")
  expect_identical(cb, suppressMessages(confint(f, method="boot", nsim=200, seed=1)))
  expect_equal(dimnames(cb), dimnames(confint(f, method="wald")))
  # percentile intervals: Monte Carlo error with 200 simulations
  w <- confint(f, method="wald")
  expect_true(all(abs(cb - w) < 0.25*(w[, 2] - w[, 1])))
})

test_that("bootstrap intervals of effective doses use the bootstrap coefficients", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  # ED50 of the four-parameter loglogistic model is b4
  e <- ed(f, respLev=c(50, 90))
  ce <- suppressMessages(confint(e, method="boot", nsim=100, seed=2))
  cb <- suppressMessages(confint(f, parm=4, method="boot", nsim=100, seed=2))
  expect_equal(unname(ce[1, ]), unname(cb[1, ]))
  expect_true(ce[2, 1] < coef(e)[2] && coef(e)[2] < ce[2, 2])

  # effective doses including the variation of the random effects: Wald intervals
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  fr <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  er <- ed(fr)
  expect_message(cr <- confint(er, method="boot"), "Wald intervals")
  expect_equal(cr, confint(er, method="wald"))
})

test_that("a stored bootstrap sample is reused by confint", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  b <- bootstrap(f, nsim=100, seed=3)
  out <- capture.output(print(b))
  expect_true(any(grepl("100 of 100 fits converged", out)))
  expect_true(any(grepl("Bias", out)))

  # the same intervals from the stored sample as from new simulations with the same seed
  ci <- suppressMessages(confint(f, method="boot", nsim=100, seed=3))
  expect_silent(cb <- confint(f, boot=b))
  expect_identical(cb, ci)
  expect_identical(confint(b), ci)
  expect_identical(confint(b, parm=4, level=0.9), confint(f, parm=4, level=0.9, boot=b))

  e <- ed(f, respLev=c(10, 50, 90))
  expect_identical(confint(e, boot=b), suppressMessages(confint(e, method="boot", nsim=100, seed=3)))

  # a bootstrap sample of another model
  f2 <- fit4(y ~ dose, data=simgauss("loglogistic", seed=2), model="loglogistic")
  expect_error(confint(f2, boot=b), "different model")
  expect_error(confint(ed(f2), boot=b), "different model")
})

test_that("bootstrap uses the iteration limits of the fit", {
  d <- simgauss("loglogistic")
  expect_warning(f <- fit4(y ~ dose, data=d, model="loglogistic", control=list(iter.max=1)), "did not converge")
  expect_message(bs <- suppressWarnings(bootstrap(f, nsim=2, seed=1, progress=FALSE)), "2 of 2 bootstrap fits did not converge")
  expect_true(all(is.na(bs$estimates)))
})
