# Likelihood profiles. Profile likelihood intervals are invariant to monotone reparameterisations,
# so they can be compared exactly with other software for equivalent models.

test_that("profile intervals of a gaussian model are likelihood ratio intervals", {
  d <- simgauss("loglogistic")
  d$dose[d$dose == 0] <- 1e-8
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  p <- profile(f)
  expect_s3_class(p, "profile.drmTMB")
  expect_named(p, names(bcoef(f)))
  ci <- confint(p)
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_true(all(ci[, 1] < bcoef(f) & bcoef(f) < ci[, 2]))

  # nls with b4 fixed at a limit: the deviance n log(RSS / RSSmin) equals the F-based cutoff of
  # small-sample inference, n log(1 + F/(n - p)), i.e. the F-test of b4 rejects at the limits
  rss0 <- sum((d$y - predict(f))^2)
  n <- nrow(d)
  fcut <- n*log(1 + qf(0.95, 1, n - 4)/(n - 4))
  for (b4 in ci[4, ]){
    n1 <- nls(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))), data=d,
              start=setNames(bcoef(f)[1:3], paste("b", 1:3, sep="")))
    expect_equal(n*log(sum(residuals(n1)^2)/rss0), fcut, tolerance=1e-3)
  }
  # the profile intervals of nls (MASS), from the t distribution of the profile t statistic
  ns <- nls(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))), data=d,
            start=setNames(bcoef(f), paste("b", 1:4, sep="")))
  expect_equal(unname(ci), unname(suppressMessages(confint(ns))), tolerance=1e-3)

  # maximum likelihood: the deviance equals the chi-square quantile
  fml <- fit4(y ~ dose, data=d, model="loglogistic", control=list(smallsample=FALSE))
  ciml <- confint(profile(fml, parm=4))
  n1 <- nls(y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4)))), data=transform(d, b4=ciml[1, 1]),
            start=setNames(bcoef(f)[1:3], paste("b", 1:3, sep="")))
  expect_equal(n*log(sum(residuals(n1)^2)/rss0), qchisq(0.95, 1), tolerance=1e-3)

  # close to the Wald intervals for a nearly linear parameter
  expect_equal(unname(ci[2, ]), unname(confint(f, method="wald")[2, ]), tolerance=0.02)
})

test_that("profile intervals with random effects agree with lme4", {
  skip_if_not_installed("lme4")
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=2)
  u <- rnorm(12, 0, 0.5)
  d$y <- rpois(nrow(d), exp(log(50) + u[d$grp] - d$dose/4))
  f <- drmTMB(y ~ dose, rform=b1 ~ 1 | grp, data=d, model="exponentialdecay", family="poisson",
              link="log", start=c(b2=-30), fix=c(b2=TRUE))
  ci <- confint(profile(f, parm=names(coef(f))))
  expect_equal(rownames(ci), c("b1:(Intercept)", "b4:(Intercept)", "sd:b1|grp"))

  # log(mu) = b1 - dose / b4: intercept, -1 / slope, and the same standard deviation
  g <- lme4::glmer(y ~ dose + (1 | grp), data=d, family=poisson)
  cg <- suppressMessages(suppressWarnings(confint(g, method="profile")))
  expect_equal(unname(ci[1, ]), unname(cg["(Intercept)", ]), tolerance=1e-4)
  expect_equal(unname(ci[2, ]), unname(sort(-1/cg["dose", ])), tolerance=1e-4)
  expect_equal(unname(ci[3, ]), unname(cg[".sig01", ]), tolerance=1e-4)
})

test_that("standard deviations at the boundary have a lower limit of zero", {
  d <- simgroups(dose=c(0, 1, 2, 4, 8, 16), ng=6, nrep=2, seed=14)
  d$y <- rpois(nrow(d), drcmean("loglogistic", d$dose, c(5, 45, -1.2, 3, 0)))
  f <- suppressWarnings(fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic", family="poisson"))
  ci <- suppressWarnings(confint(profile(f, parm="log_sd:b2|grp")))
  expect_equal(rownames(ci), "sd:b2|grp")
  expect_equal(ci[1, 1], 0)
  expect_true(is.finite(ci[1, 2]) && ci[1, 2] > 0)
})

test_that("intervals beyond the profiled level are missing with a warning", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  p <- profile(f, parm=4, level=0.8)
  expect_warning(ci <- confint(p, level=0.999), "do not reach")
  expect_true(any(is.na(ci)))
  expect_silent(confint(p, level=0.5))
})

test_that("the profiled coefficient stays within its bounds", {
  d <- simgauss("loglogistic")
  f0 <- fit4(y ~ dose, data=d, model="loglogistic")
  up <- c(Inf, Inf, Inf, bcoef(f0)[4] + 0.1, Inf)
  f <- fit4(y ~ dose, data=d, model="loglogistic", upper=up)
  p <- profile(f, parm=4)
  expect_lte(max(p[[1]]$par), up[4] + 1e-10)
  expect_warning(ci <- confint(p), "do not reach")
  expect_true(is.na(ci[1, 2]))
})

test_that("plot draws the profiles with ggplot2", {
  skip_if_not_installed("ggplot2")
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  p <- profile(f, parm=3:4)
  g <- plot(p)
  expect_s3_class(g, "ggplot")
  expect_equal(levels(g$data$parameter), names(p))
  expect_equal(nrow(g$data), sum(sapply(p, nrow)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(print(g))
})

## profiles of effective doses ---------------------------------------------------------

test_that("b4 computed from the effective dose inverts ed() for all models", {
  for (m in models){
    f <- fit4(y ~ dose, data=simgauss(m), model=m)
    e <- ed(f, respLev=c(10, 50, 90), random=FALSE)
    b <- numeric(5)
    b[!f$fix] <- bcoef(f)
    b[f$fix] <- f$start[f$fix]
    b4 <- sapply(1:3, function(i) drcTMB:::edb4(match(m, models), e$q[[i]], e$estimates[[i]], b[3], b[5]))
    expect_equal(b4, rep(b[4], 3), tolerance=1e-10, label=m)
  }
})

test_that("profile intervals of effective doses are likelihood ratio intervals", {
  for (m in c("loglogistic", "logistic")){
    d <- simgauss(m)
    d$dose[d$dose == 0] <- 1e-8
    f <- fit4(y ~ dose, data=d, model=m)
    e <- ed(f, respLev=c(10, 50, 90))
    pe <- profile(e)
    expect_s3_class(pe, "profile.edTMB")
    ci <- confint(pe)
    expect_equal(rownames(ci), names(coef(e)))
    expect_true(all(ci[, 1] < coef(e) & coef(e) < ci[, 2]))

    # nls with the effective dose fixed at a limit: b4 from the effective dose; F-based cutoff
    rss0 <- sum((d$y - predict(f))^2)
    n <- nrow(d)
    fcut <- n*log(1 + qf(0.95, 1, n - 4)/(n - 4))
    mod <- match(m, models)
    for (j in 1:3){
      q <- e$q[[j]]
      for (ED in ci[j, ]){
        fm <- if (m == "loglogistic") {
          y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(drcTMB:::edb4(mod, q, ED, b3, 0)))))
        } else {
          y ~ b1 + (b2 - b1)/(1 + exp(-b3*(dose - drcTMB:::edb4(mod, q, ED, b3, 0))))
        }
        n1 <- nls(fm, data=d, start=setNames(bcoef(f)[1:3], paste("b", 1:3, sep="")))
        expect_equal(n*log(sum(residuals(n1)^2)/rss0), fcut, tolerance=1e-4,
                     label=paste(m, names(coef(e))[j]))
      }
    }
  }
})

test_that("the ED50 profile of the four-parameter loglogistic model is the profile of b4", {
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  expect_equal(unname(confint(profile(ed(f)))), unname(confint(profile(f, parm=4))), tolerance=1e-4)
})

test_that("effective doses of groups are profiled with random effects", {
  d <- simgroups(ng=10, seed=20)
  d$trt <- factor(rep(c("A", "B"), length.out=nrow(d)))
  u <- rnorm(10, 0, 0.3)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  d$y[d$trt == "B"] <- d$y[d$trt == "B"] - 1
  f <- fit4(y ~ dose, fform=b4 ~ trt, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  # coefficients b1, b2, b3, b4:(Intercept), b4:trtB, b5
  lA <- cbind(diag(5)[, 1:4], 0, diag(5)[, 5])
  lB <- lA
  lB[4, 5] <- 1
  e <- ed(f, linfct=list(A=lA, B=lB), random=FALSE)
  ci <- confint(profile(e))
  # ED50 of group B is b4:(Intercept) + b4:trtB, profiled as a linear combination by TMB
  lc <- TMB::tmbprofile(f$obj, lincomb=c(0, 0, 0, 1, 1, 0, 0), ytol=3.5, trace=FALSE)
  expect_equal(unname(ci["50:B", ]), unname(confint(lc)[1, ]), tolerance=1e-4)
})

test_that("profiles of effective doses agree with glmer at fixed effective doses", {
  skip_if_not_installed("lme4")
  d <- simgroups(dose=c(0, 1, 2, 4, 6, 8, 10), ng=12, nrep=1, seed=2)
  u <- rnorm(12, 0, 0.8)
  d$n <- 20
  d$y <- rbinom(nrow(d), d$n, plogis(0.9*(d$dose - (5 + u[d$grp]))))
  f <- drmTMB(cbind(y, n - y) ~ dose, rform=b4 ~ 1 | grp, data=d, model="logistic", family="binomial",
              start=c(b1=0, b2=1, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE))
  e <- ed(f, respLev=c(10, 90), random=FALSE)
  ci <- confint(profile(e))
  # at ED = t for level q: logit(p) = b3 (dose - t - u) + logit(q), a GLMM with an offset
  for (j in 1:2){
    for (t in ci[j, ]){
      d$off <- qlogis(e$q[[j]])
      g <- lme4::glmer(cbind(y, n - y) ~ 0 + I(dose - t) + (1 | grp), offset=off, data=d, family=binomial)
      expect_equal(2*(-as.numeric(logLik(g)) - f$opt$objective), qchisq(0.95, 1), tolerance=1e-3)
    }
  }
})

test_that("profiles of effective doses can be plotted", {
  skip_if_not_installed("ggplot2")
  f <- fit4(y ~ dose, data=simgauss("loglogistic"), model="loglogistic")
  g <- plot(profile(ed(f, respLev=c(10, 90))))
  expect_s3_class(g, "ggplot")
})

## confint methods ------------------------------------------------------------------

test_that("confint of drmTMB gives profile intervals by default and Wald intervals on request", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")

  expect_message(ci <- confint(f, parm=c(4, 5)), "profile")
  expect_equal(ci, confint(profile(f, parm=c(4, 5), level=0.95)))
  expect_equal(rownames(ci), c("b4:(Intercept)", "sd:b2|grp"))

  w <- confint(f, method="wald", level=0.9)
  expect_equal(colnames(w), c("5 %", "95 %"))
  expect_equal(rownames(w), c(names(bcoef(f)), "sd:b2|grp", "sigma"))
  est <- coef(f)
  se <- sqrt(diag(vcov(f)))
  z <- qnorm(0.95)
  expect_equal(unname(w[1:4, ]), unname(cbind(est - z*se, est + z*se)[1:4, ]))
  # standard deviations: Wald intervals on the log scale
  expect_equal(unname(w[5, ]), unname(exp(est[5] + c(-z, z)*se[5])))
  expect_silent(confint(f, method="wald"))
})

test_that("confint of correlated random effects reports the correlation", {
  d <- simgroups(ng=20, seed=21)
  u <- matrix(rnorm(40, 0, 0.5), 20)
  d$y <- sapply(seq_len(nrow(d)), function(i){
    drcmean("loglogistic", d$dose[i], c(1 + u[d$grp[i], 1], 8 + u[d$grp[i], 2], -1.5, 3, 0))
  }) + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b1 + b2 ~ 1 | grp, data=d, model="loglogistic")
  w <- confint(f, parm="corr_par:b2,b1|grp", method="wald")
  expect_equal(rownames(w), "corr:b2,b1|grp")
  expect_true(all(w >= -1 & w <= 1))
})

test_that("confint of effective doses uses profiles for population effective doses", {
  d <- simgroups(ng=8, seed=7)
  u <- rnorm(8, 0, 0.5)
  d$y <- drcmean("loglogistic", d$dose, c(1, 8, -1.5, 3, 0)) + u[d$grp] + rnorm(nrow(d), 0, 0.3)
  f <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")

  e <- ed(f, respLev=c(10, 50), random=FALSE)
  expect_message(ci <- confint(e), "Computing profile")
  expect_equal(ci, confint(profile(e, level=0.95)))
  w <- confint(e, method="wald")
  se <- sqrt(diag(vcov(e)))
  expect_equal(unname(w), unname(cbind(coef(e) - qnorm(0.975)*se, coef(e) + qnorm(0.975)*se)))

  # with the variation of the random effects: Wald intervals
  er <- ed(f, respLev=c(10, 50))
  expect_message(cr <- confint(er), "Wald intervals")
  expect_equal(cr, confint(er, method="wald"))

  # without random effects, random=TRUE makes no difference
  f0 <- fit4(y ~ dose, data=d, model="loglogistic")
  expect_message(confint(ed(f0)), "Computing profile")
})
