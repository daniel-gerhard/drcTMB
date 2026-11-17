# residual correlation structures and variance functions, against gnls and nlme of nlme

# repeated measurements of 12 subjects at 7 doses: AR(1) residuals along the doses, a larger residual
# standard deviation at the high doses, positions in time and in the plane
set.seed(1)
ns <- 12
doses <- c(0.5, 1, 2, 4, 8, 16, 32)
dc <- expand.grid(k=seq_along(doses), subject=factor(1:ns))
dc$dose <- doses[dc$k]
dc$time <- dc$k + runif(nrow(dc), -0.3, 0.3)
dc$px <- runif(nrow(dc), 0, 10)
dc$py <- runif(nrow(dc), 0, 10)
muc <- 1 + 7/(1 + exp(1.5*(log(dc$dose) - log(4))))
ec <- unlist(lapply(1:ns, function(s) as.numeric(arima.sim(list(ar=0.6), length(doses), sd=0.5))))
dc$strat <- factor(ifelse(dc$dose > 3, "high", "low"), levels=c("low", "high"))
dc$y <- muc + ifelse(dc$strat == "high", 2, 1)*ec
ll4 <- y ~ b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4))))

fitc <- function(correlation=NULL, variance=NULL, data=dc, ...){
  drmTMB(y ~ dose, data=data, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), correlation=correlation, variance=variance, ...)
}
gnlsc <- function(f, correlation=NULL, weights=NULL){
  nlme::gnls(ll4, data=dc, start=setNames(coef(f)[1:4], c("b1", "b2", "b3", "b4")), correlation=correlation, weights=weights,
             control=nlme::gnlsControl(maxIter=500, msMaxIter=500))
}

test_that("correlation structures match gnls", {
  skip_if_not_installed("nlme")
  cases <- list(
    list(cor_ar1(~ k | subject), nlme::corAR1(form=~ k | subject)),
    list(cor_car1(~ time | subject), nlme::corCAR1(form=~ time | subject)),
    list(cor_cs(~ 1 | subject), nlme::corCompSymm(form=~ 1 | subject)),
    list(cor_exp(~ px + py | subject), nlme::corExp(form=~ px + py | subject)),
    list(cor_exp(~ time | subject, nugget=TRUE), nlme::corExp(form=~ time | subject, nugget=TRUE)),
    list(cor_gaus(~ time | subject), nlme::corGaus(form=~ time | subject)),
    list(cor_spher(~ time | subject), nlme::corSpher(3, form=~ time | subject)))
  for (cs in cases){
    f <- fitc(cs[[1]])
    g <- gnlsc(f, correlation=cs[[2]])
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance=1e-6)
    expect_equal(unname(coef(f)[1:4]), unname(coef(g)), tolerance=1e-4)
    # the correlation parameters on the natural scale
    expect_equal(unname(resnatural(f)), unname(coef(g$modelStruct$corStruct, unconstrained=FALSE)), tolerance=1e-3)
  }
  # the order of the observations is the time of the AR1 structure without positions
  expect_equal(logLik(fitc(cor_ar1(~ 1 | subject))), logLik(fitc(cor_ar1(~ k | subject))))
  # the order of the data does not matter
  o <- sample(nrow(dc))
  expect_equal(as.numeric(logLik(fitc(cor_car1(~ time | subject), data=dc[o, ]))), as.numeric(logLik(fitc(cor_car1(~ time | subject)))), tolerance=1e-8)
})

test_that("variance functions match gnls", {
  skip_if_not_installed("nlme")
  cases <- list(
    list(var_ident(~ 1 | strat), nlme::varIdent(form=~ 1 | strat)),
    list(var_power(~ dose), nlme::varPower(form=~ dose)),
    list(var_exp(~ dose), nlme::varExp(form=~ dose)))
  for (vf in cases){
    f <- fitc(variance=vf[[1]])
    g <- gnlsc(f, weights=vf[[2]])
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance=1e-6)
    expect_equal(unname(coef(f)[1:4]), unname(coef(g)), tolerance=1e-4)
  }
  fi <- fitc(variance=var_ident(~ 1 | strat))
  gi <- gnlsc(fi, weights=nlme::varIdent(form=~ 1 | strat))
  expect_equal(unname(resnatural(fi)), unname(coef(gi$modelStruct$varStruct, unconstrained=FALSE)), tolerance=1e-4)
  # together with a correlation structure
  f <- fitc(cor_car1(~ time | subject), var_ident(~ 1 | strat))
  g <- gnlsc(f, nlme::corCAR1(form=~ time | subject), nlme::varIdent(form=~ 1 | strat))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance=1e-6)
  # a power of the fitted mean: the likelihood with the standard deviation sigma |mu|^delta
  fp <- fitc(variance=var_power())
  mu <- fitted(fp)
  s <- exp(coef(fp)[["log_sigma"]])*abs(mu)^coef(fp)[["var_par:delta"]]
  expect_equal(as.numeric(logLik(fp)), sum(dnorm(dc$y, mu, s, log=TRUE)), tolerance=1e-8)
})

test_that("random effects with a correlation structure match nlme; lognormal matches gnls of log(y)", {
  skip_if_not_installed("nlme")
  set.seed(3)
  dr <- dc
  dr$y <- dr$y + rnorm(ns, 0, 0.8)[dr$subject]
  f <- fitc(cor_car1(~ time | subject), data=dr, rform=b2 ~ 1 | subject)
  # random intercepts on b2, which enters the mean linearly: the likelihood of nlme is exact
  g <- nlme::nlme(ll4, fixed=b1 + b2 + b3 + b4 ~ 1, random=b2 ~ 1 | subject, data=dr, method="ML",
                  start=setNames(coef(f)[1:4], c("b1", "b2", "b3", "b4")), correlation=nlme::corCAR1(form=~ time | subject))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance=1e-6)
  expect_equal(unname(coef(f)[1:4]), unname(nlme::fixef(g)), tolerance=1e-4)

  dl <- dc
  dl$y <- exp(dc$y/3)
  fl <- fitc(cor_ar1(~ k | subject), var_exp(~ dose), data=dl, family="lognormal")
  gl <- nlme::gnls(log(y) ~ log(b1 + (b2 - b1)/(1 + exp(-b3*(log(dose) - log(b4))))), data=dl, start=setNames(coef(fl)[1:4], c("b1", "b2", "b3", "b4")),
                   correlation=nlme::corAR1(form=~ k | subject), weights=nlme::varExp(form=~ dose))
  expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(gl)) - sum(log(dl$y)), tolerance=1e-6)
  expect_equal(unname(coef(fl)[1:4]), unname(coef(gl)), tolerance=1e-4)
})

test_that("residuals, the likelihood in R, simulation and prediction intervals", {
  f <- fitc(cor_car1(~ time | subject), var_ident(~ 1 | strat))
  # pearson residuals divided by sigma g; the normalized (quantile) residuals give the log-likelihood
  sg <- exp(coef(f)[["log_sigma"]])*ifelse(dc$strat == "high", exp(coef(f)[["var_par:high"]]), 1)
  expect_equal(residuals(f), (dc$y - fitted(f))/sg)
  z <- residuals(f, type="quantile")
  ldet <- sum(vapply(split(seq_len(nrow(dc)), dc$subject), function(i){
    sum(log(diag(chol(cormatrix(f, as.matrix(dc$time[i]))))))
  }, numeric(1)))
  expect_equal(as.numeric(logLik(f)), sum(dnorm(z, log=TRUE)) - sum(log(sg)) - ldet, tolerance=1e-8)

  # simulated residuals: lag-one correlation phi^(time difference) within subjects
  sim <- simulate(f, nsim=400, seed=1)
  e <- (as.matrix(sim) - fitted(f))/sg
  i1 <- which(dc$k < 7)
  phi <- resnatural(f)[["cor:phi"]]
  expect_equal(mean(e[i1, ]*e[i1 + 1, ]), mean(phi^(dc$time[i1 + 1] - dc$time[i1])), tolerance=0.1)
  expect_equal(sd(as.vector(e)), 1, tolerance=0.05)

  # prediction intervals with the standard deviation sigma g of new observations
  nd <- data.frame(dose=c(1, 10), strat=factor(c("low", "high")))
  p <- predict(f, newdata=nd, interval="prediction")
  s2 <- exp(coef(f)[["log_sigma"]])*c(1, exp(coef(f)[["var_par:high"]]))
  expect_equal(p$upr, p$fit + qnorm(0.975)*sqrt(p$se.fit^2 + s2^2))
  expect_error(predict(f, newdata=data.frame(dose=1), interval="prediction"))
  fp <- fitc(variance=var_power())
  pp <- predict(fp, newdata=data.frame(dose=c(1, 10)), interval="prediction")
  s3 <- exp(coef(fp)[["log_sigma"]])*abs(pp$fit)^coef(fp)[["var_par:delta"]]
  expect_equal(pp$upr, pp$fit + qnorm(0.975)*sqrt(pp$se.fit^2 + s3^2))
})

test_that("printing, intervals on the natural scale, model comparison", {
  f <- fitc(cor_car1(~ time | subject), var_ident(~ 1 | strat))
  expect_equal(names(coef(f))[6:7], c("cor_par:phi", "var_par:high"))
  expect_output(print(f), "Residual correlation: CAR1, ~time \\| subject: phi = ")
  expect_output(print(summary(f)), "Variance function: separate standard deviations, ~1 \\| strat: sd ratio high = ")
  expect_output(print(cor_exp(~ x + y | plot, nugget=TRUE)), "exponential with nugget, ~x \\+ y \\| plot")
  expect_output(print(var_power()), "power, ~fitted")
  ci <- confint(f, method="wald")
  expect_equal(rownames(ci)[6:7], c("cor:phi", "var:sd ratio high"))
  se <- sqrt(diag(vcov(f)))[["cor_par:phi"]]
  expect_equal(unname(ci["cor:phi", ]), plogis(coef(f)[["cor_par:phi"]] + c(-1, 1)*qnorm(0.975)*se))
  cp <- suppressMessages(confint(f, parm="cor_par:phi"))
  expect_true(cp[1] < resnatural(f)[["cor:phi"]] && resnatural(f)[["cor:phi"]] < cp[2])
  f0 <- fitc(variance=var_ident(~ 1 | strat))
  a <- anova(f0, f)
  expect_equal(a$Df[2], 1)
  # no small-sample inference
  expect_null(df.residual(f))
  b <- bootstrap(f, nsim=10, seed=1, progress=FALSE)
  expect_equal(rownames(confint(b, parm="cor_par:phi")), "cor:phi")
})

test_that("checks of the correlation structures and variance functions", {
  expect_error(drmTMB(y ~ dose, data=transform(dc, y=rpois(nrow(dc), 3)), model="loglogistic", family="poisson", correlation=cor_ar1()),
               "gaussian and lognormal")
  expect_error(fitc(cor_ar1(~ k | subject), weights=rep(2, nrow(dc))), "weights cannot be combined")
  expect_error(fitc(cor_car1(~ dose)), "distinct")
  expect_error(fitc(cor_ar1(~ time | subject)), "integers")
  expect_error(fitc(cor_car1(~ px + py | subject)), "single time variable")
  expect_error(fitc(cor_exp(~ 1 | subject)), "needs the positions")
  expect_error(fitc(correlation=~ time), "correlation must be a correlation structure")
  expect_error(cor_car1(y ~ time), "one-sided")
  expect_error(fitc(variance=var_ident(~ strat)), "~ 1 \\| stratum")
  expect_error(fitc(variance=var_power(~ k - 1)), "covariate of the power variance function is 0")
  expect_error(fitc(variance=var_power(~ dose | strat)), "no strata")
  fp <- fitc(variance=var_power())
  expect_error(ed(fp, 50, residual=TRUE), "variance function")
  expect_error(bmd(fp, 1, "sd"), "variance function")
  expect_error(fitc(variance=var_power(), bmdpar=list(bmr=1, def="sd")), "variance function")
})
