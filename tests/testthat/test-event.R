# time-to-event (germination) data: family="event", against drc's type="event", and random effects of
# dishes by simulation

data(chickweed)
data(germination)

fitev <- function(model="loglogistic", data=chickweed, ...){
  st <- if (model == "loglogistic") list(start=c(b5=0), fix=c(b5=TRUE)) else list()
  do.call(drmTMB, c(list(count ~ start + end, data=data, model=model, family="event"), st, list(...)))
}

test_that("chickweed: the event models of drc", {
  skip_if_not_installed("drc")
  f <- fitev()
  m <- drc::drm(count ~ start + end, data=chickweed, fct=drc::LL.3(), type="event")
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance=1e-8)
  # b2 = d on the logit scale, b3 = -b, b4 = e
  cf <- coef(f)
  expect_equal(c(plogis(cf[[1]]), cf[[2]], cf[[3]]), unname(coef(m)[c(2, 1, 3)])*c(1, -1, 1), tolerance=1e-5)
  expect_equal(unname(sqrt(diag(vcov(f)))[2:3]), unname(sqrt(diag(vcov(m)))[c(1, 3)]), tolerance=1e-4)
  e <- ed(f, c(10, 50, 90))
  em <- drc::ED(m, c(10, 50, 90), display=FALSE)
  expect_equal(unname(coef(e)), unname(em[, 1]), tolerance=1e-5)
  expect_equal(unname(sqrt(diag(vcov(e)))), unname(em[, 2]), tolerance=1e-4)
  for (pr in list(c("weibull1", "W1.3"), c("weibull2", "W2.3"), c("lognormal", "LN.3"))){
    g <- fitev(pr[1])
    mg <- drc::drm(count ~ start + end, data=chickweed, fct=getExportedValue("drc", pr[2])(), type="event")
    expect_equal(as.numeric(logLik(g)), as.numeric(logLik(mg)), tolerance=1e-6)
  }
})

test_that("absolute times, predictions, residuals and simulation", {
  f <- fitev()
  b2 <- plogis(coef(f)[[1]])
  # the time at which 10% of all seeds have germinated: the relative time of 10/b2 percent
  ea <- ed(f, 0.1, type="absolute")
  expect_equal(unname(coef(ea)), unname(coef(ed(f, 10/b2))), tolerance=1e-8)
  expect_warning(ea2 <- ed(f, 0.3, type="absolute"), "do not reach")
  expect_true(is.na(coef(ea2)))
  # the cumulative proportion b2 F(t), from the time in start alone
  p <- predict(f, newdata=data.frame(start=c(0, 196.0530, 1e6)))
  expect_equal(p, c(0, b2/2, b2), tolerance=1e-4)
  pc <- predict(f, newdata=data.frame(start=200), interval="confidence")
  expect_true(pc$lwr < pc$fit && pc$fit < pc$upr)
  expect_error(predict(f, newdata=data.frame(start=200), interval="prediction"), "not available for event data")
  # residuals against the multinomial expectations of the series total
  N <- sum(chickweed$count)
  pr <- eventprob(f)
  expect_equal(sum(pr), 1)
  expect_equal(residuals(f), (chickweed$count - N*pr)/sqrt(N*pr*(1 - pr)))
  expect_equal(residuals(f, type="response"), chickweed$count/N - pr)
  s <- summary(f)
  expect_equal(unname(s$dispersion[c("chisq", "df")]), c(sum(residuals(f)^2), nrow(chickweed) - 1 - 3))
  # simulation: multinomial counts with the series total
  sim <- simulate(f, nsim=400, seed=1)
  expect_true(all(colSums(sim) == N))
  expect_equal(rowMeans(sim), N*pr, tolerance=0.05)
  # incomplete series (the last interval left out) have no simulation, residuals or dispersion
  fi <- fitev(data=chickweed[-nrow(chickweed), ])
  expect_error(simulate(fi), "complete series")
  expect_error(residuals(fi), "complete series")
  expect_null(summary(fi)$dispersion)
})

test_that("several curves: germination of rice and wheat, as drc with curveid", {
  skip_if_not_installed("drc")
  g <- droplevels(subset(germination, species %in% c("rice", "wheat") & temp %in% c(16, 22, 28, 34)))
  g$curve <- interaction(g$species, g$temp)
  f <- drmTMB(germinated ~ start + end, fform=list(b2 ~ 0 + curve, b3 ~ 0 + curve, b4 ~ 0 + curve), data=g, model="loglogistic",
              family="event", start=c(b5=0), fix=c(b5=TRUE))
  m <- drc::drm(germinated ~ start + end, curveid=curve, data=g, fct=drc::LL.3(), type="event")
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance=1e-6)
  # the series by default: the curves
  expect_equal(nlevels(f$series), 8)
  e <- ed(f, 50, newdata=data.frame(curve=levels(g$curve)))
  em <- drc::ED(m, 50, display=FALSE)
  em <- em[paste("e:", levels(g$curve), ":50", sep=""), ]
  expect_equal(unname(coef(e)), unname(em[, 1]), tolerance=1e-4)
  expect_equal(unname(sqrt(diag(vcov(e)))), unname(em[, 2]), tolerance=1e-3)
})

test_that("random effects of dishes: estimates, series and simulation", {
  # 12 dishes of 30 seeds for each of two treatments, germinating proportions 0.8 and 0.6, b3 = 4,
  # b4 = 5 and 8 days with a log-normal dish effect of standard deviation 0.3; daily counts to day 14
  set.seed(1)
  d <- do.call(rbind, lapply(c("A", "B"), function(tr) do.call(rbind, lapply(1:12, function(k){
    t <- exp(log(if (tr == "A") 5 else 8) + rnorm(1, 0, 0.3)) * (function(u) (u/(1 - u))^(1/4))(runif(30))
    t[runif(30) > (if (tr == "A") 0.8 else 0.6)] <- Inf
    data.frame(trt=tr, dish=paste(tr, k), start=0:14, end=c(1:14, Inf), count=tabulate(findInterval(t, 0:14), nbins=15))
  }))))
  f <- drmTMB(count ~ start + end, fform=list(b2 ~ trt, b4 ~ trt), rform=b4 ~ 1 | dish, parlink=c(b4="log"), data=d,
              model="loglogistic", family="event", start=c(b5=0), fix=c(b5=TRUE))
  truth <- c(qlogis(0.8), qlogis(0.6) - qlogis(0.8), 4, log(5), log(8/5))
  se <- sqrt(diag(vcov(f)))[1:5]
  expect_true(all(abs(coef(f)[1:5] - truth) < 3*se))
  expect_equal(exp(coef(f)[["log_sd:b4|dish"]]), 0.3, tolerance=0.4)
  # the series by default: the dishes; an explicit series gives the same fit
  expect_equal(nlevels(f$series), 24)
  fs <- drmTMB(count ~ start + end, fform=list(b2 ~ trt, b4 ~ trt), rform=b4 ~ 1 | dish, parlink=c(b4="log"), data=d,
               model="loglogistic", family="event", start=c(b5=0), fix=c(b5=TRUE), series=~ dish)
  expect_equal(logLik(fs), logLik(f))
  expect_equal(as.integer(table(fs$series)), rep(15L, 24))
  # new random effects in each simulation; the dish totals are kept
  sim <- simulate(f, nsim=2, seed=2)
  expect_true(all(tapply(sim[, 1], d$dish, sum) == 30))
  # the times of the population curves
  e <- ed(f, 50, newdata=data.frame(trt=c("A", "B")), random=FALSE)
  expect_equal(unname(coef(e)), exp(c(coef(f)[[4]], coef(f)[[4]] + coef(f)[[5]])), tolerance=1e-6)
})

test_that("checks of event data", {
  expect_error(drmTMB(count ~ start + end, data=chickweed, model="logistic", family="event"), "distribution function of the time")
  expect_error(drmTMB(count ~ start, data=chickweed, model="loglogistic", family="event"), "count ~ start \\+ end")
  bad <- chickweed
  bad$end[1] <- 0
  expect_error(fitev(data=bad), "0 <= start < end")
  bad <- chickweed
  bad$count[1] <- 0.5
  expect_error(fitev(data=bad), "non-negative integers")
  expect_error(drmTMB(count ~ start + end, data=chickweed, model="loglogistic", family="event", start=c(b1=0.1, b5=0), fix=c(b1=TRUE, b5=TRUE)),
               "b1 of event data is fixed at 0")
  expect_error(fitev(bmdpar=list(bmr=0.1, def="relative")), "not available for event data")
  expect_error(bmd(fitev()), "not available for event data")
  expect_error(drmTMB(rootl ~ conc, data=drcTMB::ryegrass, model="loglogistic", series=~ conc), "only used by the event family")
})
