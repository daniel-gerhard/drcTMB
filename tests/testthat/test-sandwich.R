# robust (sandwich) covariance: against drc with package sandwich (no random effects) and glmmTMB
# (random effects), and the argument vcov. of summary, confint, predict, ed and bmd

data(ryegrass)
data(earthworms)

fitry <- function(...) drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), ...)

test_that("observation and cluster scores match drc with sandwich", {
  skip_if_not_installed("drc")
  skip_if_not_installed("sandwich")
  f <- fitry()
  m <- drc::drm(rootl ~ conc, data=ryegrass, fct=drc::LL.4())
  # drc: b = -b3; the robust covariance of the coefficients of least squares and maximum likelihood agree
  e <- ed(f, c(10, 50), vcov.=vcovrobust)
  em <- drc::ED(m, c(10, 50), vcov.=sandwich::sandwich, display=FALSE)
  expect_equal(unname(sqrt(diag(vcov(e)))), unname(em[, 2]), tolerance=1e-4)
  expect_equal(vcovrobust(f), sandwich::sandwich(f))
  expect_equal(dim(sandwich::estfun(f)), c(nrow(ryegrass), 5L))
  # clusters, without small-sample adjustment
  cl <- rep(1:6, length.out=nrow(ryegrass))
  ec <- ed(f, 50, vcov.=function(x) vcovrobust(x, cluster=cl))
  emc <- drc::ED(m, 50, vcov.=function(x) sandwich::vcovCL(x, cluster=cl, type="HC0", cadjust=FALSE), display=FALSE)
  expect_equal(unname(sqrt(diag(vcov(ec)))), unname(emc[, 2]), tolerance=1e-4)
  # HC1: G/(G - 1) (n - 1)/(n - p)
  expect_equal(vcovrobust(f, cluster=cl, type="HC1"), vcovrobust(f, cluster=cl) * 6/5 * 23/19)
  # binomial
  fb <- drmTMB(cbind(number, total - number) ~ dose, data=earthworms, family="binomial", model="loglogistic",
               start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))
  mb <- drc::drm(number/total ~ dose, weights=total, data=earthworms, fct=drc::LL.3(), type="binomial")
  eb <- ed(fb, 50, vcov.=vcovrobust)
  expect_equal(unname(sqrt(diag(vcov(eb)))), unname(drc::ED(mb, 50, vcov.=sandwich::sandwich, display=FALSE)[, 2]), tolerance=1e-4)
})

test_that("cluster scores of random effect models match glmmTMB", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("sandwich")
  set.seed(1)
  d <- data.frame(g=factor(rep(1:8, each=6)), x=rep(1:6, 8))
  d$y <- 5 + rnorm(8)[d$g] + rnorm(48, sd=0.5*(1 + d$x/3))
  d$cnt <- rpois(48, exp(1 + rnorm(8, sd=0.5)[d$g]))
  f <- drmTMB(y ~ x, rform=b1 ~ 1 | g, data=d, model="constant")
  g <- glmmTMB::glmmTMB(y ~ 1 + (1 | g), data=d)
  # glmmTMB: (Intercept), log sigma, log sd; drcTMB: b1, log sd, log sigma
  expect_equal(unname(vcovrobust(f)), unname(sandwich::vcovHC(g, full=TRUE)[c(1, 3, 2), c(1, 3, 2)]), tolerance=1e-5)
  expect_equal(nrow(sandwich::estfun(f)), 8L)
  fp <- drmTMB(cnt ~ x, rform=b1 ~ 1 | g, data=d, model="constant", family="poisson", parlink=c(b1="log"))
  gp <- glmmTMB::glmmTMB(cnt ~ 1 + (1 | g), data=d, family=poisson)
  expect_equal(unname(vcovrobust(fp)), unname(sandwich::vcovHC(gp, full=TRUE)[c(1, 2), c(1, 2)]), tolerance=1e-5)
})

# repeated measurements of 12 subjects in 4 labs with AR(1) residuals
set.seed(1)
doses <- c(0.5, 1, 2, 4, 8, 16, 32)
dc <- expand.grid(k=seq_along(doses), subject=factor(1:12))
dc$dose <- doses[dc$k]
dc$lab <- factor((as.integer(dc$subject) - 1) %/% 3)
dc$y <- 1 + 7/(1 + exp(1.5*(log(dc$dose) - log(4)))) + unlist(lapply(1:12, function(s) as.numeric(arima.sim(list(ar=0.6), 7, sd=0.5))))

test_that("series of a correlation structure, nested and crossed groups, priors", {
  f <- drmTMB(y ~ dose, data=dc, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE), correlation=cor_ar1(~ k | subject))
  # the scores of the series: as from the fits with the weights of the other series set to zero
  S <- estfun.drmTMB(f)
  expect_equal(nrow(S), 12L)
  env <- f$obj$env
  S2 <- t(sapply(levels(dc$subject), function(s){
    dat <- env$data
    dat$w <- as.numeric(dc$subject == s)
    o <- TMB::MakeADFun(dat, env$parList(f$opt$par), map=env$map, DLL=env$DLL, silent=TRUE)
    -o$gr(f$opt$par)
  }))
  expect_equal(unname(S), unname(S2), tolerance=1e-8)
  expect_equal(unname(colSums(S)), -as.vector(f$obj$gr(f$opt$par)), tolerance=1e-8)
  expect_error(vcovrobust(f, cluster=seq_len(nrow(dc))), "series of the correlation structure must be nested")
  expect_equal(dim(vcovrobust(f, cluster=~ lab)), c(6L, 6L))

  # nested random effects: the clusters are the labs by default
  fr <- drmTMB(y ~ dose, rform=b2 ~ 1 | lab/subject, data=dc, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  expect_equal(nrow(estfun.drmTMB(fr)), 4L)
  expect_error(vcovrobust(fr, cluster=~ subject), "groups of lab must be nested")
  dc$blk <- factor(dc$k %% 2)
  fc <- drmTMB(y ~ dose, rform=list(b2 ~ 1 | subject, b1 ~ 1 | blk), data=dc, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  expect_error(vcovrobust(fc), "not nested in one another")
  expect_error(vcovrobust(fitry(), cluster=rep(1, nrow(ryegrass))), "at least 2 clusters")
  fp <- fitry(prior=list(b2=normal(7, 2)))
  expect_error(vcovrobust(fp), "priors")
})

test_that("vcov. in summary, confint, predict, ed and bmd", {
  f <- fitry()
  V <- vcovrobust(f)
  se <- sqrt(diag(V))
  expect_equal(unname(summary(f, vcov.=vcovrobust)$coefficients[, 2]), unname(se[1:4]))
  expect_equal(summary(f, vcov.=V)$coefficients, summary(f, vcov.=vcovrobust)$coefficients)
  ci <- confint(f, vcov.=V)
  expect_equal(unname(ci[1:4, 2] - ci[1:4, 1]), unname(2*zquant(f, 0.95)*se[1:4]))
  expect_error(confint(f, method="profile", vcov.=V), "only used by Wald")
  # the covariance of the fit gives the default results
  expect_equal(vcov(ed(f, 50, vcov.=vcov)), vcov(ed(f, 50)))
  expect_equal(predict(f, se.fit=TRUE, vcov.=vcov), predict(f, se.fit=TRUE))
  expect_error(ed(f, 50, vcov.=V[1:4, 1:4]), "5 x 5")
  # effective doses with vcov.: Wald intervals by default
  e <- ed(f, 50, vcov.=V)
  expect_equal(confint(e), confint(e, method="wald"))
  b <- bmd(f, 0.1, def="relative", vcov.=V)
  expect_false(isTRUE(all.equal(vcov(b), vcov(bmd(f, 0.1, def="relative")))))

  # random effects: the joint covariance of the fixed and random effects for vcov. of the fit is that of
  # the joint precision of sdreport
  fs <- drmTMB(y ~ dose, rform=b2 ~ 1 | subject, data=dc, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  expect_equal(unname(jointcov(vcovfit(fs, vcov), all=TRUE)), unname(jointcov(fs, all=TRUE)), tolerance=1e-6)
  # a singular sandwich (4 clusters, 7 parameters) still gives standard errors and simulated intervals
  fr <- drmTMB(y ~ dose, rform=b2 ~ 1 | lab/subject, data=dc, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  nd <- data.frame(dose=c(1, 4), subject="1", lab="0")
  p <- predict(fr, newdata=nd, interval="confidence", vcov.=vcovrobust)
  expect_true(all(is.finite(as.matrix(p))))
  p <- predict(fr, newdata=nd, interval="prediction", vcov.=vcovrobust)
  expect_true(all(p$lwr < p$fit & p$fit < p$upr))
  expect_true(is.finite(sqrt(vcov(ed(fr, 50, vcov.=vcovrobust)))))
})

test_that("event data: the series are the clusters", {
  set.seed(1)
  d <- do.call(rbind, lapply(1:10, function(k){
    t <- exp(log(5) + rnorm(1, 0, 0.3)) * (function(u) (u/(1 - u))^(1/4))(runif(30))
    data.frame(dish=k, start=0:14, end=c(1:14, Inf), count=tabulate(findInterval(t, 0:14), nbins=15))
  }))
  f <- drmTMB(count ~ start + end, data=d, model="loglogistic", family="event", start=c(b5=0), fix=c(b5=TRUE), series=~ dish)
  expect_equal(nrow(estfun.drmTMB(f)), 10L)
  expect_equal(vcovrobust(f), vcovrobust(f, cluster=~ dish))
  expect_error(vcovrobust(f, cluster=seq_len(nrow(d))), "series of the event data")
})
