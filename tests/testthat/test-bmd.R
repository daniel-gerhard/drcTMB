# quantal data: background 0.1, loglogistic increase
simquantal <- function(seed=2){
  set.seed(seed)
  d <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16), each=3), n=20)
  p <- 0.1 + 0.9*plogis(2*(log(pmax(d$dose, 1e-9)) - log(4)))
  d$y <- rbinom(nrow(d), d$n, p)
  d
}

# grouped gaussian data with random effects on the upper asymptote and log b4
simbmdgroups <- function(seed=3){
  set.seed(seed)
  d <- simgroups(ng=12)
  u2 <- rnorm(12, 0, 0.8)
  u4 <- rnorm(12, 0, 0.5)
  d$y <- with(d, 1 + (8 + u2[grp] - 1)/(1 + exp(1.5*(log(pmax(dose, 1e-12)) - log(3) - u4[grp])))) + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("relative BMD equals the effective dose of a curve to zero", {
  d <- simgauss("loglogistic")
  f <- drmTMB(y ~ dose, data=d, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))
  b <- bmd(f, c(0.05, 0.1), def="relative")
  e <- ed(f, c(5, 10))
  expect_equal(unname(coef(b)), unname(coef(e)), tolerance=1e-8)
  expect_equal(unname(vcov(b)), unname(vcov(e)), tolerance=1e-8)
  expect_equal(unname(confint(b, parm=2, level=0.9, progress=FALSE)), unname(confint(e, parm=2, level=0.9, progress=FALSE)), tolerance=1e-6)
})

test_that("excess risk BMD equals the effective dose of a quantal curve to one", {
  d <- simquantal()
  f <- drmTMB(cbind(y, n - y) ~ dose, data=d, family="binomial", model="loglogistic", start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  b <- bmd(f, 0.1)
  e <- ed(f, 10)
  expect_equal(unname(coef(b)), unname(coef(e)), tolerance=1e-8)
  expect_equal(unname(vcov(b)), unname(vcov(e)), tolerance=1e-8)
  # additional risk: the dose at which the curve reaches the background plus bmr
  ba <- bmd(f, 0.1, def="additional")
  p0 <- coef(f)[["b1:(Intercept)"]]
  expect_equal(unname(coef(ba)), unname(coef(ed(f, p0 + 0.1, type="absolute"))), tolerance=1e-8)
  # multistage: excess risk is the level of F
  fm <- drmTMB(cbind(y, n - y) ~ dose, data=d, family="binomial", model="multistage", start=c(b2=1), fix=c(b2=TRUE))
  expect_equal(unname(coef(bmd(fm, 0.1))), unname(coef(ed(fm, 10))), tolerance=1e-8)
})

test_that("BMDs and BMDLs agree with the bmd package", {
  skip_if_not_installed("drc")
  skip_if_not_installed("bmd")
  # bmd recognises the drc model functions by their name, which must not carry the drc:: prefix
  W1.4 <- drc::W1.4
  LL.4 <- drc::LL.4
  d <- simquantal()
  d$p <- d$y/d$n
  fd <- drc::drm(p ~ dose, weights=n, data=d, type="binomial", fct=W1.4(fixed=c(NA, NA, 1, NA)))
  ft <- drmTMB(cbind(y, n - y) ~ dose, data=d, family="binomial", model="weibull1", start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))
  for (def in c("excess", "additional")){
    rd <- bmd::bmd(fd, 0.1, def=def, backgType="modelBased", display=FALSE)
    rp <- bmd::bmd(fd, 0.1, def=def, backgType="modelBased", display=FALSE, interval="profile", profileProgressInfo=FALSE)
    b <- bmd(ft, 0.1, def=def)
    expect_equal(unname(coef(b)), unname(rd$Results[1]), tolerance=1e-5)
    expect_equal(sqrt(vcov(b)[1, 1]), unname(rd$SE[1]), tolerance=1e-4)
    expect_equal(confint(b, level=0.9, progress=FALSE)[1, 1], unname(rp$Results[2]), tolerance=1e-4)
  }
  # continuous relative deviation; bmd profiles with the chi-square cutoff, not the small-sample one
  data(ryegrass, package="drcData", envir=environment())
  fd <- drc::drm(rootl ~ conc, data=ryegrass, fct=LL.4())
  ft <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  rd <- bmd::bmd(fd, 0.1, def="relative", backgType="modelBased", display=FALSE)
  rp <- bmd::bmd(fd, 0.1, def="relative", backgType="modelBased", display=FALSE, interval="profile",
                 profileProgressInfo=FALSE, profileGridSize=200)
  b <- bmd(ft, 0.1, def="relative")
  expect_equal(unname(coef(b)), unname(rd$Results[1]), tolerance=1e-5)
  expect_equal(sqrt(vcov(b)[1, 1]), unname(rd$SE[1]), tolerance=1e-4)
  expect_equal(suppressMessages(confint(b, level=0.9, method="wald"))[1, 1], unname(rd$Results[2]), tolerance=1e-5)
  ft$rdf <- NULL
  expect_equal(confint(bmd(ft, 0.1, def="relative"), level=0.9, progress=FALSE)[1, 1], unname(rp$Results[2]), tolerance=1e-4)
})

test_that("BMD in standard deviations includes the uncertainty of sigma", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  b <- bmd(f, 1, def="sd")
  mu0 <- coef(f)[["b2:(Intercept)"]]
  expect_equal(unname(coef(b)), unname(coef(ed(f, mu0 - residpar(f), type="absolute"))), tolerance=1e-8)
  # delta method by numerical derivatives of the BMD w.r.t. all parameters
  pf <- f$sdrl$par.fixed
  val <- function(p){
    pv <- f$start
    pv[!f$fix] <- p[seq_len(sum(!f$fix))]
    bmdvalue(f, b$bmd[[1]], pv, numeric(0), p[["log_sigma"]])
  }
  g <- numDeriv::grad(val, pf)
  expect_equal(sqrt(vcov(b)[1, 1]), sqrt(drop(t(g) %*% fixedcov(f) %*% g)), tolerance=1e-6)
  # profile limits: the deviance of the fit constrained to the BMD reaches the cutoff
  ci <- confint(b, level=0.9, progress=FALSE)
  pen <- function(p, x){
    v <- val(p)
    if (!is.finite(v)) return(1e10)
    f$obj$fn(p) + 1e7*(log(v) - log(x))^2
  }
  o <- optim(f$opt$par, pen, x=ci[1, 1], method="BFGS", control=list(maxit=5000, reltol=1e-14))
  o <- optim(o$par, pen, x=ci[1, 1], control=list(maxit=20000, reltol=1e-14))
  expect_equal(2*(f$obj$fn(o$par) - f$opt$objective), profcutoff(0.9, f), tolerance=1e-3)
})

test_that("marginal BMDs average the mean curve over the random effects", {
  d <- simbmdgroups()
  # random effects on the asymptote enter linearly: marginal and conditional BMDs are the same
  f2 <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=d, model="loglogistic")
  bm <- bmd(f2, 0.1, def="relative")
  expect_false(bm$integrated[[1]])
  expect_equal(coef(bm), coef(bmd(f2, 0.1, def="relative", scale="conditional")))

  f <- fit4(y ~ dose, rform=b2 + b4 ~ 1 | grp, data=d, model="loglogistic", link=c(b4="log"))
  bm <- bmd(f, c(0.05, 0.2), def="relative")
  bc <- bmd(f, c(0.05, 0.2), def="relative", scale="conditional")
  expect_true(all(bm$integrated))
  expect_true(all(coef(bm) < coef(bc)))
  # Monte Carlo integration of the mean curve
  cf <- coef(f)
  set.seed(1)
  U <- MASS::mvrnorm(2e5, c(0, 0), recov(f)[[1]]$cov)
  mm <- function(x) mean(cf[["b1:(Intercept)"]] + (cf[["b2:(Intercept)"]] + U[, 1] - cf[["b1:(Intercept)"]]) /
                           (1 + exp(-cf[["b3:(Intercept)"]]*(log(x) - cf[["b4:(Intercept)"]] - U[, 2]))))
  mc <- sapply(c(0.05, 0.2), function(r) uniroot(function(x) mm(x) - cf[["b2:(Intercept)"]]*(1 - r), c(1e-3, 50), tol=1e-10)$root)
  expect_equal(unname(coef(bm)), mc, tolerance=3e-3)
  # standard errors include the variance parameters of the random effects
  pf <- f$sdrl$par.fixed
  val <- function(p){
    pv <- f$start
    pv[!f$fix] <- p[seq_len(sum(!f$fix))]
    bmdvalue(f, bm$bmd[[1]], pv, unname(p[names(p) == "theta"]), p[["log_sigma"]])
  }
  g <- numDeriv::grad(val, pf)
  expect_equal(sqrt(vcov(bm)[1, 1]), sqrt(drop(t(g) %*% fixedcov(f) %*% g)), tolerance=1e-6)
  # more quadrature nodes do not change the estimate
  expect_equal(coef(bmd(f, c(0.05, 0.2), def="relative", nodes=30)), coef(bm), tolerance=1e-6)
  # profile intervals are not available: Wald intervals; bootstrap intervals
  expect_message(ci <- confint(bm, level=0.9), "Wald intervals")
  set.seed(1)
  cib <- suppressMessages(confint(bm, parm=1, level=0.9, method="boot", nsim=20, progress=FALSE))
  expect_true(all(is.finite(cib)) && cib[1] < coef(bm)[1] && cib[2] > coef(bm)[1])
  # BMDs of fitted groups
  bg <- bmd(f, 0.1, def="relative", newdata=data.frame(grp=c("1", "2")))
  expect_true(all(bg$conditional))
  expect_false(any(bg$integrated))
})

test_that("bmd checks its arguments", {
  d <- simgauss("loglogistic")
  f <- fit4(y ~ dose, data=d, model="loglogistic")
  expect_error(bmd(f, 0.1), "proportions")
  expect_error(bmd(f, -1, def="relative"), "positive")
  dg <- simbmdgroups()
  fr <- fit4(y ~ dose, rform=b2 ~ 1 | grp, data=dg, model="loglogistic")
  expect_error(bmd(fr, 1, def="sd"), "random effects")
  fq <- drmTMB(cbind(y, n - y) ~ dose, data=simquantal(), family="binomial", model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  expect_error(bmd(fq, 1, def="additional"), "below 1")
})
