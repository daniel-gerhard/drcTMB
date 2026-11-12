# models without an upper asymptote, as in PROAST and BMDS: exponential f(x) = b1 exp(b4 x^b3) and
# power f(x) = b1 + b4 x^b3

simexppower <- function(seed=1){
  set.seed(seed)
  d <- data.frame(dose=rep(c(0, 1, 2, 5, 10, 20), each=5))
  d$ye <- 10*exp(0.02*d$dose^1.3) + rnorm(nrow(d), 0, 0.5)
  d$yp <- 5 - 0.1*d$dose^1.5 + rnorm(nrow(d), 0, 0.5)
  d$cnt <- rpois(nrow(d), 3*exp(0.08*d$dose))
  d
}

test_that("exponential and power agree with nls", {
  d <- simexppower()
  fe <- drmTMB(ye ~ dose, data=d, model="exponential")
  ne <- nls(ye ~ a*exp(b*dose^c), data=d, start=list(a=10, b=0.02, c=1.3))
  expect_equal(unname(coef(fe)[1:3]), unname(coef(ne)[c("a", "c", "b")]), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fe)), as.numeric(logLik(ne)), tolerance=1e-6)
  fp <- drmTMB(yp ~ dose, data=d, model="power")
  np <- nls(yp ~ a + b*dose^c, data=d, start=list(a=5, b=-0.1, c=1.5))
  expect_equal(unname(coef(fp)[1:3]), unname(coef(np)[c("a", "c", "b")]), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fp)), as.numeric(logLik(np)), tolerance=1e-6)
  # nls uses the Gauss-Newton approximation of the Hessian
  expect_equal(unname(sqrt(diag(vcov(fp)))[1:3]), unname(sqrt(diag(vcov(np)))[c(1, 3, 2)]), tolerance=0.02)
  # the power b3 is bounded below by 0; b2 and b5 are not used
  expect_equal(fp$lower[2], 0)
  expect_equal(sum(!fp$fix), 3)
  expect_error(drmTMB(yp ~ dose, rform=b2 ~ 1 | dose, data=d, model="power"), "does not use b2")
  expect_error(drmTMB(yp ~ dose, data=d, model="power", link=c(b4="log")), "of either sign")
  # mean of the template and of R, at dose zero b1
  nd <- data.frame(dose=c(0, 1, 7, 30))
  expect_equal(predict(fe, nd, se.fit=TRUE)$fit, predict(fe, nd))
  expect_equal(predict(fe, nd), unname(predict(ne, nd)), tolerance=1e-5)
  expect_equal(predict(fp, nd), unname(predict(np, nd)), tolerance=1e-5)
  expect_true(is.finite(predict(fp, nd, se.fit=TRUE)$se.fit[1]))
  expect_error(ed(fe, 50), "no upper asymptote")
})

test_that("exponential with b3 = 1 (exp2) is a log-linear model, and anova tests the power", {
  d <- simexppower()
  f2 <- drmTMB(cnt ~ dose, data=d, model="exponential", family="poisson", start=c(b3=1), fix=c(b3=TRUE), link=c(b1="log"))
  g <- glm(cnt ~ dose, family=poisson, data=d)
  expect_equal(unname(coef(f2)), unname(coef(g)), tolerance=1e-6)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(g)), tolerance=1e-8)
  f3 <- drmTMB(cnt ~ dose, data=d, model="exponential", family="poisson", link=c(b1="log"))
  a <- anova(f2, f3)
  expect_equal(a$Df[2], 1)
  # lognormal: linear in log(y)
  set.seed(2)
  d$yl <- 10*exp(0.02*d$dose^1.3)*exp(rnorm(nrow(d), 0, 0.1))
  fl <- drmTMB(yl ~ dose, data=d, model="exponential", family="lognormal", link=c(b1="log"))
  nl <- nls(log(yl) ~ la + b*dose^c, data=d, start=list(la=2.3, b=0.02, c=1.3))
  expect_equal(unname(coef(fl)[1:3]), unname(coef(nl)[c("la", "c", "b")]), tolerance=1e-5)
})

test_that("absolute effective doses and BMDs of exponential and power in closed form", {
  d <- simexppower()
  fe <- drmTMB(ye ~ dose, data=d, model="exponential")
  b <- unname(coef(fe)[1:3])
  e <- ed(fe, c(12, 15), type="absolute")
  expect_equal(unname(coef(e)), (log(c(12, 15)/b[1])/b[3])^(1/b[2]))
  # responses below the response at dose zero are not reached by an increasing curve
  expect_warning(e0 <- ed(fe, 9, type="absolute"), "do not reach")
  expect_true(is.na(coef(e0)))
  bm <- bmd(fe, c(0.05, 0.1), def="relative")
  expect_equal(unname(coef(bm)), (log(1 + c(0.05, 0.1))/b[3])^(1/b[2]))
  # delta method standard errors
  g <- numDeriv::jacobian(function(p) (log(1.1)/p[3])^(1/p[2]), b)
  expect_equal(vcov(bm)[2, 2], (g %*% vcov(fe)[1:3, 1:3] %*% t(g))[1, 1], tolerance=1e-6)
  # profile BMDLs, close to the Wald limits
  ci <- suppressMessages(confint(bm, level=0.9))
  cw <- confint(bm, level=0.9, method="wald")
  expect_true(all(ci[, 1] < coef(bm) & ci[, 2] > coef(bm)))
  expect_equal(ci, cw, tolerance=0.05)

  # decreasing power curve: the relative BMD is a decrease, and BMDs in standard deviations
  fp <- drmTMB(yp ~ dose, data=d, model="power")
  b <- unname(coef(fp)[1:3])
  expect_lt(b[3], 0)
  bp <- bmd(fp, 0.1, def="relative")
  expect_equal(unname(coef(bp)), (-0.1*b[1]/b[3])^(1/b[2]))
  sig <- residpar(fp)
  bs <- bmd(fp, 1, def="sd")
  expect_equal(unname(coef(bs)), (-sig/b[3])^(1/b[2]), tolerance=1e-6)
  cs <- suppressMessages(confint(bs, level=0.9))
  expect_true(cs[1] < coef(bs) && cs[2] > coef(bs))
  ea <- ed(fp, 3, type="absolute")
  expect_equal(unname(coef(ea)), ((3 - b[1])/b[3])^(1/b[2]))
  ca <- suppressMessages(confint(ea))
  expect_equal(ca, confint(ea, method="wald"), tolerance=0.05)
})

test_that("marginal BMDs of the exponential model with a random coefficient b4", {
  set.seed(3)
  d <- expand.grid(dose=c(0, 1, 2, 5, 10, 20), rep=1:2, grp=factor(1:10))
  u <- rnorm(10, 0, 0.01)
  d$y <- 10*exp((0.03 + u[d$grp])*d$dose^1.2)*exp(rnorm(nrow(d), 0, 0.05))
  f <- drmTMB(y ~ dose, rform=b4 ~ 1 | grp, data=d, model="exponential", family="lognormal")
  bm <- bmd(f, 0.1, def="relative")
  bc <- bmd(f, 0.1, def="relative", scale="conditional")
  expect_true(bm$integrated[[1]])
  cf <- coef(f)
  b <- cf[["b1:(Intercept)"]]
  expect_equal(unname(coef(bc)), (log(1.1)/cf[["b4:(Intercept)"]])^(1/cf[["b3:(Intercept)"]]))
  # Monte Carlo integration of the median curve
  set.seed(1)
  U <- rnorm(2e5, 0, recov(f)[[1]]$sd)
  mm <- function(x) mean(b*exp((cf[["b4:(Intercept)"]] + U)*x^cf[["b3:(Intercept)"]]))
  mc <- uniroot(function(x) mm(x) - 1.1*b, c(1e-3, 50), tol=1e-10)$root
  expect_equal(unname(coef(bm)), mc, tolerance=3e-3)
  # standard errors include the variance of the random effects
  pf <- f$sdrl$par.fixed
  val <- function(p){
    pv <- f$start
    pv[!f$fix] <- p[seq_len(sum(!f$fix))]
    bmdvalue(f, bm$bmd[[1]], pv, unname(p[names(p) == "theta"]), p[["log_sigma"]])
  }
  g <- numDeriv::grad(val, pf)
  expect_equal(sqrt(vcov(bm)[1, 1]), sqrt(drop(t(g) %*% fixedcov(f) %*% g)), tolerance=1e-5)
  # BMDs of fitted groups
  bg <- bmd(f, 0.1, def="relative", newdata=data.frame(grp=c("1", "2")))
  pg <- predict(f, data.frame(dose=1, grp=c("1", "2")), type="parameters")
  expect_equal(unname(coef(bg)), (log(1.1)/pg[, "b4"])^(1/pg[, "b3"]), tolerance=1e-6)
})

test_that("model averaging with the exponential and power models", {
  d <- simexppower()
  fe <- drmTMB(ye ~ dose, data=d, model="exponential")
  fp <- drmTMB(ye ~ dose, data=d, model="power")
  ma <- bmdma(fe, fp, bmr=0.1, def="relative", average="bmd", nsim=20, seed=1, progress=FALSE)
  w <- ma$models$weight
  expect_equal(unname(coef(ma)), sum(w*c(coef(bmd(fe, 0.1, def="relative")), coef(bmd(fp, 0.1, def="relative")))))
  mc <- bmdma(fe, fp, bmr=0.1, def="relative", nsim=20, seed=1, progress=FALSE)
  expect_true(is.finite(coef(mc)))
})
