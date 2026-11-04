# Hormesis models braincousens and cedergreen. drc's BC.5/BC.4 and CRS.5/CRS.4 map to b1 = c,
# b2 = d, b3 = -b, b4 = e, b5 = f. The fits start at the drc estimates, as the likelihood of these
# models is flat in places and the optimisers stop at slightly different points otherwise. drc does
# not reach the maximum likelihood for CRS.4b, CRS.4c and CRS.5a on lettuce; those are not compared.

skip_if_not_installed("drc")
skip_if_not_installed("drcData")
data(lettuce, package="drcData")

drccases <- list(
  BC.5   = list(fct=function() drc::BC.5(),   model="braincousens", alpha=NULL, b1fix=FALSE),
  BC.4   = list(fct=function() drc::BC.4(),   model="braincousens", alpha=NULL, b1fix=TRUE),
  CRS.4a = list(fct=function() drc::CRS.4a(), model="cedergreen",   alpha=1,    b1fix=TRUE)
)

for (nm in names(drccases)){
  test_that(paste("hormesis model", nm, "agrees with drc"), {
    cs <- drccases[[nm]]
    f2 <- drc::drm(weight ~ conc, data=lettuce, fct=cs$fct(), control=drc::drmc(relTol=1e-12))
    p <- coef(f2)
    names(p) <- sub(":.*", "", names(p))
    st <- c(b1=if (cs$b1fix) 0 else p[["c"]], b2=p[["d"]], b3=-p[["b"]], b4=p[["e"]], b5=p[["f"]])
    f1 <- drmTMB(weight ~ conc, data=lettuce, model=cs$model, alpha=cs$alpha, start=st, fix=c(b1=cs$b1fix))
    free <- if (cs$b1fix) 2:5 else 1:5

    expect_equal(unname(bcoef(f1)), unname(st[free]), tolerance=1e-3)
    expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-6)
    # drc's vcov is in the order b, c, d, e, f, unnamed
    ord <- match(c("c", "d", "b", "e", "f"), c("b", if (!cs$b1fix) "c", "d", "e", "f"))
    expect_equal(unname(sqrt(diag(vcov(f1)))[seq_along(free)]), unname(sqrt(diag(vcov(f2))))[ord[!is.na(ord)]], tolerance=1e-2)

    e1 <- ed(f1, respLev=c(10, 50))
    e2 <- drc::ED(f2, c(10, 50), display=FALSE, upper=1e6)
    expect_equal(unname(coef(e1)), unname(e2[, 1]), tolerance=2e-3)
    expect_equal(unname(sqrt(diag(vcov(e1)))), unname(e2[, 2]), tolerance=1e-2)
  })
}

test_that("anova tests for hormesis as drc", {
  f0 <- drmTMB(weight ~ conc, data=lettuce, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))
  f1 <- drmTMB(weight ~ conc, data=lettuce, model="braincousens")
  a <- anova(f0, f1, test="F")
  d0 <- drc::drm(weight ~ conc, data=lettuce, fct=drc::LL.4())
  d1 <- drc::drm(weight ~ conc, data=lettuce, fct=drc::BC.5())
  ad <- anova(d0, d1, details=FALSE)
  expect_equal(a$F[2], ad[2, "F value"], tolerance=1e-3)
  expect_equal(a[["Pr(>F)"]][2], ad[2, "p value"], tolerance=1e-3)
})

test_that("the mean functions of the template and of R agree for the hormesis models", {
  for (m in list(list("braincousens", NULL), list("cedergreen", 0.25))){
    # NA/NaN evaluations of the objective on the way to the estimates
    f <- suppressWarnings(drmTMB(weight ~ conc, data=lettuce, model=m[[1]], alpha=m[[2]]))
    b <- coef(f)
    nd <- data.frame(conc=c(0, 0.1, 1, 5, 50))
    mu <- drcTMB:::drcmeanR(drcTMB:::modelindex(m[[1]]), nd$conc, b[[1]], b[[2]], b[[3]], b[[4]], b[[5]], m[[2]])
    p <- predict(f, newdata=nd, se.fit=TRUE)
    expect_equal(p$fit, mu)
    expect_equal(predict(f), as.vector(fitted(f)))
    # dose zero: the response of the control
    expect_equal(mu[1], b[[2]])
  }
})

test_that("effective doses of the hormesis models have the delta method standard errors", {
  skip_if_not_installed("numDeriv")
  for (m in list(list("braincousens", NULL), list("cedergreen", 1))){
    f <- drmTMB(weight ~ conc, data=lettuce, model=m[[1]], alpha=m[[2]])
    mod <- drcTMB:::modelindex(m[[1]])
    b <- coef(f)[1:5]
    e <- ed(f, c(10, 50))
    edR <- function(b, q) drcTMB:::edfun(mod, q, b[3], b[4], b[5], b[5]/(b[2] - b[1]), m[[2]])
    g <- t(sapply(1:2, function(i) numDeriv::grad(edR, b, q=e$q[[i]])))
    expect_equal(unname(coef(e)), sapply(1:2, function(i) edR(b, e$q[[i]])))
    expect_equal(unname(vcov(e)), unname(g %*% vcov(f)[1:5, 1:5] %*% t(g)), tolerance=1e-5)

    # the effective dose is above the peak: the response there is the level of the effective dose,
    # and the curve is decreasing
    x <- coef(e)[[2]]
    mu <- function(x) drcTMB:::drcmeanR(mod, x, b[[1]], b[[2]], b[[3]], b[[4]], b[[5]], m[[2]])
    expect_equal(mu(x), (b[[1]] + b[[2]])/2)
    expect_lt(mu(1.01*x), mu(x))

    # absolute effective doses at the same response
    ea <- ed(f, mu(x), type="absolute")
    expect_equal(coef(ea)[[1]], x, tolerance=1e-8)

    # b4 from the effective dose, as in the profiles
    q <- e$q[[2]]
    expect_equal(drcTMB:::edb4(mod, q, x, b[[3]], b[[5]], b[[5]]/(b[[2]] - b[[1]]), m[[2]]), b[[4]])
  }
})

test_that("profile and bootstrap intervals of hormesis effective doses contain the estimates", {
  f <- drmTMB(weight ~ conc, data=lettuce, model="braincousens")
  e <- ed(f, 10)
  ci <- suppressMessages(confint(e, progress=FALSE))
  expect_true(ci[1, 1] < coef(e) && coef(e) < ci[1, 2])
  set.seed(1)
  cb <- suppressMessages(confint(e, method="boot", nsim=30, progress=FALSE))
  expect_true(cb[1, 1] < coef(e) && coef(e) < cb[1, 2])
})

test_that("hormesis mixed models recover the simulated curve", {
  set.seed(3)
  d <- simgroups(dose=c(0, 0.1, 0.3, 1, 3, 10, 30), ng=12)
  u2 <- rnorm(12, 0, 0.5)
  u4 <- rnorm(12, 0, 0.3)
  d$y <- drcTMB:::drcmeanR(11, d$dose, 1, 8 + u2[d$grp], -2, 3*exp(u4[d$grp]), 2) + rnorm(nrow(d), 0, 0.4)
  f <- drmTMB(y ~ dose, rform=b2 + b4 ~ 1 | grp, data=d, model="braincousens", link=c(b4="log"))
  expect_equal(f$opt$convergence, 0)
  cf <- coef(f)
  expect_equal(unname(cf[1:5]), c(1, 8, -2, log(3), 2), tolerance=0.15)

  # effective doses of fitted groups: the curves with their predicted random effects
  e <- ed(f, 50, newdata=data.frame(grp=c("1", "2")))
  re <- predict(f, newdata=data.frame(dose=coef(e), grp=c("1", "2")))
  expect_equal(re, (cf[[1]] + predict(f, newdata=data.frame(dose=0, grp=c("1", "2"))))/2, tolerance=1e-6)
})

test_that("alpha is checked", {
  expect_error(drmTMB(weight ~ conc, data=lettuce, model="cedergreen"), "alpha")
  expect_error(drmTMB(weight ~ conc, data=lettuce, model="cedergreen", alpha=-1), "positive")
  expect_error(drmTMB(weight ~ conc, data=lettuce, model="braincousens", alpha=1), "only used by the cedergreen")
  f <- suppressWarnings(drmTMB(weight ~ conc, data=lettuce, model="cedergreen", alpha=0.25))
  expect_output(print(f), "cedergreen \\(alpha = 0.25\\)")
})

# drcTMB fit of a hormesis model of drc, starting at the drc estimates
drcstart <- function(f2, b1fix){
  p <- coef(f2)
  names(p) <- sub(":.*", "", names(p))
  c(b1=if (b1fix) 0 else p[["c"]], b2=p[["d"]], b3=-p[["b"]], b4=p[["e"]], b5=p[["f"]])
}

test_that("the peak agrees with MAX of drc", {
  f2 <- drc::drm(weight ~ conc, data=lettuce, fct=drc::BC.5(), control=drc::drmc(relTol=1e-12))
  f1 <- drmTMB(weight ~ conc, data=lettuce, model="braincousens", start=drcstart(f2, FALSE))
  pk <- peak(f1)
  expect_equal(names(coef(pk)), c("dose:lf1", "response:lf1"))
  expect_equal(unname(coef(pk)), unname(drc::MAX(f2)[1, ]), tolerance=1e-3)

  # CRS.4a on lettuce: the stimulation does not lift the curve above the control, no peak
  f2 <- drc::drm(weight ~ conc, data=lettuce, fct=drc::CRS.4a(), control=drc::drmc(relTol=1e-12))
  f1 <- drmTMB(weight ~ conc, data=lettuce, model="cedergreen", alpha=1, start=drcstart(f2, TRUE), fix=c(b1=TRUE))
  expect_warning(pk <- peak(f1), "no peak")
  expect_true(all(is.na(coef(pk))))
})

test_that("the peak and the return to control have the delta method standard errors", {
  skip_if_not_installed("numDeriv")
  for (m in list(list("braincousens", NULL), list("cedergreen", 0.25))){
    f <- suppressWarnings(drmTMB(weight ~ conc, data=lettuce, model=m[[1]], alpha=m[[2]]))
    mod <- drcTMB:::modelindex(m[[1]])
    b <- coef(f)[1:5]
    mu <- function(x, b) drcTMB:::drcmeanR(mod, x, b[[1]], b[[2]], b[[3]], b[[4]], b[[5]], m[[2]])
    pk <- peak(f)
    # peak as the root of the analytic derivative of the curve on the log-dose scale
    val <- function(b, j){
      D <- function(s) drcTMB:::hormderivs(f, b, exp(s))[1]
      s <- uniroot(D, log(coef(pk)[[1]]) + c(-0.5, 0.5), tol=1e-14)$root
      if (j == 1) exp(s) else mu(exp(s), b)
    }
    g <- rbind(numDeriv::grad(val, b, j=1), numDeriv::grad(val, b, j=2))
    expect_equal(unname(coef(pk)), c(val(b, 1), val(b, 2)), tolerance=1e-8)
    expect_equal(unname(vcov(pk)), unname(g %*% vcov(f)[1:5, 1:5] %*% t(g)), tolerance=1e-5)
    expect_gt(coef(pk)[[2]], b[[2]])
    expect_equal(drcTMB:::peakb4(mod, coef(pk)[[1]], b[[3]], b[[5]]/(b[[2]] - b[[1]]), m[[2]]), b[[4]])

    # return to control: the dose above the peak with the response of the control
    e0 <- ed(f, 0)
    x0 <- coef(e0)[[1]]
    expect_gt(x0, coef(pk)[[1]])
    expect_equal(mu(x0, b), b[[2]])
    edR <- function(b) drcTMB:::edfun(mod, 1, b[3], b[4], b[5], b[5]/(b[2] - b[1]), m[[2]])
    g0 <- numDeriv::grad(edR, b)
    expect_equal(vcov(e0)[1, 1], as.numeric(g0 %*% vcov(f)[1:5, 1:5] %*% g0), tolerance=1e-5)
  }
})

test_that("intervals of the peak", {
  f <- drmTMB(weight ~ conc, data=lettuce, model="braincousens")
  pk <- peak(f)
  expect_message(ci <- confint(pk, progress=FALSE), "Wald intervals for the response")
  expect_true(all(ci[, 1] < coef(pk) & coef(pk) < ci[, 2]))
  expect_equal(ci[2, ], confint(pk, parm=2, method="wald")[1, ])
  set.seed(1)
  cb <- suppressMessages(confint(pk, method="boot", nsim=30, progress=FALSE))
  expect_true(all(cb[, 1] < coef(pk) & coef(pk) < cb[, 2]))
  expect_error(peak(drmTMB(weight ~ conc, data=lettuce, model="loglogistic", start=c(b5=0), fix=c(b5=TRUE))), "hormesis and substrate inhibition models")
})

test_that("peaks of fitted groups are the maxima of their predicted curves", {
  set.seed(3)
  d <- simgroups(dose=c(0, 0.1, 0.3, 1, 3, 10, 30), ng=12)
  u2 <- rnorm(12, 0, 0.5)
  u4 <- rnorm(12, 0, 0.3)
  d$y <- drcTMB:::drcmeanR(11, d$dose, 1, 8 + u2[d$grp], -2, 3*exp(u4[d$grp]), 2) + rnorm(nrow(d), 0, 0.4)
  f <- drmTMB(y ~ dose, rform=b2 + b4 ~ 1 | grp, data=d, model="braincousens", link=c(b4="log"))
  pk <- peak(f, newdata=data.frame(grp=c("1", "2")))
  for (g in 1:2){
    o <- optimise(function(x) predict(f, newdata=data.frame(dose=x, grp=as.character(g))), c(0.1, 3), maximum=TRUE, tol=1e-10)
    expect_equal(coef(pk)[[paste("dose:grp=", g, sep="")]], o$maximum, tolerance=1e-6)
    expect_equal(coef(pk)[[paste("response:grp=", g, sep="")]], o$objective, tolerance=1e-8)
  }
})
