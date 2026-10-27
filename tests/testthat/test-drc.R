# Models without random effects maximise the same likelihood as drc::drm; estimates, ED and
# log-likelihoods agree up to the convergence tolerance of the optimisers.
#
# drc parameterisation, mapped to b1-b5:
#   LL.4, L.4, W1.4, G.4: b1 = c, b2 = d, b3 = -b, b4 = e
#   W2.4, LN.4:           b1 = c, b2 = d, b3 =  b, b4 = e
#   LL.5, L.5:            as LL.4, b5 = log(f)
#   MM.3, AR.3:           b1 = c, b2 = d, b4 = e
#   EXD.3:                b1 = d, b2 = c, b4 = e
# drc estimates the residual variance by RSS / (n - p), drcTMB by maximum likelihood, so the
# gaussian standard errors of drc are larger by sqrt(n / (n - p)).

skip_if_not_installed("drc")
skip_if_not_installed("numDeriv")

drcmap <- list(
  loglogistic          = list(fct=function() drc::LL.4(), map=function(p) c(p["c"], p["d"], -p["b"], p["e"])),
  logistic             = list(fct=function() drc::L.4(),  map=function(p) c(p["c"], p["d"], -p["b"], p["e"])),
  weibull1             = list(fct=function() drc::W1.4(), map=function(p) c(p["c"], p["d"], -p["b"], p["e"])),
  weibull2             = list(fct=function() drc::W2.4(), map=function(p) c(p["c"], p["d"],  p["b"], p["e"])),
  lognormal            = list(fct=function() drc::LN.4(), map=function(p) c(p["c"], p["d"],  p["b"], p["e"])),
  gompertz             = list(fct=function() drc::G.4(),  map=function(p) c(p["c"], p["d"], -p["b"], p["e"])),
  michaelismenten      = list(fct=function() drc::MM.3(), map=function(p) c(p["c"], p["d"], p["e"])),
  asymptoticregression = list(fct=function() drc::AR.3(), map=function(p) c(p["c"], p["d"], p["e"])),
  exponentialdecay     = list(fct=function() drc::EXD.3(), map=function(p) c(p["d"], p["c"], p["e"]))
)

# drc coefficients named b, c, d, e, f
drccoef <- function(fit){
  p <- coef(fit)
  names(p) <- sub(":.*", "", names(p))
  p
}

# drc covariance matrix in the order of the drcTMB coefficients
drcvcov <- function(fit, map){
  p <- drccoef(fit)
  ord <- match(names(map(p)), names(p))
  vcov(fit)[ord, ord]
}

for (m in names(drcmap)){
  test_that(paste("gaussian", m, "agrees with drc"), {
    d <- simgauss(m)
    f1 <- fit4(y ~ dose, data=d, model=m)
    f2 <- drc::drm(y ~ dose, data=d, fct=drcmap[[m]]$fct())
    map <- drcmap[[m]]$map

    expect_equal(unname(bcoef(f1)), unname(map(drccoef(f2))), tolerance=1e-3)
    expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-5)

    # small-sample inference as drc: unbiased residual variance with n - p degrees of freedom
    np <- length(bcoef(f1))
    se1 <- sqrt(diag(vcov(f1)))[seq_len(np)]
    se2 <- sqrt(diag(drcvcov(f2, map)))
    expect_equal(unname(se1), unname(se2), tolerance=1e-2)
    expect_equal(residpar(f1), summary(f2)$rseMat[1, 1], tolerance=1e-4)

    e1 <- ed(f1, respLev=c(10, 50, 90))
    e2 <- drc::ED(f2, c(10, 50, 90), display=FALSE)
    expect_equal(unname(e1$estimates), unname(e2[, 1]), tolerance=1e-3)
    # drc::ED gives too large standard errors for G.4; compare with the delta method instead
    if (m != "gompertz"){
      expect_equal(unname(sqrt(diag(e1$variance))), unname(e2[, 2]), tolerance=1e-2)
    } else {
      edg <- function(b, q) b[4] - log(-log(q))/b[3]
      se <- sapply(c(0.1, 0.5, 0.9), function(q){
        g <- numDeriv::grad(edg, bcoef(f1), q=q)
        sqrt(drop(t(g) %*% vcov(f1)[1:4, 1:4] %*% g))
      })
      expect_equal(unname(sqrt(diag(e1$variance))), se, tolerance=1e-4)
    }
  })
}

test_that("five-parameter loglogistic and logistic agree with drc (b5 = log f)", {
  d <- simgauss("loglogistic")
  f1 <- drmTMB(y ~ dose, data=d, model="loglogistic")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::LL.5())
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(c(p["c"], p["d"], -p["b"], p["e"], log(p["f"]))), tolerance=1e-3)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-5)
  expect_equal(unname(ed(f1, respLev=c(10, 50, 90))$estimates),
               unname(drc::ED(f2, c(10, 50, 90), display=FALSE)[, 1]), tolerance=1e-3)

  d <- simgauss("logistic")
  f1 <- drmTMB(y ~ dose, data=d, model="logistic")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::L.5())
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(c(p["c"], p["d"], -p["b"], p["e"], log(p["f"]))), tolerance=1e-3)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-5)
})

test_that("binomial loglogistic agrees with drc", {
  set.seed(4)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=4), n=50)
  d$y <- rbinom(nrow(d), d$n, drcmean("loglogistic", d$dose, c(0.15, 0.85, 2, 2.5, 0)))

  # four parameters; the binomial log-likelihood constants of drc differ, so no logLik comparison
  f1 <- fit4(cbind(y, n - y) ~ dose, data=d, model="loglogistic", family="binomial")
  f2 <- drc::drm(y/n ~ dose, weights=n, data=d, fct=drc::LL.4(), type="binomial")
  map <- drcmap$loglogistic$map
  expect_equal(unname(bcoef(f1)), unname(map(drccoef(f2))), tolerance=5e-3)
  expect_equal(unname(sqrt(diag(vcov(f1)))), unname(sqrt(diag(drcvcov(f2, map)))), tolerance=1e-2)
  expect_equal(unname(ed(f1, respLev=c(10, 50, 90))$estimates),
               unname(drc::ED(f2, c(10, 50, 90), display=FALSE)[, 1]), tolerance=5e-3)

  # two parameters: asymptotes fixed at 0 and 1
  d$y <- rbinom(nrow(d), d$n, drcmean("loglogistic", d$dose, c(0, 1, 2, 2.5, 0)))
  f1 <- drmTMB(cbind(y, n - y) ~ dose, data=d, model="loglogistic", family="binomial",
               start=c(b1=0, b2=1, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE))
  f2 <- drc::drm(y/n ~ dose, weights=n, data=d, fct=drc::LL.2(), type="binomial")
  p <- drccoef(f2)
  expect_equal(unname(bcoef(f1)), unname(c(-p["b"], p["e"])), tolerance=1e-5)
  expect_equal(f1$opt$objective, -as.numeric(logLik(f2)), tolerance=1e-6)
  expect_equal(unname(sqrt(diag(vcov(f1)))), unname(sqrt(diag(vcov(f2)))), tolerance=1e-4)
})

test_that("poisson loglogistic agrees with drc", {
  set.seed(4)
  d <- data.frame(dose=rep(c(0, 0.5, 1, 2, 4, 8, 16), each=4))
  d$y <- rpois(nrow(d), drcmean("loglogistic", d$dose, c(5, 45, -1.2, 3, 0)))
  f1 <- fit4(y ~ dose, data=d, model="loglogistic", family="poisson")
  f2 <- drc::drm(y ~ dose, data=d, fct=drc::LL.4(), type="Poisson")
  map <- drcmap$loglogistic$map
  # drc stops slightly short of the optimum: compare the likelihood at the drc estimates
  expect_lte(f1$opt$objective, f1$obj$fn(map(drccoef(f2))) + 1e-8)
  expect_equal(unname(bcoef(f1)), unname(map(drccoef(f2))), tolerance=1e-2)
  expect_equal(unname(sqrt(diag(vcov(f1)))), unname(sqrt(diag(drcvcov(f2, map)))), tolerance=1e-2)
})
