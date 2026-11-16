# BMD parameterisation (bmdpar): b4 replaced by the log BMD of the curve of the intercepts

library(drcData)
data(ryegrass)

ctl <- list(smallsample=FALSE)

# the fit in the BMD parameterisation has the likelihood of the standard fit, and its BMD and the
# standard error of the log BMD are those of bmd() (conditional, maximum likelihood)
expect_bmdpar_equal <- function(args, bmr, def, tol=1e-5){
  f0 <- do.call(drmTMB, c(args, list(control=ctl)))
  fb <- do.call(drmTMB, c(args, list(control=ctl, bmdpar=list(bmr=bmr, def=def))))
  b0 <- bmd(f0, bmr, def, scale="conditional")
  expect_equal(as.numeric(logLik(fb)), as.numeric(logLik(f0)), tolerance=1e-8)
  expect_equal(exp(coef(fb)[["log_bmd"]]), unname(coef(b0)), tolerance=tol)
  expect_equal(sqrt(vcov(fb)["log_bmd", "log_bmd"]), sqrt(vcov(b0)[1, 1])/unname(coef(b0)), tolerance=1e-3)
  invisible(list(f0=f0, fb=fb))
}

test_that("the BMD parameterisation reproduces the standard fit and bmd()", {
  base <- list(form=rootl ~ conc, data=ryegrass)
  r <- expect_bmdpar_equal(c(base, list(model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))), 0.05, "relative")
  expect_equal(names(coef(r$fb))[3], "log_bmd")
  expect_equal(r$fb$bmdpar$inc, -1)
  expect_null(r$fb$rdf)
  expect_output(print(r$fb), "log benchmark dose")
  expect_output(print(summary(r$fb)), "log benchmark dose")
  expect_bmdpar_equal(c(base, list(model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))), 1, "sd")
  expect_bmdpar_equal(c(base, list(model="weibull2", start=c(b1=0), fix=c(b1=TRUE))), 0.1, "relative")
  expect_bmdpar_equal(c(base, list(model="gamma", start=c(b2=0), fix=c(b2=TRUE))), 0.1, "relative")
  expect_bmdpar_equal(c(base, list(model="nec", start=c(b2=0, b5=0), fix=c(b2=TRUE, b5=TRUE))), 0.1, "relative")
  expect_bmdpar_equal(c(base, list(model="exponential")), 0.1, "relative")
  expect_bmdpar_equal(c(base, list(model="power")), 0.1, "relative")
  expect_bmdpar_equal(c(base, list(model="loglogistic", family="lognormal", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE))), 1, "sd")
})

test_that("the BMD parameterisation for quantal responses", {
  set.seed(1)
  dq <- data.frame(dose=rep(c(0, 1, 2, 4, 8, 16), each=4), n=50)
  p <- 0.05 + 0.95/(1 + exp(-1.5*(log(pmax(dq$dose, 1e-8)) - log(5))))
  dq$y <- rbinom(nrow(dq), dq$n, p)
  q <- list(form=cbind(y, n - y) ~ dose, data=dq, family="binomial")
  expect_bmdpar_equal(c(q, list(model="loglogistic", parlink=c(b1="logit"), start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE))), 0.1, "excess")
  expect_bmdpar_equal(c(q, list(model="weibull1", start=c(b2=1), fix=c(b2=TRUE))), 0.05, "additional")
  # multistage: the computed coefficient b4 keeps its lower bound 0
  r <- expect_bmdpar_equal(c(q, list(model="multistage", start=c(b2=1), fix=c(b2=TRUE))), 0.1, "excess")
  expect_gte(r$fb$obj$report(r$fb$opt$par)$f4bmd, 0)
})

test_that("the BMD parameterisation with random effects, and a prior for the BMD", {
  data(vinclozolin)
  r <- expect_bmdpar_equal(list(form=effect ~ conc, rform=b1 + b2 ~ 1 | exper, data=vinclozolin, model="loglogistic",
                                start=c(b5=0), fix=c(b5=TRUE)), 0.1, "relative")
  expect_equal(dim(simulate(r$fb, nsim=2)), c(nrow(vinclozolin), 2))

  fp <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
               bmdpar=list(bmr=0.05, def="relative"), prior=list(bmd=pert(0.1, 0.5, 3)))
  f0 <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
               bmdpar=list(bmr=0.05, def="relative"))
  # the prior is the density of the BMD: the log Jacobian of the log BMD is added
  lb <- coef(fp)[["log_bmd"]]
  expect_equal(fp$nlprior, -(lpdfprior(pert(0.1, 0.5, 3), exp(lb)) + lb), tolerance=1e-10)
  expect_equal(fp$opt$objective, f0$obj$fn(fp$opt$par) + fp$nlprior, tolerance=1e-10)
  expect_equal(unname(fp$lower[3]), log(0.1))
  expect_lt(exp(lb), exp(coef(f0)[["log_bmd"]]))
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
                      bmdpar=list(bmr=0.05, def="relative"), prior=list(b4=normal(3, 1))), "prior for the BMD")
})

test_that("methods and models without the BMD parameterisation", {
  fb <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", start=c(b1=0, b5=0), fix=c(b1=TRUE, b5=TRUE),
               bmdpar=list(bmr=0.05, def="relative"))
  expect_error(ed(fb, 50), "BMD parameterisation")
  expect_error(bmd(fb, 0.1), "BMD parameterisation")
  expect_error(predict(fb), "BMD parameterisation")
  expect_true(all(is.finite(suppressMessages(confint(fb)))))
  bp <- list(bmr=0.1, def="relative")
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="logistic", bmdpar=bp), "logistic model")
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="braincousens", bmdpar=bp), "braincousens model")
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", rform=b4 ~ 1 | conc, bmdpar=bp), "random effects on b4")
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", bmdpar=list(bmr=0.1)), "list with")
  expect_error(drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic", bmdpar=list(bmr=0.1, def="excess")), "proportions")
  data(S.alba)
  expect_error(drmTMB(DryMatter ~ Dose, fform=b2 ~ Herbicide, data=S.alba, model="loglogistic", bmdpar=bp), "without covariates")
})
