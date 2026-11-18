# drcTMB 0.1.0

First release.

## Models and data

* `drmTMB()` fits nonlinear dose-response models with the Template Model Builder: log-logistic,
  logistic, Weibull (two types), log-normal, Gompertz, Michaelis-Menten, asymptotic regression,
  gamma, multistage, the threshold model nec, the hormesis models of Brain-Cousens and Cedergreen,
  the substrate inhibition model of Haldane, the exponential and power models of PROAST and BMDS,
  and a constant model.
* Response distributions: gaussian, lognormal, binomial, beta-binomial, beta, poisson, negative
  binomial (nbinom1, nbinom2), and time-to-event data (`family = "event"`).
* Covariates on any parameter (`fform`), random effects on any parameter with the syntax of lme4
  (`rform`), parameter links (`parlink`) and links of the whole curve (`curvelink`).
* Correlated residuals (`cor_ar1()`, `cor_car1()`, `cor_cs()`, `cor_exp()`, `cor_gaus()`,
  `cor_spher()`) and variance functions (`var_power()`, `var_exp()`, `var_ident()`) for the gaussian
  and lognormal families.
* Priors on all parameters (`prior`, see `?priors`), with the estimates as posterior modes.
* Starting values by linearised models or a grid search (`drmTMBstart()`); choice of the optimiser
  and a restart at the estimates (`drmTMBcontrol()`).
* Data sets from the package drcData, which is not on CRAN: chickweed, daphnids, earthworms,
  germination, lettuce, liver.tumor, ryegrass, S.alba, selenium, TCDD and vinclozolin.

## Inference

* Wald, profile likelihood and parametric bootstrap intervals (`confint()`, `profile()`,
  `bootstrap()`); small-sample inference with t quantiles and F-based profiles for gaussian and
  lognormal models, as in drc.
* Likelihood ratio, F and parametric bootstrap tests with `anova()`.
* Robust (sandwich) covariance, also for clustered data (`vcovrobust()`, argument `vcov.`).
* Confidence and prediction intervals, and the parameters of the curves, with `predict()`;
  randomised quantile residuals and diagnostic plots (`plot()`).

## Effective doses and benchmark doses

* Effective doses of one or several curves (`ed()`), also with random effects, and the peak of
  hormesis and substrate inhibition curves (`peak()`).
* Benchmark doses with the excess, additional, relative and standard deviation definitions
  (`bmd()`), conditional and marginal BMDs with random effects, and model averaging with AIC or
  BIC weights and bootstrap BMDLs (`bmdma()`).
* Bayesian benchmark doses: the BMD parameterisation (`bmdpar`) with a prior on the BMD, marginal
  likelihoods (`marglik()`), Bayesian model averaging (`bmdbayes()`), the default priors and models
  of EFSA (2022) (`efsapriors()`, `bmdefsa()`), and informative PERT priors from earlier studies
  (`fitpert()`).
