# Changelog

## drcTMB 0.1.0

First release.

### Models and data

- [`drmTMB()`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)
  fits nonlinear dose-response models with the Template Model Builder:
  log-logistic, logistic, Weibull (two types), log-normal, Gompertz,
  Michaelis-Menten, asymptotic regression, gamma, multistage, the
  threshold model nec, the hormesis models of Brain-Cousens and
  Cedergreen, the substrate inhibition model of Haldane, the exponential
  and power models of PROAST and BMDS, and a constant model.
- Response distributions: gaussian, lognormal, binomial, beta-binomial,
  beta, poisson, negative binomial (nbinom1, nbinom2), and time-to-event
  data (`family = "event"`).
- Covariates on any parameter (`fform`), random effects on any parameter
  with the syntax of lme4 (`rform`), parameter links (`parlink`) and
  links of the whole curve (`curvelink`).
- Correlated residuals
  ([`cor_ar1()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`cor_car1()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`cor_cs()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`cor_exp()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`cor_gaus()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`cor_spher()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md))
  and variance functions
  ([`var_power()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`var_exp()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
  [`var_ident()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md))
  for the gaussian and lognormal families.
- Priors on all parameters (`prior`, see
  [`?priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)),
  with the estimates as posterior modes.
- Starting values by linearised models or a grid search
  ([`drmTMBstart()`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBstart.md));
  choice of the optimiser and a restart at the estimates
  ([`drmTMBcontrol()`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBcontrol.md)).
- Data sets from the package drcData, which is not on CRAN: chickweed,
  daphnids, earthworms, germination, lettuce, liver.tumor, ryegrass,
  S.alba, selenium, TCDD and vinclozolin.

### Inference

- Wald, profile likelihood and parametric bootstrap intervals
  ([`confint()`](https://rdrr.io/r/stats/confint.html),
  [`profile()`](https://rdrr.io/r/stats/profile.html),
  [`bootstrap()`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md));
  small-sample inference with t quantiles and F-based profiles for
  gaussian and lognormal models, as in drc.
- Likelihood ratio, F and parametric bootstrap tests with
  [`anova()`](https://rdrr.io/r/stats/anova.html).
- Robust (sandwich) covariance, also for clustered data
  ([`vcovrobust()`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md),
  argument `vcov.`).
- Confidence and prediction intervals, and the parameters of the curves,
  with [`predict()`](https://rdrr.io/r/stats/predict.html); randomised
  quantile residuals and diagnostic plots
  ([`plot()`](https://rdrr.io/r/graphics/plot.default.html)).

### Effective doses and benchmark doses

- Effective doses of one or several curves
  ([`ed()`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)),
  also with random effects, and the peak of hormesis and substrate
  inhibition curves
  ([`peak()`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md)).
- Benchmark doses with the excess, additional, relative and standard
  deviation definitions
  ([`bmd()`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)),
  conditional and marginal BMDs with random effects, and model averaging
  with AIC or BIC weights and bootstrap BMDLs
  ([`bmdma()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)).
- Bayesian benchmark doses: the BMD parameterisation (`bmdpar`) with a
  prior on the BMD, marginal likelihoods
  ([`marglik()`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md)),
  Bayesian model averaging
  ([`bmdbayes()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)),
  the default priors and models of EFSA (2022)
  ([`efsapriors()`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md),
  [`bmdefsa()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)),
  and informative PERT priors from earlier studies
  ([`fitpert()`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md)).
