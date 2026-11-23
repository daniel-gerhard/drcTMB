# Package index

## Model fitting

- [`drmTMB()`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)
  : Fit dose-response models with TMB
- [`drmTMBcontrol()`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBcontrol.md)
  : Control parameters of drmTMB
- [`drmTMBstart()`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBstart.md)
  : Starting values for a drmTMB model
- [`cor_ar1()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`cor_car1()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`cor_cs()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`cor_exp()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`cor_gaus()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`cor_spher()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`var_power()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`var_exp()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  [`var_ident()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  : Residual correlation structures and variance functions
- [`negate()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  [`normal()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  [`student_t()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  [`lognormal()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  [`pert()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  [`half_normal()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  [`half_student_t()`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)
  : Prior distributions for drmTMB models

## Methods for fitted models

- [`print(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`summary(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`print(`*`<summary.drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`logLik(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`nobs(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`df.residual(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`fitted(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`residuals(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`predict(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`coef(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`vcov(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  [`simulate(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)
  : Methods for fitted drmTMB models
- [`anova(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
  : Likelihood ratio tests for drmTMB models
- [`confint(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
  [`confint(`*`<edTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
  : Confidence intervals for drmTMB models and effective doses
- [`plot(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/plot.drmTMB.md)
  : Diagnostic plots for drmTMB models
- [`profile(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md)
  [`confint(`*`<profile.drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md)
  [`plot(`*`<profile.drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md)
  : Likelihood profiles of a drmTMB model
- [`vcovrobust()`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md)
  [`estfun(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md)
  [`bread(`*`<drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md)
  : Robust (sandwich) covariance of a drmTMB fit
- [`bootstrap()`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
  [`print(`*`<bootstrap.drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
  [`confint(`*`<bootstrap.drmTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
  : Parametric bootstrap of a drmTMB model
- [`marglik()`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md)
  : Marginal likelihood of a drmTMB model with priors

## Effective doses and peaks

- [`ed()`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md) :
  Effective dose estimation
- [`peak()`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md)
  : Peak of hormesis and substrate inhibition curves
- [`print(`*`<edTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/edTMB-methods.md)
  [`coef(`*`<edTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/edTMB-methods.md)
  [`vcov(`*`<edTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/edTMB-methods.md)
  : Methods for effective doses, peaks and benchmark doses
- [`profile(`*`<edTMB>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/profile.edTMB.md)
  : Likelihood profiles of effective doses, peak doses and benchmark
  doses

## Benchmark doses

- [`bmd()`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md) :
  Benchmark dose estimation
- [`bmdma()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)
  [`print(`*`<bmdma>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)
  [`coef(`*`<bmdma>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)
  [`confint(`*`<bmdma>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)
  : Model-averaged benchmark dose
- [`bmdbayes()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)
  [`coef(`*`<bmdbayes>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)
  [`confint(`*`<bmdbayes>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)
  [`bmdbayesmodels()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)
  [`print(`*`<bmdbayes>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)
  : Bayesian model-averaged benchmark dose
- [`bmdefsa()`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)
  : Bayesian model-averaged BMD with the EFSA defaults
- [`efsapriors()`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md)
  [`print(`*`<efsapriors>`*`)`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md)
  : Default priors of the EFSA Bayesian benchmark dose modelling
- [`fitpert()`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md)
  : A PERT prior from a posterior distribution

## Data sets

- [`chickweed`](https://daniel-gerhard.github.io/drcTMB/reference/chickweed.md)
  : Germination of common chickweed
- [`daphnids`](https://daniel-gerhard.github.io/drcTMB/reference/daphnids.md)
  : Immobility of daphnids
- [`earthworms`](https://daniel-gerhard.github.io/drcTMB/reference/earthworms.md)
  : Earthworm toxicity test
- [`germination`](https://daniel-gerhard.github.io/drcTMB/reference/germination.md)
  : Germination of three crops
- [`lettuce`](https://daniel-gerhard.github.io/drcTMB/reference/lettuce.md)
  : Hormesis in lettuce plants
- [`liver.tumor`](https://daniel-gerhard.github.io/drcTMB/reference/liver.tumor.md)
  : Liver tumour incidence after exposure to TCDD (NTP 2006)
- [`ryegrass`](https://daniel-gerhard.github.io/drcTMB/reference/ryegrass.md)
  : Effect of ferulic acid on the root growth of ryegrass
- [`S.alba`](https://daniel-gerhard.github.io/drcTMB/reference/S.alba.md)
  : Potency of two herbicides in white mustard
- [`selenium`](https://daniel-gerhard.github.io/drcTMB/reference/selenium.md)
  : Toxicity of four types of selenium
- [`TCDD`](https://daniel-gerhard.github.io/drcTMB/reference/TCDD.md) :
  Liver tumour incidence after exposure to TCDD (Kociba et al. 1978)
- [`vinclozolin`](https://daniel-gerhard.github.io/drcTMB/reference/vinclozolin.md)
  : Vinclozolin in an androgen receptor reporter gene assay
