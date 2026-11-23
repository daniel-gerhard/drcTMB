# drcTMB: Fitting Dose-Response Curves Using the Template Model Builder

Fits non-linear dose-response models (log-logistic, Weibull, log-normal,
Gompertz, hormesis, threshold, Michaelis-Menten, and the gamma,
multistage, exponential and power models of the benchmark dose software
of the US EPA and of 'PROAST'), optionally with covariates and random
effects on the parameters, using the Template Model Builder ('TMB',
Kristensen et al. (2016)
[doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05) )
with the Laplace approximation for the random effects. The interface
follows the 'drc' package (Ritz et al. (2015)
[doi:10.1371/journal.pone.0146021](https://doi.org/10.1371/journal.pone.0146021)
). The response can be Gaussian, lognormal, binomial, beta-binomial,
beta, Poisson, negative binomial or time-to-event counts, and Gaussian
and lognormal residuals can be correlated or heteroscedastic. Estimation
is by maximum likelihood, or by the posterior mode with priors.
Effective doses and benchmark doses come with standard errors by
automatic differentiation and with profile likelihood, bootstrap or
robust sandwich intervals; benchmark doses can be averaged over models
with information criteria or by Bayesian model averaging as recommended
by the European Food Safety Authority (EFSA Scientific Committee (2022)
[doi:10.2903/j.efsa.2022.7584](https://doi.org/10.2903/j.efsa.2022.7584)
).

## Details

The main functions are:

- [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md):
  fit a dose-response model; methods in
  [`drmTMB-methods`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
  diagnostics with
  [`plot.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/plot.drmTMB.md);

- [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
  [`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md):
  effective doses and the peak of hormesis curves;

- [`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md),
  [`bmdma`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md):
  benchmark doses and their model averaging;

- [`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
  [`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md):
  Bayesian model averaging of benchmark doses;

- [`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md),
  [`profile.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md),
  [`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md),
  [`vcovrobust`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md):
  intervals and standard errors;

- [`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md):
  model comparison.

[`vignette("drcTMB", package = "drcTMB")`](https://daniel-gerhard.github.io/drcTMB/articles/drcTMB.md)
introduces the package with examples.

## Author

**Maintainer**: Daniel Gerhard <daniel.gerhard@canterbury.ac.nz>
([ORCID](https://orcid.org/0000-0002-9336-3454))

Authors:

- Daniel Gerhard <daniel.gerhard@canterbury.ac.nz>
  ([ORCID](https://orcid.org/0000-0002-9336-3454))
