# Model-averaged benchmark dose

Benchmark doses (see
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md))
averaged over several dose-response models fitted to the same data, with
AIC or BIC weights, and bootstrap confidence intervals: the BMDL is the
lower limit of the two-sided 90% interval, `confint(bmdma(...))`.

## Usage

``` r
bmdma(
  ...,
  bmr = 0.1,
  def = c("excess", "additional", "relative", "sd"),
  newdata = NULL,
  scale = c("marginal", "conditional"),
  average = c("curve", "bmd"),
  weights = c("AIC", "BIC"),
  nodes = NULL,
  nsim = 500,
  seed = NULL,
  progress = interactive()
)

# S3 method for class 'bmdma'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'bmdma'
coef(object, ...)

# S3 method for class 'bmdma'
confint(object, parm = NULL, level = 0.9, method = c("boot", "buckland"), ...)
```

## Arguments

- ...:

  drmTMB objects, fitted to the same response, doses and weights; the
  models and families may differ (e.g. gaussian and lognormal, whose
  AICs are comparable). Models are labelled by their argument names
  (e.g. LL = fit1) or, if unnamed, by the expressions. The fits are
  maximum likelihood fits without priors; for fits with priors, see
  [`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md).

- bmr, def, newdata, scale, nodes:

  the benchmark response and the definition, the covariates of the
  curves, and the BMD of the marginal or conditional curve with random
  effects, as in
  [`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md);
  newdata cannot give fitted groups of random effects, and `def="sd"`
  with `average="curve"` needs gaussian families

- average:

  "curve" (default) for the BMD of the model-averaged curve, the
  weighted mean of the curves of the models, or "bmd" for the weighted
  mean of the BMDs of the models

- weights:

  the information criterion of the model weights, "AIC" (default) or
  "BIC"

- nsim:

  number of bootstrap samples; 0 for no bootstrap (then there is no BMDL
  for average="curve")

- seed:

  an optional seed for the random number generator

- progress:

  logical, show a progress bar

- x, object:

  a bmdma object

- digits:

  number of significant digits in the printed output

- parm:

  BMDs, names or positions in `coef`; default all

- level:

  the confidence level; 0.90 by default, for the BMDL and BMDU

- method:

  "boot" for bootstrap percentile intervals, or "buckland" for Wald
  intervals with the standard error of Buckland et al. (1997), for
  `average="bmd"` only; by default "boot" if the bootstrap was run (nsim
  \> 0), and "buckland" otherwise

## Value

bmdma returns a bmdma object, a list with the averaged BMDs
(`estimates`), the table of the models (`models`: model, family, logLik,
df, AIC or BIC, weight), the
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)
results of the models (`bmds`), the averaged BMDs of the bootstrap
samples (`boot`, with NA rows for samples left out), their weights
(`bootweights`, NA for failed fits) and the models the samples were
simulated from (`source`), and the arguments. coef returns the averaged
BMDs and confint their confidence intervals.

## Details

AIC weights were recommended by EFSA (2017); the current guidance of
EFSA (2022) uses Bayesian model averaging, see
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md)
and
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md).

Model k gets the weight exp(-IC_k/2), relative to the sum over the
models, with the information criterion IC_k, AIC_k = -2 logLik_k + 2
df_k or BIC_k = -2 logLik_k + log(n) df_k, with the number of
observations n (`nobs`), as in Ritz, Gerhard and Hothorn (2013). With
`average="curve"`, the BMD is that of the averaged curve, the weighted
mean of the curves (the marginal mean curves with `scale="marginal"` and
random effects), with the benchmark response relative to the averaged
response at dose zero, as in EFSA (2017); with `def="sd"`, the standard
deviation is the weighted mean of the residual standard deviations. With
`average="bmd"`, it is the weighted mean of the BMDs of the models;
models whose curve does not reach the benchmark response are left out,
with the weights of the others rescaled.

The bootstrap samples the responses from the mixture of the models: for
each sample, a model is drawn with its weight and a response simulated
from it, with new random effects (see
[`simulate.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)).
All models are fitted to each sample, starting at their estimates, and
the weights and the averaged BMD are computed again. Models whose fit
fails are left out of a sample, with the weights of the others rescaled;
samples without any fit are left out. Models of the lognormal family get
the weight 0 in samples with non-positive responses (e.g. simulated from
a gaussian model), whose likelihood is zero under the lognormal family.
The confidence intervals are percentile intervals of the averaged BMDs.

## References

Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model
selection: an integral part of inference. Biometrics 53, 603-618.

Ritz, C., Gerhard, D. and Hothorn, L. A. (2013). A unified framework for
benchmark dose estimation applied to mixed models and model averaging.
Statistics in Biopharmaceutical Research 5, 79-90.

EFSA Scientific Committee (2017). Update: use of the benchmark dose
approach in risk assessment. EFSA Journal 15(1), 4658.

## See also

[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md),
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)

## Examples

``` r
fit <- function(model, family = "gaussian")
  drmTMB(rootl ~ conc, data = ryegrass, model = model, family = family,
         start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
LL <- fit("loglogistic")
W1 <- fit("weibull1")
W2 <- fit("weibull2")
# AIC weights and the BMD of the averaged curve, without bootstrap
bmdma(LL = LL, W1 = W1, W2 = W2, bmr = 0.05, def = "relative", nsim = 0)
#> Model-averaged benchmark doses, relative deviation, averaged curves:
#> 
#>          model   family logLik df   AIC   weight    BMD
#> LL loglogistic gaussian -18.60  4 45.21 0.138258 0.9909
#> W1    weibull1 gaussian -22.22  4 52.44 0.003717 0.7538
#> W2    weibull2 gaussian -16.78  4 41.56 0.858025 1.3520
#> 
#>           BMD
#> 0.05:lf1 1.31
# \donttest{
# with a bootstrap BMDL
ma <- bmdma(LL = LL, W1 = W1, W2 = W2, bmr = 0.05, def = "relative", nsim = 100, seed = 1)
ma
#> Model-averaged benchmark doses, relative deviation, averaged curves:
#> 
#>          model   family logLik df   AIC   weight    BMD
#> LL loglogistic gaussian -18.60  4 45.21 0.138258 0.9909
#> W1    weibull1 gaussian -22.22  4 52.44 0.003717 0.7538
#> W2    weibull2 gaussian -16.78  4 41.56 0.858025 1.3520
#> 
#>           BMD   BMDL  BMDU
#> 0.05:lf1 1.31 0.8222 1.492
#> 
#> Bootstrap percentile intervals (100 samples), level 90%
confint(ma)
#>                5 %     95 %
#> 0.05:lf1 0.8221607 1.491713
# }
```
