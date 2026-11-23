# Confidence intervals for drmTMB models and effective doses

Profile likelihood intervals (default), Wald intervals, or parametric
bootstrap intervals for the coefficients and variance parameters of a
drmTMB model, or for effective doses estimated by
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md).

## Usage

``` r
# S3 method for class 'drmTMB'
confint(
  object,
  parm = NULL,
  level = 0.95,
  method = c("profile", "wald", "boot"),
  nsim = 500,
  seed = NULL,
  boot = NULL,
  vcov. = NULL,
  ...
)

# S3 method for class 'edTMB'
confint(
  object,
  parm = NULL,
  level = if (identical(object$type, "bmd")) 0.9 else 0.95,
  method = c("profile", "wald", "boot"),
  nsim = 500,
  seed = NULL,
  boot = NULL,
  ...
)
```

## Arguments

- object:

  a drmTMB or edTMB object

- parm:

  parameters or effective doses, names or positions in `coef(object)`;
  default all

- level:

  the confidence level: 0.95 by default, and 0.90 for the benchmark
  doses of
  [`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md),
  whose lower and upper limits of the two-sided 90% interval are the
  BMDL and BMDU of EFSA and the US EPA (as for
  [`bmdma`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)
  and
  [`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md))

- method:

  "profile" for profile likelihood intervals, "wald" for Wald intervals,
  or "boot" for percentile intervals of a parametric bootstrap:
  responses are simulated from the fitted model (with new random
  effects, see
  [`simulate.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md))
  and the model is fitted again to each, starting at the estimates. Fits
  that do not converge are left out, with a message.

- nsim:

  number of bootstrap simulations

- seed:

  an optional seed for the random number generator of the bootstrap

- boot:

  a bootstrap sample of the model from
  [`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md),
  used instead of new simulations; implies `method="boot"`

- vcov.:

  for drmTMB objects, the covariance of the estimated parameters for
  Wald intervals: a matrix as vcov() of the fit, or a function of the
  fit that returns one, e.g. the sandwich covariance
  [`vcovrobust`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md);
  implies `method="wald"`. Effective doses computed with vcov. (see
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)) get
  Wald intervals by default.

- ...:

  further arguments passed to
  [`profile.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md),
  [`profile.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.edTMB.md)
  or
  [`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md),
  e.g. `progress=FALSE`

## Value

a matrix with the lower and upper confidence limits

## Details

Standard deviations of the random effects, the residual standard
deviation and the precision are reported on their natural scale; Wald
intervals for them are computed on the log scale. The correlation of a
term with two random effects is reported as a correlation; the
parameters of larger correlation matrices on the unconstrained scale of
`coef`.

For effective doses, peaks and benchmark doses, the available methods
are:

- profile: effective doses of the population curves (`random=FALSE` in
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md), or
  models without random effects), peak doses, and BMDs that are not
  averaged over random effects (see
  [`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md) for
  the exceptions);

- bootstrap: as profile, and also the curves of fitted groups (`newdata`
  in [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)),
  with the random effects of the groups drawn from their conditional
  distribution at each bootstrap fit;

- Wald only: effective doses with `random=TRUE` for random effect terms
  without a fitted group, or with `residual=TRUE`, and the response at a
  peak.

Where a method is not available, Wald intervals are returned with a
message. Wald intervals of effective doses and BMDs are symmetric on the
dose scale and can have negative lower limits; with small-sample
inference (see
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md))
they use t quantiles with n - p degrees of freedom. Bootstrap intervals
need many fits of the model; with random effects this can take several
minutes. Computing profiles can take a while for models with random
effects; see
[`profile.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md)
and
[`profile.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.edTMB.md),
whose results can be reused with `confint`.

## See also

[`profile.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md),
[`profile.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.edTMB.md),
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md),
[`vcovrobust`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md),
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
confint(mod)
#> Computing profile confidence intervals ...
#>                    2.5 %     97.5 %
#> b2:(Intercept)  7.436904  8.2868123
#> b3:(Intercept) -3.398753 -1.8910294
#> b4:(Intercept)  2.879209  3.7038016
#> sigma           0.405628  0.7178603
confint(mod, method = "wald")
#>                     2.5 %     97.5 %
#> b2:(Intercept)  7.4303933  8.2804619
#> b3:(Intercept) -3.1808940 -1.7597546
#> b4:(Intercept)  2.8549066  3.6718047
#> sigma           0.3958726  0.6970697
confint(ed(mod, respLev = 50), level = 0.9)
#> Computing profile confidence intervals ...
#>             5 %     95 %
#> 50:lf1 2.942286 3.622696
```
