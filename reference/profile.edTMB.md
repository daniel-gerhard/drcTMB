# Likelihood profiles of effective doses, peak doses and benchmark doses

Profile likelihood of the estimates of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md) and
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md), used
by `confint(method="profile")`; the result can be reused with `confint`
and plotted.

## Usage

``` r
# S3 method for class 'edTMB'
profile(
  fitted,
  parm = NULL,
  level = 0.99,
  maxit = 100,
  progress = interactive(),
  ...
)
```

## Arguments

- fitted:

  an edTMB object

- parm:

  effective doses to profile, names or positions in `coef(fitted)`;
  default all

- level:

  the highest confidence level the profiles are computed for

- maxit:

  maximum number of steps in each direction

- progress:

  logical, show a progress bar over the profiled effective doses

- ...:

  ignored

## Value

a profile.edTMB object, which has the confint and plot methods of
[`profile.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md)

## Details

Profiles the (marginal) log-likelihood for the effective doses. At each
value of an effective dose, the model is reparameterised: the effective
dose is a function of b3, b4 and b5 that can be solved for b4, so one
coefficient of b4 is computed from the effective dose and the other
coefficients, which are re-optimised.

The dose at the peak of a hormesis curve
([`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md)) is
profiled in the same way, as b4 follows from the peak dose and the other
parameters; profiles are not available for the response at the peak.

The profiles are for the effective doses of the population curves,
without the random effects; the variation between groups (`random=TRUE`
in [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)) and
the residual variation are not included. For the models on the log-dose
scale and the models without b3, the effective doses are profiled on the
log scale. Gaussian and lognormal models with small-sample inference
(without random effects, see
[`drmTMBcontrol`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBcontrol.md))
use the F-based cutoff of
[`profile.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.drmTMB.md).

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
pe <- profile(ed(mod, respLev = c(10, 50)))
confint(pe)
#>            2.5 %   97.5 %
#> 10:lf1 0.9813209 1.774205
#> 50:lf1 2.8792091 3.703799
plot(pe)
```
