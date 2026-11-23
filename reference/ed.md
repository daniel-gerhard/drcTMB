# Effective dose estimation

Estimates effective doses (ED p, e.g. ED50 or EC50) of fitted
dose-response curves, for one or several curves, with standard errors by
the delta method; confidence intervals follow with `confint`. Not
available for fits in the BMD parameterisation (bmdpar of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)).

## Usage

``` r
ed(
  x,
  respLev = 50,
  linfct = NULL,
  newdata = NULL,
  random = TRUE,
  residual = FALSE,
  type = c("relative", "absolute"),
  vcov. = NULL
)
```

## Arguments

- x:

  a drmTMB object

- respLev:

  effective dose level in percent, default is 50. The effective dose ED
  p is the dose at which the response has changed by p percent of the
  difference between the asymptote at low doses (the control) and the
  asymptote at high doses; for the hormesis models, the dose above the
  peak of the stimulation, see Details. For `type="absolute"`, values of
  the response.

- linfct:

  a matrix of linear functions to identify the parameters of the
  dose-response curve, with 5 rows (b1-b5) and one column for each model
  coefficient, or a list of such matrices; default is the identity
  matrix, which only works without covariates in fform; otherwise give
  newdata. Each matrix is paired with an entry of respLev.

- newdata:

  a data frame with the values of the covariates in fform, one row for
  each dose-response curve, as an alternative to linfct; the dose is not
  needed. With grouping variables of the random effects, the curves of
  fitted groups, see Details. The effective doses are computed at each
  level of respLev for each row, named by the covariate and group values
  (or by the row names of newdata, if set).

- random:

  logical. TRUE (default): rows of newdata with a fitted group give the
  effective dose of the curve of that group; for the other random effect
  terms, the standard errors include the variance of the random effects,
  i.e. they describe the effective dose of a new, randomly drawn group.
  FALSE: the effective doses of the population curve (random effects at
  zero), with the uncertainty of the fixed effects only; profile and
  bootstrap intervals need FALSE (or fitted groups). Without random
  effects, TRUE and FALSE are the same.

- residual:

  logical, if TRUE, the residual variance is included in the standard
  error, for the dose at which a single new observation reaches the
  response level; gaussian family only, and not with a curvelink or a
  variance function

- type:

  "relative" (default) for effective doses relative to the asymptotes,
  or "absolute" for the doses at which the curves reach the responses in
  respLev, on the scale of the response (proportions for binomial,
  beta-binomial and beta). Absolute effective doses depend on the
  asymptotes, and on their random effects; responses outside the range
  of a curve give NA, with a warning. The exponential and power models,
  without an upper asymptote, have absolute effective doses only.

- vcov.:

  the covariance of the estimated parameters for the standard errors: a
  matrix as vcov() of the fit, or a function of the fit that returns
  one, e.g. the sandwich covariance
  [`vcovrobust`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md);
  NULL (default) for vcov() of the fit.
  [`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
  then gives Wald intervals by default.

## Value

An object of class edTMB, a list with the estimates (estimates), their
covariance matrix (variance), the fit and the linear functions of the
curves. The entries are named "level:curve", e.g. "50:lf1" or
"50:time=24h", and "y=value:curve" for absolute effective doses. Use
`coef`, `vcov` and `confint`, see
[`edTMB-methods`](https://daniel-gerhard.github.io/drcTMB/reference/edTMB-methods.md).

## Details

Without random effects, or with `random=FALSE`, the effective doses are
those of the population curves, with all random effects at zero, and
their standard errors include the uncertainty of the fixed effects only.
With `random=TRUE`, the standard errors also include the variance of the
random effects, the variation of the effective dose between groups.

If a row of newdata gives a level of a grouping factor that is in the
fitted data, the effective dose of that row is for the curve of the
group, with its predicted random effects (the conditional modes), as in
[`predict.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md).
These are shrunken towards the population curve, more so for groups with
few observations. Their standard errors include the uncertainty of the
fixed effects and the prediction error of the random effects, from the
joint precision matrix of
[`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). Random effect
terms without a grouping variable in newdata, or with a level that is
not in the fitted data, are treated as above. Profile intervals are not
available for the effective doses of fitted groups, see
[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md).

For the hormesis models braincousens and cedergreen, the effective dose
ED p is the dose above the peak of the stimulation at which the response
has decreased from the response at dose zero (b2) by p percent of the
difference to the lower asymptote b1, as in drc; with `type="absolute"`,
the dose above the peak at which the curve reaches the response. It has
no closed form: it is found numerically, and its derivatives, for the
standard errors, follow from the implicit function theorem. The models
need a decreasing curve, with b3 \< 0. With `respLev=0`, the effective
dose is the dose above the peak at which the response returns to the
response at dose zero, the limit of the stimulation (the zero equivalent
point); see
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md) for
the peak itself.

The effective dose ED p is relative to the range between the asymptotes,
whereas the relative BMD of
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)
(def="relative") is relative to the response at dose zero; they coincide
only for a lower asymptote of 0. For the event family, the relative
effective dose is the time by which p percent of the units with an event
have had it (b2 is the proportion with an event); the time by which a
proportion q of all units has had an event is given by `type="absolute"`
with respLev = q.

## See also

[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md),
[`profile.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/profile.edTMB.md),
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md),
[`edTMB-methods`](https://daniel-gerhard.github.io/drcTMB/reference/edTMB-methods.md),
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md),
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
# ED10, ED50 and ED90, with profile likelihood intervals
(e <- ed(mod, respLev = c(10, 50, 90)))
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 10:lf1     1.34     0.189
#> 50:lf1     3.26     0.196
#> 90:lf1     7.94     1.055
confint(e)
#> Computing profile confidence intervals ...
#>            2.5 %    97.5 %
#> 10:lf1 0.9813186  1.774243
#> 50:lf1 2.8791625  3.703885
#> 90:lf1 6.1400175 10.652850
# doses at which the root length is 2 and 4 cm
ed(mod, respLev = c(2, 4), type = "absolute")
#> Effective dose estimates:
#> 
#>         Estimate Std.Error
#> y=2:lf1     5.04     0.414
#> y=4:lf1     3.22     0.173

# effective doses of several curves, given by the covariates in newdata
dmod <- drmTMB(cbind(no, total - no) ~ dose, data = daphnids, model = "loglogistic",
               family = "binomial", fform = b3 + b4 ~ 0 + time,
               start = c(b1 = 0, b2 = 1, b4 = 1000, b5 = 0),
               fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
ed(dmod, respLev = 50, newdata = data.frame(time = c("24h", "48h")))
#> Effective dose estimates:
#> 
#>             Estimate Std.Error
#> 50:time=24h     5134      1057
#> 50:time=48h     1509       188
```
