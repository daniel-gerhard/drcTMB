# Benchmark dose estimation

The benchmark dose (BMD) is the dose at which the response has changed
from the response at dose zero by a specified amount, the benchmark
response (BMR). Its lower confidence limit, the BMDL, is the reference
point of regulatory risk assessment; with the convention of EFSA and the
US EPA, the lower limit of a two-sided 90% interval, `confint(bmd(fit))`
(confint uses the level 0.90 for BMDs). Not available for fits in the
BMD parameterisation (bmdpar of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)).

## Usage

``` r
bmd(
  x,
  bmr = 0.1,
  def = c("excess", "additional", "relative", "sd"),
  linfct = NULL,
  newdata = NULL,
  scale = c("marginal", "conditional"),
  nodes = NULL,
  vcov. = NULL
)
```

## Arguments

- x:

  a drmTMB object

- bmr:

  the benchmark response, see `def`; a vector gives a BMD for each value

- def:

  the definition of the benchmark response, see Details: "excess"
  (extra) risk or "additional" (added) risk for quantal responses,
  "relative" deviation, or "sd", the change in standard deviations of
  the response (gaussian and lognormal families, models without random
  effects and variance function). The default "excess" is for quantal
  responses only; for continuous responses choose "relative" (the EFSA
  default, e.g. bmr=0.05) or "sd" (the US EPA default, bmr=1).

- linfct:

  a matrix of linear functions to identify the parameters of the
  dose-response curve, or a list of such matrices, as in
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md); each
  matrix is paired with an entry of bmr

- newdata:

  a data frame with the covariates of the curves, as an alternative to
  linfct, as in
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md); rows
  with a level of a grouping factor of the fitted data give the BMD of
  the curve of that group

- scale:

  with random effects, "marginal" (default) for the BMD of the marginal
  (population-averaged) mean curve, the mean response over the
  distribution of the random effects, or "conditional" for the BMD of
  the typical curve, with the random effects at zero

- nodes:

  number of Gauss-Hermite quadrature nodes per random coefficient for
  the marginal mean; the default decreases with the number of random
  coefficients, from 30 for one to 5 for five or more

- vcov.:

  the covariance of the estimated parameters, as in
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)

## Value

an edTMB object, as from
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md), with
the BMDs and their covariance matrix, entries named "bmr:curve"; print
shows the estimates and their standard errors (delta method), and
`confint` the BMDL and BMDU (level 0.90 by default).

## Details

With the response mu0 at dose zero, the BMD is the dose at which the
curve reaches the response

- excess (extra) risk: mu0 + bmr (1 - mu0) for an increasing curve, and
  mu0 (1 - bmr) for a decreasing curve, i.e. the change relative to the
  proportion that can still change;

- additional (added) risk: mu0 + bmr, or mu0 - bmr for a decreasing
  curve;

- relative deviation: mu0 (1 + bmr), or mu0 (1 - bmr) for a decreasing
  curve;

- sd: mu0 + bmr sigma, or mu0 - bmr sigma for a decreasing curve, with
  the residual standard deviation sigma (for weights of 1); for the
  lognormal family, mu0 exp(bmr sigma), or mu0 exp(-bmr sigma), with the
  standard deviation sigma of log(y), as in BMDS.

The definitions excess and additional are for responses that are
proportions (binomial, beta-binomial and beta families). For excess
risk, the BMD of a curve from zero to one is the effective dose
`ed(x, 100*bmr)`. The response at dose zero, and sigma, are estimated,
and their uncertainty is included in the standard errors and confidence
intervals. Curves that do not reach the response give NA, with a
warning. The direction of the change is that of the fitted curve. With
small-sample inference (see
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)),
sigma for def="sd" is the unbiased residual standard deviation.

The relative BMD is relative to the response at dose zero, whereas the
effective dose ED p of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md) is
relative to the range between the asymptotes; they coincide only for a
lower asymptote of 0. BMDs are not available for the hormesis models and
for event data (see
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)).

For the lognormal family, mu0 is the median at dose zero; as the mean is
the median times a constant, the relative BMD is the same for the median
and the mean.

With random effects, the BMD depends on the curve it refers to. The
conditional BMD is that of the typical curve, with all random effects at
zero, as the effective doses of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md) with
`random=FALSE`. The marginal BMD is that of the mean response of the
population of groups, the mean curve over the distribution of the random
effects, found by Gauss-Hermite quadrature. Both are the same if the
random effects enter the mean linearly, i.e. on b1 and b2 with the
identity link; otherwise they differ, and the marginal BMD also depends
on the variances of the random effects, whose uncertainty is included in
its standard errors. The standard errors are for the uncertainty of the
estimates, not the variation between groups. Rows of newdata with fitted
groups give the BMDs of the curves of these groups, with their predicted
random effects; with `scale="marginal"`, the terms without a fitted
group are averaged over.

Confidence intervals: see
[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md).
Profile likelihood intervals are available for BMDs without averaging
over random effects, except for the models logistic and gompertz and the
smoothed threshold of nec, whose response at dose zero depends on b4;
otherwise Wald intervals are returned, with a message, and bootstrap
intervals are available.

## References

EFSA Scientific Committee (2022). Guidance on the use of the benchmark
dose approach in risk assessment. EFSA Journal 20(10), 7584.

US EPA (2012). Benchmark Dose Technical Guidance. EPA/100/R-12/001, Risk
Assessment Forum, Washington, DC.

## See also

[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md),
[`bmdma`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md),
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md),
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)

## Examples

``` r
# 10% extra risk of immobility of daphnids after 48 hours
d48 <- subset(daphnids, time == "48h")
dmod <- drmTMB(cbind(no, total - no) ~ dose, data = d48, model = "loglogistic",
               family = "binomial", start = c(b1 = 0, b2 = 1, b4 = 1000, b5 = 0),
               fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
(b <- bmd(dmod, bmr = 0.1, def = "excess"))
#> Benchmark doses, excess risk:
#> 
#>         Estimate Std.Error
#> 0.1:lf1      460       103
# the BMDL: lower limit of the two-sided 90% interval
confint(b)
#> Computing profile confidence intervals ...
#>              5 %     95 %
#> 0.1:lf1 300.1731 638.5276

# 5% decrease of the root length of ryegrass, and a decrease by one standard deviation
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
bmd(mod, bmr = 0.05, def = "relative")
#> Benchmark doses, relative deviation:
#> 
#>          Estimate Std.Error
#> 0.05:lf1    0.991     0.178
bmd(mod, bmr = 1, def = "sd")
#> Benchmark doses, standard deviations:
#> 
#>       Estimate Std.Error
#> 1:lf1     1.16     0.205
```
