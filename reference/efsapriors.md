# Default priors of the EFSA Bayesian benchmark dose modelling

The default (weakly informative) priors of the Bayesian benchmark dose
modelling of EFSA (2022), as in its R package BMABMDR (Kremer et al.
2026), for a model in the BMD parameterisation of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)
(`bmdpar`): PERT priors for the response at dose zero (background), the
response at high doses (maximum response) and the BMD, lognormal priors
for the steepness and for the residual standard deviation. The result
gives the arguments of `drmTMB`;
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)
fits a set of models with these priors and averages them.

## Usage

``` r
efsapriors(
  form,
  data,
  model,
  family = c("gaussian", "lognormal", "binomial"),
  bmr = NULL,
  def = NULL,
  background = NULL,
  maxresponse = NULL,
  bmd = NULL,
  dprior = c("default", "alternative"),
  extended = TRUE
)

# S3 method for class 'efsapriors'
print(x, ...)
```

## Arguments

- form:

  formula of the response and the dose, as in
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)

- data:

  a data.frame

- model:

  one of "weibull2", "weibull1", "loglogistic", "lognormal", "gamma",
  "multistage", and for the binomial family "probit" and "logit"

- family:

  "gaussian", "lognormal" or "binomial"

- bmr, def:

  the benchmark response and its definition (see
  [`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)); by
  default bmr=0.05 with def="relative" for gaussian and lognormal, and
  bmr=0.1 with def="excess" for binomial, the defaults of EFSA (2022)

- background, maxresponse, bmd:

  informative priors: `c(min, mode, max)` or a prior (see
  [`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md));
  NULL for the defaults; maxresponse only for continuous responses

- dprior:

  "default" or "alternative" prior of the steepness

- extended:

  logical, the upper bound of the default BMD prior is three times the
  highest dose (TRUE) or the highest dose

- x:

  an efsapriors object

- ...:

  not used

## Value

an efsapriors object, a list with the EFSA model (name), the arguments
model, family, curvelink, start, fix, prior and bmdpar of `drmTMB`, the
direction of the curve (increasing) and whether it levels off (flat)

## Details

The models are those of EFSA (2022) available in drcTMB: weibull2
(exponential), weibull1 (inverse exponential), loglogistic (Hill),
lognormal (log-normal), gamma, and multistage of degree 2 (quadratic
exponential, LMS two-stage), and for quantal responses also probit and
logit (the probit and logistic models of EFSA, pnorm(b1 + b4 x^b3) and
plogis(b1 + b4 x^b3): the power model with `curvelink="probit"` or
`"logit"`). Families: "gaussian" (normal distribution), "lognormal"
(fitted with `curvelink="log"`: the curve is the log of the median, as
the lognormal models of EFSA) and "binomial" (quantal responses, with b2
fixed at 1).

Default priors, with the response means at dose zero (y0) and at the
highest dose (yN) (geometric means for lognormal), the benchmark
response q and the highest dose dN:

- background, increasing curve: pert(0.001, y0, 2 y0); decreasing:
  pert(0.5 y0, y0, 2 y0); quantal: pert(max(L/2, 1/(10 n0)), max(p0,
  1/(5 n0)), min(3 U/2, 1 - 1/(10 n0))) with the observed proportion p0
  of the n0 controls and its confidence interval (L, U) of prop.test
  (lower bound 1e-8 without events in the controls); shape 4. For the
  probit and logit models the prior is on the background probability
  pnorm(b1) or plogis(b1).

- maximum response (continuous): increasing pert(y0 (1.01 + q), yN, 2
  yN), or pert(y0 (1.01 + q), 3 yN, 6 yN) if the curve does not level
  off at the high doses; decreasing pert(0.5 yN, yN, y0 (1 - q - 0.01)),
  or pert(0.1 yN, 0.5 yN, y0 (1 - q - 0.01)); shape 4. Whether the curve
  levels off follows from a fractional polynomial of degree 2 of the
  means: its last slope on the log-dose scale is below half its largest
  slope. EFSA puts this prior on the fold change, the maximum response
  divided by the mode of the background; here it is on the maximum
  response (the asymptote) itself.

- BMD: pert(0, dN/2, 3 dN) (`extended=TRUE`) or pert(0, dN/2, dN), with
  shape 0.0001, nearly uniform.

- steepness d: \|b3\| lognormal(1, 1) (EFSA truncates log d at 5, which
  removes a probability of 3e-5 and is left out) or, with
  `dprior="alternative"`, lognormal(0.4, sqrt(0.5)), for a sensitivity
  analysis (the power b3 of the probit and logit models has the same
  prior); for decreasing curves of the log-dose models b3 is negative,
  with the prior `negate(lognormal(...))`. The multistage model of
  degree 2 in the BMD parameterisation is the quadratic exponential
  model of EFSA, with the technical parameter log_d (in place of b3),
  with the prior normal(0, 1) of EFSA.

- residual standard deviation: lognormal(log(1.5 s), 0.5), with the mean
  s of the standard deviations within the doses (of log(y) for
  lognormal), i.e. log(1/sigma^2) normal with mean -2 log(1.5 s) and
  standard deviation 1.

For the lognormal family, the priors of the background and the maximum
response are priors of the median, exp(b1) or exp(b2). Informative
priors are given as `c(min, mode, max)` (mode NA for the middle of the
range) or as a prior, e.g. from
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md).

The data need a control group at dose zero and at least two other doses,
replicates at the doses (for the prior of the residual standard
deviation), a positive response for the continuous families, and, for
quantal responses, cbind(affected, unaffected) increasing with the dose.

## References

EFSA Scientific Committee (2022). Guidance on the use of the benchmark
dose approach in risk assessment. EFSA Journal 20(10), 7584.

Kremer, C., Saenen, N. D., Bijnens, K., Faes, C., De Boever, P., Smeets,
K. and Aerts, M. (2026). Compilation of NTP and published dose-response
data: study of informative priors in the EFSA platform for Bayesian
benchmark dose analysis. EFSA Supporting Publications, EN-9902.

## See also

[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md),
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md),
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)

## Examples

``` r
# default priors of EFSA (2022) for the log-logistic model of the ryegrass data
efsapriors(rootl ~ conc, data = ryegrass, model = "loglogistic", family = "gaussian")
#> EFSA priors for the loglogistic model, gaussian family, decreasing curve, levelling off
#> BMD: relative, bmr = 0.05
#> 
#>   b1 ~ pert(0.151667, 0.303333, 7.28438)
#>   b2 ~ pert(3.87467, 7.74935, 15.4987)
#>   b3 ~ negate(lognormal(1, 1))
#>   bmd ~ pert(0, 15, 90, shape=0.0001)
#>   sigma ~ lognormal(-0.374202, 0.5)
# an informative prior of the BMD, by its minimum, mode and maximum
efsapriors(rootl ~ conc, data = ryegrass, model = "loglogistic", family = "gaussian",
           bmd = c(0.1, 1, 5))
#> EFSA priors for the loglogistic model, gaussian family, decreasing curve, levelling off
#> BMD: relative, bmr = 0.05
#> 
#>   b1 ~ pert(0.151667, 0.303333, 7.28438)
#>   b2 ~ pert(3.87467, 7.74935, 15.4987)
#>   b3 ~ negate(lognormal(1, 1))
#>   bmd ~ pert(0.1, 1, 5, shape=0.0001)
#>   sigma ~ lognormal(-0.374202, 0.5)
# the model with these priors, in the BMD parameterisation
e <- efsapriors(rootl ~ conc, data = ryegrass, model = "loglogistic", family = "gaussian")
fit <- drmTMB(rootl ~ conc, data = ryegrass, model = e$model, family = e$family,
              curvelink = e$curvelink, start = e$start, fix = e$fix, prior = e$prior,
              bmdpar = e$bmdpar)
summary(fit)
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, start = e$start, 
#>     fix = e$fix, family = e$family, model = e$model, curvelink = e$curvelink, 
#>     prior = e$prior, bmdpar = e$bmdpar)
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 24
#> Log-likelihood: -16.18 (df = 5), AIC: 42.36, BIC: 48.25
#> 
#> Priors:
#>   b1 ~ pert(0.151667, 0.303333, 7.28438)
#>   b2 ~ pert(3.87467, 7.74935, 15.4987)
#>   b3 ~ negate(lognormal(1, 1))
#>   bmd ~ pert(0, 15, 90, shape=0.0001)
#>   sigma ~ lognormal(-0.374202, 0.5)
#> 
#> Coefficients (posterior modes):
#>                Estimate Std. Error
#> b1:(Intercept)   0.4648      0.196
#> b2:(Intercept)   7.7913      0.176
#> b3:(Intercept)  -2.9634      0.425
#> log_bmd          0.1506      0.151
#> log_bmd is the log benchmark dose (relative, bmr = 0.05, decreasing curve) of the curve of the intercepts: BMD = 1.163
#> Fixed: b5:(Intercept) = 0
#> 
#> Residual standard deviation: 0.4886
#> 
#> Optimiser: nlminb: relative convergence (4), maximum absolute gradient 7.2e-10
```
