# Starting values for a drmTMB model

Computes starting values for the fixed-effect parameters of a
dose-response model, which can be inspected, modified, and passed to the
`start` argument of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

## Usage

``` r
drmTMBstart(
  form,
  fform = NULL,
  data,
  start = NULL,
  fix = NULL,
  weights = NULL,
  family = "gaussian",
  model = "logistic",
  alpha = NULL,
  smooth = NULL,
  parlink = NULL,
  curvelink = "identity",
  method = "lm",
  gsl = NULL,
  gsu = NULL,
  gsn = NULL
)
```

## Arguments

- form:

  formula with name of response and dose scale

- fform:

  fixed effects for the model parameters: a formula of the form b2 + b4
  ~ Stress, or a list of such formulas, e.g. list(b2 ~ Stress, b4 ~
  Stress + Time). Parameters not listed have a single intercept (~ 1). A
  list of five right-hand side formulas, one for each of b1-b5 (X1 = ~
  1, ..., X5 = ~ 1), is also accepted.

- data:

  a data.frame object. Observations with missing values in the variables
  of the model or in the weights are left out.

- start:

  starting values, as in
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md):
  named by parameter, e.g. c(b5=0), or one entry for each coefficient.
  Only NA entries are replaced.

- fix:

  logical, the coefficients fixed at their starting values, as in
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md):
  named by parameter, e.g. c(b5=TRUE), or one entry for each
  coefficient. Parameters not used by the model are always fixed, at 0
  unless given in start, and b1 of event data at 0.

- weights:

  a vector of weights, one for each row of data.

- family:

  the distribution of the response, see
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

- model:

  the dose-response model, see
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

- alpha:

  the exponent of the hormesis term of the cedergreen model, see
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

- smooth:

  smoothing of the threshold of the nec model, relative to the smallest
  interval between doses, see
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md);
  without random effects, the default is 0, an exact threshold.

- parlink:

  links of the parameters b1-b5, see
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

- curvelink:

  link of the curve, see
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

- method:

  "lm" (default), fast starting values from linearised models, or
  "grid", a grid search of the likelihood, see Details. Random effects
  are not used; drmTMB starts their standard deviations separately.

- gsl:

  a vector of lower bounds for the grid search, one entry for each model
  coefficient, on the scale of the coefficients (log scale for b4 with
  the log link); NA (or NULL) for the default range.

- gsu:

  a vector of upper bounds for the grid search, one entry for each model
  coefficient; NA (or NULL) for the default range.

- gsn:

  number of grid points for each parameter, a single number or a vector
  with one entry for each model coefficient (default 10). A parameter
  with gsn = 1 stays at its "lm" starting value.

## Value

a vector of starting values named by the coefficients, as `coef` of the
fit without the variance parameters, including fixed and unused
parameters; it can be changed and passed to the argument start of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

## Details

Method "lm" takes the asymptotes from the mean response at the lowest
and highest dose, transforms the response to the relative scale between
the asymptotes, and estimates b3 and b4 by linear regression on a
linearised version of the model function. The asymmetry parameter b5
starts at 0. Models without such a linearisation use other strategies:
for asymptoticregression, exponentialdecay and gamma, b4 is profiled
over a grid, with the asymptotes estimated by least squares for each
value; for haldane, b4 and the inhibition b5; for nec, thresholds within
each interval between the doses; for multistage, the polynomial by least
squares; and for exponential and power, the power b3 over a grid. With
covariates in `fform`, the pooled starting value of each parameter is
projected onto its design matrix, giving the intercept or each group the
pooled value and contrasts a value of zero.

Method "grid" evaluates the negative log-likelihood (without random
effects) at all combinations of a grid of parameter values and returns
the best grid point; parameters not included in the grid start at the
values of method "lm". Without bounds in `gsl` and `gsu`, default ranges
are used for b3 (the sign of the "lm" starting value, with 1/4 to 4
times its steepness, log-spaced), b4 (the dose range extended at both
ends; log-spaced for models on the log-dose scale, and on the log scale
with the log link for b4), and b5 (-2 to 2, logistic, loglogistic and
nec; for braincousens and cedergreen, from 0 to four times the hormesis
effect that equals the difference of the asymptotes at dose b4). Default
ranges are only used for parameters with a single coefficient in
`fform`; the asymptotes b1 and b2 are only included in the grid with
user-specified bounds.

## See also

[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)

## Examples

``` r
# starting values of the log-logistic model, as used by drmTMB
(sv <- drmTMBstart(rootl ~ conc, data = ryegrass, model = "loglogistic"))
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept) b5:(Intercept) 
#>    -0.06896726     8.14026086    -1.71197137     3.74376044     0.00000000 
# changed and passed to drmTMB
sv[c("b1:(Intercept)", "b5:(Intercept)")] <- 0
drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic", start = sv,
       fix = c(TRUE, FALSE, FALSE, FALSE, TRUE))
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, start = sv, fix = c(TRUE, 
#>     FALSE, FALSE, FALSE, TRUE), model = "loglogistic")
#> 
#> Model:  loglogistic
#> Family: gaussian, identity link for the asymptotes
#> 
#> Coefficients:
#> b2:(Intercept)  b3:(Intercept)  b4:(Intercept)  
#>          7.855          -2.470           3.263  
#> Fixed: b1:(Intercept) = 0, b5:(Intercept) = 0
#> 
#> Residual standard deviation: 0.5616
```
