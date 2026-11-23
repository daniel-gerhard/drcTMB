# Prior distributions for drmTMB models

Prior distributions for the `prior` argument of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).
With priors, the estimates are the posterior mode (maximum a
posteriori), and the standard errors are from the curvature of the log
posterior, see
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

## Usage

``` r
negate(p)

normal(mean = 0, sd = 1)

student_t(df = 3, location = 0, scale = 1)

lognormal(meanlog = 0, sdlog = 1)

pert(min, mode, max, shape = 4)

half_normal(scale = 1)

half_student_t(df = 3, scale = 1)
```

## Arguments

- p:

  a prior, for negate

- mean, sd:

  mean and standard deviation of the normal distribution

- df:

  degrees of freedom of the Student t distribution

- location, scale:

  location and scale of the Student t distribution, or scale of the half
  normal and half Student t distributions

- meanlog, sdlog:

  mean and standard deviation of the log of a lognormal parameter

- min, mode, max:

  minimum, mode and maximum of the PERT distribution

- shape:

  the shape (weight of the mode) of the PERT distribution

## Value

a drcprior object, for the argument prior of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)

## Details

The priors are given by name in `prior=list(...)` of drmTMB: b1-b5,
coefficient names, "bmd" (with bmdpar), "sd" for the standard deviations
of the random effects, "sigma", "phi" or "theta" for the residual
parameter, and the parameters of correlation structures and variance
functions; see the argument prior of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).
The prior `lognormal()` is not related to `family="lognormal"` or
`model="lognormal"`.

- `normal(mean, sd)`: normal distribution.

- `student_t(df, location, scale)`: Student t distribution with df
  degrees of freedom, shifted by location and scaled by scale.

- `lognormal(meanlog, sdlog)`: lognormal distribution of a positive
  parameter, with log(parameter) normal with mean meanlog and standard
  deviation sdlog.

- `pert(min, mode, max, shape)`: the PERT distribution of EFSA (2022)
  for benchmark dose modelling, a beta distribution on (min, max) with
  mode `mode` and the shape parameters 1 + shape (mode - min)/(max -
  min) and 1 + shape (max - mode)/(max - min); shape = 4 by default.

- `half_normal(scale)`, `half_student_t(df, scale)`: normal and Student
  t distributions folded at 0, for standard deviations.

- `negate(p)`: the prior p of the negative of the parameter, e.g.
  `b3=negate(lognormal(1, 1))` for a negative steepness b3 (a decreasing
  curve) whose absolute value is lognormal.

A prior whose support is wider than the range of the parameter, e.g. a
normal prior for a parameter bounded below by 0, is truncated to that
range and renormalised; the normalisation does not change the estimates,
but enters the marginal likelihood of
[`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md).

## See also

[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md),
[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md),
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md),
[`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md)

## Examples

``` r
normal(8, 1)
#> normal(8, 1) 
pert(0.5, 3, 10)
#> pert(0.5, 3, 10) 
half_normal(1)
#> half_normal(1) 
# posterior modes with priors on the upper asymptote, the ED50 and the residual standard deviation
pmod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
               start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE),
               prior = list(b2 = normal(8, 1), b4 = pert(0.5, 3, 10), sigma = half_normal(1)))
summary(pmod)
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, start = c(b1 = 0, 
#>     b5 = 0), fix = c(b1 = TRUE, b5 = TRUE), model = "loglogistic", 
#>     prior = list(b2 = normal(8, 1), b4 = pert(0.5, 3, 10), sigma = half_normal(1)))
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 24
#> Log-likelihood: -18.61 (df = 4), AIC: 45.22, BIC: 49.93
#> 
#> Priors:
#>   b2 ~ normal(8, 1)
#>   b4 ~ pert(0.5, 3, 10)
#>   sigma ~ half_normal(1)
#> 
#> Coefficients (posterior modes):
#>                Estimate Std. Error
#> b2:(Intercept)    7.862      0.191
#> b3:(Intercept)   -2.467      0.323
#> b4:(Intercept)    3.259      0.185
#> Fixed: b1:(Intercept) = 0, b5:(Intercept) = 0
#> 
#> Residual standard deviation: 0.5333
#> 
#> Optimiser: nlminb: both X-convergence and relative convergence (5), maximum absolute gradient 3.4e-07
```
