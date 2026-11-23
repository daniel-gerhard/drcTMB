# Robust (sandwich) covariance of a drmTMB fit

Robust (sandwich) standard errors for drmTMB fits, robust to a
misspecified variance and, with clusters, to correlation within the
clusters. vcovrobust is passed as the argument vcov. of
[`summary.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
[`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
(Wald intervals),
[`predict.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md) and
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md).

## Usage

``` r
vcovrobust(object, cluster = NULL, type = c("HC0", "HC1"))

# S3 method for class 'drmTMB'
estfun(x, cluster = NULL, ...)

# S3 method for class 'drmTMB'
bread(x, cluster = NULL, ...)
```

## Arguments

- object, x:

  a drmTMB object

- cluster:

  the clusters: a one-sided formula, e.g. ~ plate, evaluated in the data
  of the fit, or a vector with an entry for each observation of the fit
  (the rows of the data without missing values). NULL for the default
  clusters, see Details.

- type:

  "HC0" (default) without, or "HC1" with the small-sample adjustment.

- ...:

  anything else

## Value

vcovrobust and bread return a covariance matrix of the estimated
parameters, named as `vcov`; estfun returns a matrix of scores with a
row for each cluster (or observation) and a column for each parameter.

## Details

vcovrobust computes the sandwich covariance V M V of all estimated
parameters (the coefficients and the variance parameters, as `vcov`),
with the inverse Hessian V of the fit and the outer products M of the
scores of independent clusters, the gradients of their contributions to
the log-likelihood. It is robust to a misspecified variance (e.g.
heteroscedasticity, overdispersion) and, with clusters, to correlation
within the clusters. The result can be given as the argument vcov. of
[`summary.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
[`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
(Wald intervals),
[`predict.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md) and
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md), e.g.
`ed(fit, 50, vcov.=vcovrobust)` or
`ed(fit, 50, vcov.=function(f) vcovrobust(f, cluster=~plate))`.

Without random effects, the scores of the observations are the gradients
of their contributions to the log-likelihood, by automatic
differentiation. With random effects, the clusters must contain whole
groups: the score of a cluster is the gradient of the Laplace
approximation of its marginal log-likelihood, from a fit of the model
with the weights of all other observations set to zero (as in glmmTMB).
The groups of all random effect terms, the series of a correlation
structure
([`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md))
and the series of event data must be nested in the clusters; by default,
the clusters are the groups of the random effect term that all others
are nested in (e.g. lab for 1 \| lab/plate), or the series, and
otherwise the observations.

The sandwich covariance is not available for fits with priors, whose
prior densities are not contributions of independent clusters. Small
numbers of clusters give too small standard errors; type="HC1"
multiplies the covariance by G/(G - 1) (n - 1)/(n - p), for G clusters,
n observations with positive weights and p estimated parameters, as
type="HC1" of
[`sandwich::vcovCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html).

estfun and bread methods for the generics of package sandwich are
provided, with the scores of the default clusters as estfun and n V as
bread for n clusters, so that `sandwich::sandwich(fit)` equals
`vcovrobust(fit)`, and
[`sandwich::vcovCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
can aggregate the scores of the observations of a fit without random
effects or correlation structures. With random effects, estfun has a row
for each cluster.

## References

Zeileis, A. (2006). Object-oriented computation of sandwich estimators.
Journal of Statistical Software 16(9), 1-16.

Zeileis, A., Koell, S. and Graham, N. (2020). Various versatile
variances: an object-oriented implementation of clustered covariances in
R. Journal of Statistical Software 95(1), 1-36.

## See also

[`vcovCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)

## Examples

``` r
fit <- drmTMB(rootl ~ conc, data=ryegrass, model="loglogistic")
vcovrobust(fit)
#>                b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept)
#> b1:(Intercept)    0.021658774   -0.001528977     0.20600311     0.08653237
#> b2:(Intercept)   -0.001528977    0.023866632     0.07467818     0.01120339
#> b3:(Intercept)    0.206003110    0.074678183     5.14137568     1.82341424
#> b4:(Intercept)    0.086532369    0.011203386     1.82341424     0.69364723
#> b5:(Intercept)    0.099749710    0.021308337     2.08504735     0.76915311
#> log_sigma         0.002460853    0.001430827     0.08405124     0.02973294
#>                b5:(Intercept)   log_sigma
#> b1:(Intercept)     0.09974971 0.002460853
#> b2:(Intercept)     0.02130834 0.001430827
#> b3:(Intercept)     2.08504735 0.084051242
#> b4:(Intercept)     0.76915311 0.029732941
#> b5:(Intercept)     0.87340946 0.032708796
#> log_sigma          0.03270880 0.030198062
summary(fit, vcov.=vcovrobust)
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, model = "loglogistic")
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 24
#> Log-likelihood: -15.88 (df = 6), AIC: 43.76, BIC: 50.82
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b1:(Intercept)   0.3215      0.147
#> b2:(Intercept)   7.7605      0.154
#> b3:(Intercept)  -3.9301      2.267
#> b4:(Intercept)   2.2148      0.833
#> b5:(Intercept)  -0.7594      0.935
#> Standard errors from the covariance vcov.
#> The asymmetry parameter b5 is on the log scale.
#> 
#> Residual standard deviation: 0.527 on 19 degrees of freedom
#> 
#> Optimiser: nlminb: relative convergence (4), maximum absolute gradient 3.1e-12
ed(fit, c(10, 50), vcov.=vcovrobust)
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 10:lf1     1.56     0.186
#> 50:lf1     3.02     0.331
```
