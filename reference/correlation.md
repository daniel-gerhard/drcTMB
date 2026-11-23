# Residual correlation structures and variance functions

Correlated residuals within series (e.g. repeated measurements of a
subject over time, or plots in a field) and residual standard deviations
that depend on a covariate or on the fitted mean, for the gaussian and
lognormal families of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)
(arguments `correlation` and `variance`), as the corStruct and varFunc
classes of nlme.

## Usage

``` r
cor_ar1(form = ~1)

cor_car1(form)

cor_cs(form = ~1)

cor_exp(form, nugget = FALSE)

cor_gaus(form, nugget = FALSE)

cor_spher(form, nugget = FALSE)

var_power(form = ~fitted)

var_exp(form = ~fitted)

var_ident(form)
```

## Arguments

- form:

  a one-sided formula. For the correlation structures, ~ position \|
  series; cor_ar1 and cor_cs default to ~ 1, the order of the rows, with
  all observations one series unless a series is given, e.g. ~ 1 \|
  subject. For var_power and var_exp, ~ fitted (default, the fitted
  mean) or a covariate, e.g. ~ dose; for var_ident, ~ 1 \| stratum
  (required).

- nugget:

  logical; if TRUE, a nugget parameter in (0, 1) reduces the correlation
  at positive distances, see Details

## Value

a correlation structure (class drccor) or a variance function (class
drcvar)

## Details

The residuals e (y minus the mean; for the lognormal family log(y) minus
the log median) of a series are multivariate normal with covariance D R
D: R is the correlation matrix of the structure, and D the diagonal
matrix of the residual standard deviations, sigma times the variance
function g. Series are independent, and observations without a series
are one series.

The formula of a correlation structure is `~ position | series`: the
position of the observations in time or space (one or more numeric
variables; Euclidean distances), and an optional grouping factor of the
series (several factors with `a:b` or `a/b`). The correlation of two
residuals of a series at distance d:

- `cor_ar1`: rho^k for observations k positions apart, -1 \< rho \< 1;
  the position is an integer time (or, with `~ 1 | series`, the order of
  the observations in the data), as corAR1.

- `cor_car1`: phi^d for times d apart, 0 \< phi \< 1 (continuous-time
  AR1, corCAR1).

- `cor_cs`: rho for all pairs (compound symmetry, corCompSymm),
  -1/(m - 1) \< rho \< 1 for series of up to m observations; for rho \>
  0 the same as a random intercept of the series.

- `cor_exp`: exp(-d/range) (corExp); `cor_gaus`: exp(-(d/range)^2)
  (corGaus); `cor_spher`: 1 - 1.5 d/range + 0.5 (d/range)^3 for d \<
  range, 0 beyond (corSpher). With `nugget=TRUE`, the correlation is
  (1 - nugget) times that, for d \> 0.

Positions must be distinct within a series.

Variance functions give the residual standard deviation sigma g(v) for a
covariate v, the fitted mean (the median for lognormal, `~ fitted`, the
default) or a variable of the data (e.g. `~ dose`):

- `var_power`: g = \|v\|^delta (varPower).

- `var_exp`: g = exp(delta v) (varExp).

- `var_ident`: g = exp(lambda_s) for the strata s of `~ 1 | stratum`,
  with lambda 0 for the first stratum (varIdent): a separate standard
  deviation for each stratum.

With the fitted mean as the covariate, the variance depends on the
coefficients of the curve, and the likelihood is maximised over both,
unlike the iterative estimation of gnls in nlme, whose estimates differ
slightly.

The parameters are estimated on an unconstrained scale (`coef` of the
fit): for the correlation, cor_par:rho (atanh(rho) for AR1, the logit of
(rho + 1/(m - 1))/(1 + 1/(m - 1)) for compound symmetry), cor_par:phi
(logit(phi)), cor_par:range (log(range)) and cor_par:nugget
(logit(nugget)); for the variance functions var_par:delta (delta itself)
and `var_par:STRATUM` (lambda, the log ratio of the standard deviation
of the stratum to that of the first). `summary` prints them on their
natural scale, and `confint` gives their intervals on that scale, named
cor:rho, cor:phi, cor:range, cor:nugget, var:delta and
`var:ratio:STRATUM`, with the name of the stratum, e.g. var:ratio:B.
Priors (argument `prior` of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md))
are given by these names, on the natural scale, e.g.
`prior=list("cor:phi"=pert(0, 0.5, 1))`; the marginal likelihood
([`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md))
needs priors for all of them.

## References

Pinheiro, J. C. and Bates, D. M. (2000). Mixed-Effects Models in S and
S-PLUS. Springer, New York.

## See also

[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md),
[`corClasses`](https://rdrr.io/pkg/nlme/man/corClasses.html),
[`varClasses`](https://rdrr.io/pkg/nlme/man/varClasses.html)

## Examples

``` r
# standard deviation a power of the fitted mean, with random effects of the experiments
vmod <- drmTMB(effect ~ conc, data = vinclozolin, model = "loglogistic",
               rform = b1 + b2 ~ 1 | exper, start = c(b5 = 0), fix = c(b5 = TRUE),
               lower = c(0, 0, -Inf, 0, -Inf), variance = var_power())
summary(vmod)
#> Call:
#> drmTMB(form = effect ~ conc, rform = b1 + b2 ~ 1 | exper, data = vinclozolin, 
#>     start = c(b5 = 0), fix = c(b5 = TRUE), lower = c(0, 0, -Inf, 
#>         0, -Inf), model = "loglogistic", variance = var_power())
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 53, groups: exper 6
#> Log-likelihood: -335.9 (df = 9), AIC: 689.9, BIC: 707.6
#> 
#> Coefficients:
#>                  Estimate Std. Error
#> b1:(Intercept)  252.72312     78.637
#> b2:(Intercept) 1924.39964    258.759
#> b3:(Intercept)   -0.72717      0.103
#> b4:(Intercept)    0.05837      0.013
#> Fixed: b5:(Intercept) = 0
#> 
#> Random effects:
#>  Groups Name Std.Dev. Corr
#>  exper  b1   156.8        
#>         b2   566.6    0.96
#> 
#> Residual standard deviation: 0.2197
#> 
#> Variance function: power, ~fitted: delta = 0.9453
#> 
#> Optimiser: nlminb: relative convergence (4), maximum absolute gradient 0.00017
# AR(1) correlation of neighbouring concentrations within each experiment
vmodc <- update(vmod, correlation = cor_ar1(~ 1 | exper))
anova(vmod, vmodc)
#> Likelihood ratio tests of drmTMB models
#> 
#>       npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#> vmod     9 689.86 707.59 -335.93   671.86                     
#> vmodc   10 691.72 711.42 -335.86   671.72 0.1422  1     0.7061
# a separate standard deviation for each experiment
update(vmod, variance = var_ident(~ 1 | exper))
#> Warning: NaNs produced
#> Call:
#> drmTMB(form = effect ~ conc, rform = b1 + b2 ~ 1 | exper, data = vinclozolin, 
#>     start = c(b5 = 0), fix = c(b5 = TRUE), lower = c(0, 0, -Inf, 
#>         0, -Inf), model = "loglogistic", variance = var_ident(~1 | 
#>         exper))
#> 
#> Model:  loglogistic
#> Family: gaussian, identity link for the asymptotes
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  b4:(Intercept)  
#>      246.50255      1857.18772        -0.70407         0.07133  
#> Fixed: b5:(Intercept) = 0
#> 
#> Random effects:
#>  Groups Name Std.Dev. Corr
#>  exper  b1   133.3        
#>         b2   618.9    1.00
#> 
#> Residual standard deviation: 142.5
#> 
#> Variance function: separate standard deviations, ~1 | exper: ratio:10821 = 0.6035, ratio:10828 = 0.8605, ratio:10904 = 1.5230, ratio:11023 = 0.4860, ratio:11106 = 1.6844
```
