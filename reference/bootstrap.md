# Parametric bootstrap of a drmTMB model

Simulates responses from the fitted model, with new random effects (see
[`simulate.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)),
and fits the model again to each simulated response, with the bounds,
priors and optimiser settings of the fit (`control` in
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)),
starting at the estimates. The bootstrap sample can be passed as the
argument boot to
[`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
(for the coefficients, and for the results of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md) and
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md)) and
to
[`predict.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
so that the model is not fitted again.

## Usage

``` r
bootstrap(object, nsim = 500, seed = NULL, progress = interactive())

# S3 method for class 'bootstrap.drmTMB'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'bootstrap.drmTMB'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  a drmTMB object

- nsim:

  number of simulations

- seed:

  an optional seed for the random number generator

- progress:

  logical, show a progress bar

- x:

  a bootstrap.drmTMB object

- digits:

  number of significant digits in the printed output

- ...:

  further arguments, not used

- parm:

  parameters, names or positions in `coef`; default all

- level:

  the confidence level

## Value

bootstrap returns a bootstrap.drmTMB object, a list with the estimates
of `coef(object)` for each simulation (`estimates`, a matrix with NA
rows for fits that did not converge), the fitted model (`fit`), and
`nsim` and `seed`. confint returns percentile intervals of the
parameters, on the scale of
[`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md).
print shows the estimate, the mean of the bootstrap estimates, the bias
and the bootstrap standard error of each parameter.

## See also

[`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md),
[`predict.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
b <- bootstrap(mod, nsim = 100, seed = 1)
b
#> Parametric bootstrap: 100 of 100 fits converged
#> 
#>                Estimate   Mean     Bias Std. Error
#> b2:(Intercept)   7.8554  7.829 -0.02618     0.1838
#> b3:(Intercept)  -2.4703 -2.544 -0.07400     0.3279
#> b4:(Intercept)   3.2634  3.280  0.01673     0.1806
#> log_sigma       -0.6438 -0.692 -0.04821     0.1448
confint(b)
#>                     2.5 %     97.5 %
#> b2:(Intercept)  7.5016059  8.1870409
#> b3:(Intercept) -3.1500315 -1.9681389
#> b4:(Intercept)  2.9866350  3.6290124
#> sigma           0.3575887  0.6327268
predict(mod, newdata = data.frame(conc = c(1, 3, 10)), interval = "confidence", boot = b)
#>         fit    se.fit       lwr       upr
#> 1 7.4541188 0.1691204 7.1312504 7.7435829
#> 2 4.3344621 0.2503276 3.9180362 4.7813573
#> 3 0.4648086 0.1604718 0.1921772 0.7888688
```
