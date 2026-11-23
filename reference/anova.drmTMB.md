# Likelihood ratio tests for drmTMB models

Compares two or more models fitted to the same data by likelihood ratio
tests, e.g. a model with and without covariates in `fform`, with and
without random effects, or a poisson and a negative binomial model. A
model with `model="constant"` as the smaller model tests for a dose
effect.

## Usage

``` r
# S3 method for class 'drmTMB'
anova(
  object,
  ...,
  test = c("Chisq", "F", "boot"),
  nsim = 1000,
  seed = NULL,
  progress = interactive()
)
```

## Arguments

- object:

  a drmTMB object

- ...:

  further drmTMB objects

- test:

  "Chisq" (default) for likelihood ratio tests, "F" for F-tests of
  gaussian or lognormal models without random effects, priors,
  correlation structure and variance function, or "boot" for parametric
  bootstrap likelihood ratio tests

- nsim:

  number of simulations for `test="boot"`

- seed:

  an optional seed for the random number generator, for `test="boot"`

- progress:

  logical, show a progress bar for `test="boot"`

## Value

a data frame of class anova with the number of parameters (npar), AIC,
BIC, log-likelihood and deviance (minus twice the log-likelihood) of
each model, and the likelihood ratio statistic (Chisq), its degrees of
freedom (Df) and p-value (Pr(\>Chisq)) for the comparison with the
previous model. For `test="F"`, the residual degrees of freedom (Res.Df)
and residual sum of squares (RSS) instead of the deviance, and the F
statistic (F) and its p-value (Pr(\>F)). For `test="boot"`, additionally
the bootstrap p-value (Pr(boot)), with the simulated likelihood ratio
statistics of each comparison in the attribute `"boot"`.

## Details

The models are ordered by their number of parameters, and each model is
compared with the previous one: the likelihood ratio statistic is twice
the difference of their log-likelihoods (the Laplace approximation of
the marginal log-likelihoods for models with random effects), with the
difference in the number of parameters as degrees of freedom. The test
assumes that the smaller model is nested in the larger one, which is not
checked. For parameters on the boundary of their range under the smaller
model, e.g. a standard deviation of random effects of zero, or an
infinite size of the negative binomial for the poisson model, the
chi-square distribution is conservative: the p-value is approximately
twice the correct value for a single standard deviation.

For gaussian and lognormal (on the log scale) models without random
effects, priors, correlation structure and variance function, `test="F"`
gives F-tests, as drc does for continuous responses: the reduction of
the (weighted) residual sum of squares per additional parameter, divided
by the residual mean square of the largest model, with n - p residual
degrees of freedom for p coefficients b1-b5. The F-test is exact for
linear models and accounts for the estimation of the residual variance,
which the likelihood ratio test does not, so it is usually preferable
for small experiments.

`test="boot"` gives parametric bootstrap likelihood ratio tests, as
`PBmodcomp` of the package pbkrtest: for each comparison, `nsim`
responses are simulated from the smaller model, with new random effects,
both models are fitted to each of them, and the p-value `Pr(boot)` is
the proportion of the simulated likelihood ratio statistics at least as
large as the observed one, (1 + \#{LR\* \>= LR}) / (1 + nsim). It does
not rely on the chi-square distribution, which is only approximate in
small samples, for parameters on the boundary (standard deviations of
random effects of zero, the poisson model within the negative binomial),
and for parameters that are not identified under the smaller model (e.g.
the steepness and location of the curve compared with
`model="constant"`). The chi-square p-value is shown for comparison.
Each model is fitted twice to each simulated response, starting at its
estimates and at the starting values of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md),
and the better fit is used; simulations without a finite log-likelihood
of both models are left out.

## References

Halekoh, U. and Hojsgaard, S. (2014). A Kenward-Roger approximation and
parametric bootstrap methods for tests in linear mixed models - the R
package pbkrtest. Journal of Statistical Software 59(9), 1-32.

## See also

[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)
(model="constant" for a test of a dose effect),
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
# a free lower asymptote b1
mod4 <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
               start = c(b5 = 0), fix = c(b5 = TRUE))
anova(mod, mod4)
#> Likelihood ratio tests of drmTMB models
#> 
#>      npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)  
#> mod     4 45.208 49.920 -18.604   37.208                      
#> mod4    5 42.310 48.201 -16.155   32.310 4.898  1    0.02689 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(mod, mod4, test = "F")
#> F-tests of drmTMB models
#> 
#>      npar    AIC    BIC  logLik Res.Df    RSS     F Df  Pr(>F)  
#> mod     4 45.208 49.920 -18.604     21 6.6228                   
#> mod4    5 42.310 48.201 -16.155     20 5.4002 4.528  1 0.04597 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# test for any dose effect, against a constant mean
mod0 <- drmTMB(rootl ~ conc, data = ryegrass, model = "constant")
anova(mod0, mod)
#> Likelihood ratio tests of drmTMB models
#> 
#>      npar     AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
#> mod0    2 128.188 130.54 -62.094  124.188                        
#> mod     4  45.208  49.92 -18.604   37.208 86.98  2  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
