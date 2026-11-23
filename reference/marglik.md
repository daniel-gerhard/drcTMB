# Marginal likelihood of a drmTMB model with priors

The log marginal likelihood log p(y), the likelihood integrated over the
prior of all parameters, for Bayesian model comparison and model
averaging: the posterior probability of a model is proportional to its
prior probability times its marginal likelihood.

## Usage

``` r
marglik(object, method = c("laplace", "is"), n = 5000, df = 4)
```

## Arguments

- object:

  a drmTMB object, fitted with priors for all parameters

- method:

  "laplace" (default) or "is" (importance sampling)

- n:

  number of draws of the final round of importance sampling

- df:

  degrees of freedom of the multivariate t proposal of importance
  sampling

## Value

the log marginal likelihood, a number of class "drcmarglik" with the
attributes method, the number of parameters k, and for importance
sampling n, the effective sample size ess, the Monte Carlo standard
error se, the number of draws whose objective could not be evaluated
(failed, given weight 0), and the Laplace approximation laplace.

## Details

The fit needs a proper prior for every estimated parameter (see
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)):
all free coefficients, the standard deviations of the random effects,
the residual parameter (sigma, phi or theta) and the parameters of a
residual correlation structure or variance function (see
[`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md));
with a flat prior, the marginal likelihood is not defined. Correlations
of random effects have no priors, so that random effects of several
parameters or with slopes need independent terms (`||` or separate
formulas). The random effects are integrated out by the Laplace
approximation of the fit.

`method="laplace"` approximates the log posterior of the parameters (on
the scale of their optimisation, with the log Jacobians of the priors)
by a normal distribution at the posterior mode: log p(y) = -objective +
k/2 log(2 pi) - 1/2 log det H, with k parameters and the Hessian H of
the objective, as in the Bayesian benchmark dose modelling of EFSA
(2022). The random effects are integrated out first, so that this is a
nested Laplace approximation. It is poor if the posterior mode is close
to a bound, e.g. of a pert prior or a multistage coefficient at 0, which
gives a warning. For small dose-response data sets, the Laplace
approximation can be noticeably biased, in particular with posterior
modes close to the bounds of pert priors, which changes the posterior
model probabilities by a few percentage points; check it with
`method="is"`.

`method="is"` estimates the marginal likelihood by adaptive importance
sampling. The proposal is a multivariate t distribution with `df`
degrees of freedom, first centred at the posterior mode with the
covariance of the fit as its scale matrix; three rounds of n/2 draws (at
least 10 per parameter) move it to the weighted mean and covariance of
the draws, which follows skewed or curved posteriors (e.g. a strong
correlation of the log BMD with the slope), and a final round of n draws
gives the estimate. The proposal of each round is a mixture, with
probability 0.2 of the initial proposal, so that the weights stay
bounded in its tails. Draws outside the bounds of the parameters have
weight 0. The result reports the effective sample size of the weights of
the final round and the Monte Carlo standard error of the log marginal
likelihood; a small effective sample size (e.g. below 100) means the
estimate is unreliable. With random effects, each draw needs an inner
optimisation, so that importance sampling is slower.

## See also

[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md),
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)

## Examples

``` r
# proper priors for all estimated parameters, in the BMD parameterisation
pr <- list(bmd = pert(0.1, 1, 5), b2 = normal(8, 1), b3 = normal(-2, 2), sigma = half_normal(1))
fit <- function(model)
  drmTMB(rootl ~ conc, data = ryegrass, model = model,
         start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE),
         bmdpar = list(bmr = 0.05, def = "relative"), prior = pr)
LL <- fit("loglogistic")
W2 <- fit("weibull2")
marglik(LL)
#> Log marginal likelihood: -27.05 
#> Laplace approximation, 4 parameters
marglik(W2)
#> Log marginal likelihood: -25.77 
#> Laplace approximation, 4 parameters
# importance sampling as a check of the Laplace approximation
set.seed(1)
marglik(LL, method = "is", n = 2000)
#> Log marginal likelihood: -26.83 
#> Importance sampling, 2000 draws, effective sample size 1174
#> Monte Carlo standard error 0.019; Laplace approximation -27.05
```
