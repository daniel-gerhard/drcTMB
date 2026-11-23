# Bayesian model-averaged benchmark dose

The posterior distribution of the benchmark dose averaged over
dose-response models fitted in the BMD parameterisation with priors, as
in the Bayesian benchmark dose modelling of EFSA (2022): the models are
weighted by their posterior probabilities, from their marginal
likelihoods
([`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md)),
and the BMD and its credible interval are the median and quantiles of
the mixture of the posterior distributions of the BMD of the models.

## Usage

``` r
bmdbayes(
  ...,
  prior.weights = NULL,
  method = c("laplace", "is"),
  n = 5000,
  df = 4
)

# S3 method for class 'bmdbayes'
coef(object, ...)

# S3 method for class 'bmdbayes'
confint(object, parm = NULL, level = 0.9, ...)

bmdbayesmodels(object, level = 0.9)

# S3 method for class 'bmdbayes'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- ...:

  drmTMB objects fitted in the BMD parameterisation (`bmdpar` of
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md))
  with the same benchmark response and definition, the same prior for
  the BMD (`prior=list(bmd=...)`) and proper priors for all estimated
  parameters, to the same doses, responses and weights; the models and
  families may differ. A single model is allowed. Models are labelled by
  their argument names (e.g. LL = fit1) or, if unnamed, by the
  expressions. The fits use
  `drmTMB(..., bmdpar=list(bmr=, def=), prior=list(bmd=..., ...))`, in
  which b4 is replaced by the log BMD (coefficient log_bmd), so that the
  BMD has a prior;
  [`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md)
  gives the default priors of EFSA, and
  [`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)
  does all steps.

- prior.weights:

  prior probabilities of the models, in the order of `...`; equal by
  default.

- method:

  "laplace" (default) for the Laplace approximation of the marginal
  likelihoods and a normal posterior distribution of the log BMD of each
  model, or "is" for importance sampling of both (see
  [`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md));
  importance sampling is recommended when posterior modes are close to
  the bounds of PERT priors

- n, df:

  number of draws of the final round and degrees of freedom of the
  proposal of importance sampling (see
  [`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md))

- parm:

  not used; there is one BMD

- level:

  the probability of the credible interval; 0.90 by default, for the
  BMDL and BMDU

- x, object:

  a bmdbayes object

- digits:

  number of significant digits in the printed output

## Value

bmdbayes returns a bmdbayes object, a list with the model-averaged BMD
(`estimates`), the table of the models (`models`: model, family, log
marginal likelihood, prior and posterior probabilities, and with
importance sampling the effective sample size), the posterior
distribution of the log BMD of each model (`posterior`: the mean and
standard deviation of the normal distribution, or the draws and their
weights), the lowest dose above zero, and the arguments. coef returns
the BMD, confint the credible interval, and `bmdbayesmodels` the BMDs
and credible intervals of the models, a data frame with the posterior
median BMD, BMDL and BMDU of each model.

## Details

The posterior probability of model k is proportional to its prior
probability times its marginal likelihood p(y \| model k). With
`method="laplace"`, the posterior distribution of the log BMD of each
model is normal with the posterior mode of log_bmd as the mean and its
standard error (the curvature of the log posterior) as the standard
deviation, so that the median BMD of a model is its posterior mode; the
model-averaged posterior distribution of the log BMD is the mixture of
these normal distributions with the posterior probabilities of the
models as weights. With `method="is"`, the importance sampling draws of
each model, weighted by their importance weights and the posterior
probability of the model, give both the marginal likelihoods and the
mixture, without the normal approximation. The BMD is the median of the
mixture, and the BMDL and BMDU the quantiles (1 - level)/2 and (1 +
level)/2, the 5% and 95% quantiles of a 90% credible interval.

With random effects, the BMD of the BMD parameterisation is that of the
curve of the intercepts, with the random effects at zero (the
conditional BMD of
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)).

A model whose marginal likelihood cannot be computed (e.g. the fit did
not converge, or its Hessian is not positive definite) gets the weight
0, with a warning. `print` notes the criteria of EFSA (2022) for a
credible interval that is too wide, with the 90% interval: BMD/BMDL
above 20, BMDU/BMDL above 50, and a BMDL below a tenth of the lowest
dose above zero.

## References

EFSA Scientific Committee (2022). Guidance on the use of the benchmark
dose approach in risk assessment. EFSA Journal 20(10), 7584.

## See also

[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md),
[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md),
[`marglik`](https://daniel-gerhard.github.io/drcTMB/reference/marglik.md),
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md),
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md),
[`bmdma`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)

## Examples

``` r
pr <- list(bmd = pert(0.1, 1, 5), b2 = normal(8, 1), b3 = normal(-2, 2), sigma = half_normal(1))
fit <- function(model)
  drmTMB(rootl ~ conc, data = ryegrass, model = model,
         start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE),
         bmdpar = list(bmr = 0.05, def = "relative"), prior = pr)
(ba <- bmdbayes(LL = fit("loglogistic"), W1 = fit("weibull1"), W2 = fit("weibull2")))
#> Bayesian model-averaged benchmark dose, relative deviation 0.05, prior pert(0.1, 1, 5) for the BMD
#> Laplace approximation
#> 
#>          model   family  logML  prior   weight    BMD   BMDL  BMDU
#> LL loglogistic gaussian -27.05 0.3333 0.215716 1.0102 0.7724 1.321
#> W1    weibull1 gaussian -30.57 0.3333 0.006396 0.8108 0.5484 1.199
#> W2    weibull2 gaussian -25.77 0.3333 0.777888 1.3590 1.1538 1.601
#> 
#>                 BMD   BMDL  BMDU
#> model average 1.316 0.8869 1.582
#> 
#> Posterior median and 90% credible interval
coef(ba)
#>      BMD 
#> 1.315526 
confint(ba)
#>           5 %    95 %
#> BMD 0.8869296 1.58189
```
