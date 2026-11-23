# Bayesian model-averaged BMD with the EFSA defaults

Fits the models of the Bayesian benchmark dose modelling of EFSA (2022)
available in drcTMB, in the BMD parameterisation with the priors of
[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md),
for each model and family, and averages them by
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md).

## Usage

``` r
bmdefsa(
  form,
  data,
  bmr = NULL,
  def = NULL,
  models = NULL,
  family = c("gaussian", "lognormal"),
  ...,
  prior.weights = NULL,
  method = c("is", "laplace"),
  n = 5000,
  df = 4,
  control = list()
)
```

## Arguments

- form, data:

  formula and data, as in
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)

- bmr, def:

  the benchmark response and its definition (see
  [`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md))

- models:

  the models (see
  [`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md)),
  by default all: six for continuous responses, and eight, including
  probit and logit, for quantal responses

- family:

  the families: for continuous responses "gaussian" and "lognormal"
  (both by default), for quantal responses "binomial"

- ...:

  further arguments of
  [`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md):
  background, maxresponse, bmd, dprior, extended

- prior.weights, method, n, df:

  see
  [`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md);
  prior weights in the order of the fits, for each model the families,
  labelled model_N (gaussian), model_LN (lognormal) and model_Q
  (binomial), e.g. weibull2_N, weibull2_LN, weibull1_N, .... The default
  method is importance sampling ("is"): with the PERT priors, the
  posterior modes are often close to a bound of a prior, e.g. a
  background close to 0, where the Laplace approximation is poor.

- control:

  a list of control parameters of
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)

## Value

a bmdbayes object, with the fits (`fits`) and the priors (`efsa`) of the
models; models whose fit fails are left out with a warning

## Details

The data requirements are those of
[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md).

## References

EFSA Scientific Committee (2022). Guidance on the use of the benchmark
dose approach in risk assessment. EFSA Journal 20(10), 7584.

## See also

[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md),
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md),
[`TCDD`](https://daniel-gerhard.github.io/drcTMB/reference/TCDD.md),
[`liver.tumor`](https://daniel-gerhard.github.io/drcTMB/reference/liver.tumor.md)

## Examples

``` r
# \donttest{
# quantal: liver tumours after exposure to TCDD
set.seed(1)
(ba <- bmdefsa(cbind(incidence, total - incidence) ~ conc, data = TCDD, family = "binomial"))
#> Bayesian model-averaged benchmark dose, excess risk 0.1, prior pert(0, 19.28, 115.68, shape=0.0001) for the BMD
#> Importance sampling
#> 
#>                     model   family  logML prior  weight  ess    BMD  BMDL  BMDU
#> weibull2_Q       weibull2 binomial -13.24 0.125 0.13971 3045  9.517 4.004 20.93
#> weibull1_Q       weibull1 binomial -13.84 0.125 0.07622 3679  7.835 3.969 15.79
#> loglogistic_Q loglogistic binomial -13.05 0.125 0.16843 3454  8.948 3.873 19.70
#> lognormal_Q     lognormal binomial -13.74 0.125 0.08462 3323  8.628 4.022 18.80
#> gamma_Q             gamma binomial -13.17 0.125 0.14976 2799  9.947 4.204 20.99
#> multistage_Q   multistage binomial -12.51 0.125 0.28948 4044 12.766 8.055 21.82
#> probit_Q           probit binomial -14.22 0.125 0.05226 1553 12.198 5.272 26.92
#> logit_Q             logit binomial -14.50 0.125 0.03953 1832 12.965 5.363 26.20
#> 
#>                 BMD  BMDL  BMDU
#> model average 10.56 4.522 21.56
#> 
#> Posterior median and 90% credible interval
# }
```
