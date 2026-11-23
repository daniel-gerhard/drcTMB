# A PERT prior from a posterior distribution

Fits a PERT distribution to the posterior distribution of a parameter,
e.g. the model-averaged BMD of
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
as an informative prior for the analysis of a new study (Kremer et al.
2026, EFSA EN-9902): the posterior is summarised by 1000 quantiles; the
minimum and maximum of the PERT are their smallest and largest value
(the 0.05% and 99.95% quantiles), so that the prior excludes values
outside them; the mode is the mode of their kernel density, and the
shape is estimated by maximum likelihood (between 0 and 20). Several
posteriors, e.g. of several historical studies, are mixed with equal
weights, the recommended way to combine studies. With discounting, the
shape is reduced so that the variance of the PERT doubles
(approximately, as the shape is at least 1e-4).

## Usage

``` r
fitpert(..., weights = NULL, discount = FALSE)
```

## Arguments

- ...:

  bmdbayes objects (the posterior of the BMD), or numeric vectors of
  posterior draws

- weights:

  for a single numeric vector of draws, their weights (e.g. of
  importance sampling)

- discount:

  logical, discount the prior (double its variance)

## Value

a prior (see
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)),
e.g. for `prior=list(bmd=...)` of `drmTMB` or the bmd argument of
[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md)

## References

Kremer, C., Saenen, N. D., Bijnens, K., Faes, C., De Boever, P., Smeets,
K. and Aerts, M. (2026). Compilation of NTP and published dose-response
data: study of informative priors in the EFSA platform for Bayesian
benchmark dose analysis. EFSA Supporting Publications, EN-9902.

## See also

[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md),
[`bmdbayes`](https://daniel-gerhard.github.io/drcTMB/reference/bmdbayes.md),
[`efsapriors`](https://daniel-gerhard.github.io/drcTMB/reference/efsapriors.md),
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)

## Examples

``` r
# a PERT prior from posterior draws of a BMD
set.seed(1)
draws <- exp(rnorm(2000, log(2), 0.3))
fitpert(draws)
#> pert(0.762757, 1.83879, 5.95985, shape=6.20185) 

# from the model-averaged posterior of the BMD of bmdbayes (or bmdefsa), e.g. of a historical
# study, as an informative prior of the BMD for a current study
pr <- list(bmd = pert(0.1, 1, 5), b2 = normal(8, 1), b3 = normal(-2, 2), sigma = half_normal(1))
fit <- function(model)
  drmTMB(rootl ~ conc, data = ryegrass, model = model,
         start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE),
         bmdpar = list(bmr = 0.05, def = "relative"), prior = pr)
ba <- bmdbayes(LL = fit("loglogistic"), W2 = fit("weibull2"))
fitpert(ba)
#> pert(0.636195, 1.34061, 1.87419, shape=4.20974) 
```
