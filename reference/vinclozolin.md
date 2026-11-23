# Vinclozolin in an androgen receptor reporter gene assay

Luminescence of Chinese hamster ovary cells exposed to concentrations of
vinclozolin for 22 hours; the experiment was repeated six times on
different days.

## Usage

``` r
vinclozolin
```

## Format

A data frame with 53 rows and 3 variables:

- exper:

  experiment (day), a factor with 6 levels

- conc:

  concentration of vinclozolin

- effect:

  luminescence

## Source

Nellemann, C., Dalgaard, M., Lam, H. R. and Vinggaard, A. M. (2003) The
combined effects of vinclozolin and procymidone do not deviate from
expected additivity in vitro and in vivo. *Toxicological Sciences*,
**71**, 251–262; from the package drcData.

## Examples

``` r
vmod <- drmTMB(effect ~ conc, data = vinclozolin, model = "loglogistic",
               rform = b1 + b2 ~ 1 | exper, start = c(b5 = 0), fix = c(b5 = TRUE),
               lower = c(0, 0, -Inf, 0, -Inf))  # lower bounds of b1-b5
# ED50 of the population curve, with profile likelihood interval
confint(ed(vmod, 50, random = FALSE))
#> Computing profile confidence intervals ...
#>             2.5 %    97.5 %
#> 50:lf1 0.03225858 0.0922736
```
