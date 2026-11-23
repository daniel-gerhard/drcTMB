# Immobility of daphnids

The number of immobile daphnids out of 20, for several concentrations of
a toxic substance; the same daphnids were counted after 24 and after 48
hours.

## Usage

``` r
daphnids
```

## Format

A data frame with 16 rows and 4 variables:

- dose:

  concentration of the toxic substance

- no:

  number of immobile daphnids

- total:

  total number of daphnids

- time:

  time of the count, a factor with levels `24h` and `48h`

## Source

Nina Cedergreen, Faculty of Life Sciences, University of Copenhagen;
from the package drcData.

## Examples

``` r
dmod <- drmTMB(cbind(no, total - no) ~ dose, data = daphnids, model = "loglogistic",
               family = "binomial", fform = b3 + b4 ~ 0 + time,
               start = c(b1 = 0, b2 = 1, b4 = 1000, b5 = 0),
               fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
ed(dmod, 50, newdata = data.frame(time = c("24h", "48h")))
#> Effective dose estimates:
#> 
#>             Estimate Std.Error
#> 50:time=24h     5134      1057
#> 50:time=48h     1509       188
```
