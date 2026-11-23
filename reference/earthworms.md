# Earthworm toxicity test

The number of earthworms remaining in a container contaminated with a
toxic substance, i.e. not migrating to a neighbouring uncontaminated
container. At dose 0, about half of the earthworms are expected in each
container, so the upper asymptote is about 0.5 rather than 1.

## Usage

``` r
earthworms
```

## Format

A data frame with 35 rows and 3 variables:

- dose:

  dose of the toxic substance

- number:

  number of earthworms remaining in the contaminated container

- total:

  total number of earthworms in the containers

## Source

Nina Cedergreen, Faculty of Life Sciences, University of Copenhagen;
from the package drcData.

## Examples

``` r
emod <- drmTMB(cbind(number, total - number) ~ dose, data = earthworms, model = "loglogistic",
               family = "binomial", start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
emod
#> Call:
#> drmTMB(form = cbind(number, total - number) ~ dose, data = earthworms, 
#>     start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE), 
#>     family = "binomial", model = "loglogistic")
#> 
#> Model:  loglogistic
#> Family: binomial, identity link for the asymptotes
#> 
#> Coefficients:
#> b2:(Intercept)  b3:(Intercept)  b4:(Intercept)  
#>         0.6049         -1.5056          0.2924  
#> Fixed: b1:(Intercept) = 0, b5:(Intercept) = 0
```
