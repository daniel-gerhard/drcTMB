# Toxicity of four types of selenium

The number of dead flies for concentrations of four types of selenium.

## Usage

``` r
selenium
```

## Format

A data frame with 25 rows and 4 variables:

- type:

  type of selenium (1 to 4)

- conc:

  total selenium concentration

- total:

  number of flies

- dead:

  number of dead flies

## Source

Jeske, D. R., Xu, H. K., Blessinger, T., Jensen, P. and Trumble, J.
(2009) Testing for the equality of EC50 values in the presence of
unequal slopes with application to toxicity of selenium types. *Journal
of Agricultural, Biological, and Environmental Statistics*, **14**,
469–483; from the package drcData.

## Examples

``` r
semod <- drmTMB(cbind(dead, total - dead) ~ conc, data = selenium, model = "loglogistic",
                family = "binomial", fform = b3 + b4 ~ factor(type),
                start = c(b2 = 1, b5 = 0), fix = c(b2 = TRUE, b5 = TRUE))
ed(semod, 50, newdata = data.frame(type = 1:4))
#> Effective dose estimates:
#> 
#>           Estimate Std.Error
#> 50:type=1      266     14.61
#> 50:type=2      408     45.20
#> 50:type=3      133      7.20
#> 50:type=4       84      5.74
```
