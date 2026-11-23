# Hormesis in lettuce plants

The shoot biomass of lettuce (*Lactuca sativa*) plants after 21 days in
a nutrient solution with isobutyl alcohol. Low concentrations stimulate
the growth (hormesis).

## Usage

``` r
lettuce
```

## Format

A data frame with 14 rows and 2 variables:

- conc:

  concentration of isobutyl alcohol (mg/l)

- weight:

  shoot biomass (g)

## Source

van Ewijk, P. H. and Hoekstra, J. A. (1993) Calculation of the EC50 and
its confidence interval when subtoxic stimulus is present.
*Ecotoxicology and Environmental Safety*, **25**, 25–32; from the
package drcData.

## References

van Ewijk, P. H. and Hoekstra, J. A. (1994) Curvature measures and
confidence intervals for the linear logistic model. *Applied
Statistics*, **43**, 477–487.

## Examples

``` r
hmod <- drmTMB(weight ~ conc, data = lettuce, model = "braincousens")
peak(hmod)
#> Peak of the hormesis curves, dose and response:
#> 
#>              Estimate Std.Error
#> dose:lf1         0.72    0.2670
#> response:lf1     1.20    0.0694
```
