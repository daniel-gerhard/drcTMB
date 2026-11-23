# Effect of ferulic acid on the root growth of ryegrass

The root length of perennial ryegrass (*Lolium perenne* L.) for several
concentrations of ferulic acid; part of a study of the joint action of
phenolic acids on root growth inhibition.

## Usage

``` r
ryegrass
```

## Format

A data frame with 24 rows and 2 variables:

- rootl:

  root length (cm)

- conc:

  concentration of ferulic acid (mM)

## Source

Inderjit, Streibig, J. C. and Olofsdotter, M. (2002) Joint action of
phenolic acid mixtures and its significance in allelopathy research.
*Physiologia Plantarum*, **114**, 422–428; from the package drcData.

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
ed(mod, c(10, 50, 90))
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 10:lf1     1.34     0.189
#> 50:lf1     3.26     0.196
#> 90:lf1     7.94     1.055
```
