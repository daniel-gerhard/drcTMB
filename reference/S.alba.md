# Potency of two herbicides in white mustard

The dry matter of white mustard (*Sinapis alba*) for doses of the
herbicides glyphosate and bentazone. The asymptotes of the two
herbicides can be assumed equal, while their slopes and ED50s differ.

## Usage

``` r
S.alba
```

## Format

A data frame with 68 rows and 3 variables:

- Dose:

  dose (g/ha)

- Herbicide:

  a factor with levels `Bentazone` and `Glyphosate`

- DryMatter:

  dry matter (g/pot)

## Source

Christensen, M. G., Teicher, H. B. and Streibig, J. C. (2003) Linking
fluorescence induction curve and biomass in herbicide screening. *Pest
Management Science*, **59**, 1303–1310; from the package drcData.

## Examples

``` r
smod <- drmTMB(DryMatter ~ Dose, data = S.alba, model = "loglogistic",
               fform = b3 + b4 ~ Herbicide, parlink = c(b4 = "log"),
               start = c(b5 = 0), fix = c(b5 = TRUE))
ed(smod, 50, newdata = data.frame(Herbicide = c("Bentazone", "Glyphosate")))
#> Effective dose estimates:
#> 
#>                         Estimate Std.Error
#> 50:Herbicide=Bentazone      28.6      2.04
#> 50:Herbicide=Glyphosate     66.9      5.98
```
