# Liver tumour incidence after exposure to TCDD (Kociba et al. 1978)

Liver tumour incidence in Sprague-Dawley rats exposed to
2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD); the historical study of the
EFSA example on informative priors, see
[`liver.tumor`](https://daniel-gerhard.github.io/drcTMB/reference/liver.tumor.md)
for the current study, and
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)
and
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md)
for the analysis.

## Usage

``` r
TCDD
```

## Format

A data frame with 4 rows and 3 variables:

- conc:

  dose metric of TCDD as used in the EFSA example (EN-9902), not the
  administered dose; the doses of the study were 0, 0.001, 0.01 and 0.1
  micrograms/kg body weight per day, in the diet

- total:

  number of rats

- incidence:

  number of rats with liver tumours

## Source

Kociba, R., Keyes, D., Beyer, J., Carreon, R., Wade, C., Dittenber, D.,
Kalnins, R., Frauson, L., Park, C., Barnard, S., Hummel, R. and
Humiston, C. (1978) Results of a two-year chronic toxicity and
oncogenicity study of 2,3,7,8-tetrachlorodibenzo-p-dioxin in rats.
*Toxicology and Applied Pharmacology*, **46**, 279–303; from the package
drcData.

## Examples

``` r
tmod <- drmTMB(cbind(incidence, total - incidence) ~ conc, data = TCDD, model = "loglogistic",
               family = "binomial", start = c(b2 = 1, b5 = 0), fix = c(b2 = TRUE, b5 = TRUE))
bmd(tmod, bmr = 0.1, def = "excess")
#> Benchmark doses, excess risk:
#> 
#>         Estimate Std.Error
#> 0.1:lf1     6.13      2.99
```
