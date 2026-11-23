# Liver tumour incidence after exposure to TCDD (NTP 2006)

Liver tumour incidence in female Harlan Sprague-Dawley rats exposed to
2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD); the current study of the
EFSA example on informative priors, see
[`TCDD`](https://daniel-gerhard.github.io/drcTMB/reference/TCDD.md) for
the historical study, and
[`bmdefsa`](https://daniel-gerhard.github.io/drcTMB/reference/bmdefsa.md)
and
[`fitpert`](https://daniel-gerhard.github.io/drcTMB/reference/fitpert.md)
for the analysis with an informative prior from the historical study.

## Usage

``` r
liver.tumor
```

## Format

A data frame with 6 rows and 3 variables:

- conc:

  dose metric of TCDD as used in the EFSA example (EN-9902), not the
  administered dose; the doses of the study were 0, 3, 10, 22, 46 and
  100 ng/kg body weight per day, by gavage

- total:

  number of rats

- incidence:

  number of rats with liver tumours

## Source

National Toxicology Program (2006) NTP technical report on the
toxicology and carcinogenesis studies of
2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD) (CAS No. 1746-01-6) in female
Harlan Sprague-Dawley rats (gavage studies). *National Toxicology
Program Technical Report Series*, **521**, 4–232; from the package
drcData.

## Examples

``` r
# no tumours in the controls: background risk b1 fixed at 0
lmod <- drmTMB(cbind(incidence, total - incidence) ~ conc, data = liver.tumor,
               model = "loglogistic", family = "binomial", start = c(b1 = 0, b2 = 1, b5 = 0),
               fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
bmd(lmod, bmr = 0.1, def = "excess")
#> Benchmark doses, excess risk:
#> 
#>         Estimate Std.Error
#> 0.1:lf1     23.9      1.98
```
