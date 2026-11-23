# Germination of three crops

Germination of mungbean, rice and wheat seeds at temperatures between 10
and 40 degrees Celsius, recorded daily for up to 18 days: time-to-event
data with several curves.

## Usage

``` r
germination
```

## Format

A data frame with 192 rows and 5 variables:

- temp:

  temperature (degrees Celsius)

- species:

  a factor with levels `mungbean`, `rice` and `wheat`

- start:

  left endpoint of the monitoring interval (days)

- end:

  right endpoint of the monitoring interval (days)

- germinated:

  number of seeds germinated in the interval

## Source

From the package drcData.

## Details

For each species and temperature, 20 seeds were placed in a Petri dish
in a dark climate cabinet (10, 16, 22, 28, 34 or 40 degrees Celsius; not
all temperatures for all species), and the germinated seeds were counted
and removed daily.

## References

Ritz, C., Pipper, C. B. and Streibig, J. C. (2013) Analysis of
germination data from agricultural experiments. *European Journal of
Agronomy*, **45**, 1–6.

## Examples

``` r
wheat <- subset(germination, species == "wheat" & temp == 22)
wmod <- drmTMB(germinated ~ start + end, data = wheat, model = "loglogistic", family = "event",
               start = c(b5 = 0), fix = c(b5 = TRUE))
ed(wmod, 50)
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 50:lf1     3.83     0.216
```
