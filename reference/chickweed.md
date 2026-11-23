# Germination of common chickweed

Germination of 200 chickweed (*Stellaria media*) seeds in one
germination test, from a study of chlorsulfuron resistant and sensitive
biotypes, recorded in time intervals: an example of time-to-event data
(`family = "event"` of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md)).

## Usage

``` r
chickweed
```

## Format

A data frame with 35 rows and 3 variables:

- start:

  left endpoint of the monitoring interval (hours)

- end:

  right endpoint of the monitoring interval (hours)

- count:

  number of seeds germinated in the interval

## Source

Lisa Borggaard (formerly Faculty of Life Sciences, University of
Copenhagen); from the package drcData.

## Details

200 seeds were placed in Petri dishes in a dark growth cabinet at 5
degrees Celsius, and the germinated seeds were counted and removed at 34
inspection times. Chickweed has dormant seeds, so that the proportion of
germinated seeds stays well below 1. The last row, with `end = Inf`,
holds the seeds that had not germinated by the end of the experiment.

## References

Ritz, C., Pipper, C. B. and Streibig, J. C. (2013) Analysis of
germination data from agricultural experiments. *European Journal of
Agronomy*, **45**, 1–6.

## Examples

``` r
gmod <- drmTMB(count ~ start + end, data = chickweed, model = "loglogistic", family = "event",
               start = c(b5 = 0), fix = c(b5 = TRUE))
# times by which 10%, 50% and 90% of the germinating seeds have germinated
ed(gmod, c(10, 50, 90))
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 10:lf1      176      3.49
#> 50:lf1      196      2.51
#> 90:lf1      218      4.27
# time by which 10% of all seeds have germinated
ed(gmod, 0.1, type = "absolute")
#> Effective dose estimates:
#> 
#>           Estimate Std.Error
#> y=0.1:lf1      196      3.66
```
