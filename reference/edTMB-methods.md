# Methods for effective doses, peaks and benchmark doses

Print, coef and vcov methods for objects of class edTMB, returned by
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md) and
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md);
confidence intervals are given by
[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md).

## Usage

``` r
# S3 method for class 'edTMB'
print(x, ...)

# S3 method for class 'edTMB'
coef(object, ...)

# S3 method for class 'edTMB'
vcov(object, ...)
```

## Arguments

- x, object:

  an edTMB object

- ...:

  further arguments, not used

## Value

print returns x invisibly; coef a named vector of the estimates; vcov
their covariance matrix (delta method).

## Examples

``` r
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
e <- ed(mod, respLev = c(10, 50))
e
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 10:lf1     1.34     0.189
#> 50:lf1     3.26     0.196
coef(e)
#>   10:lf1   50:lf1 
#> 1.340856 3.263356 
vcov(e)
#>            10:lf1     50:lf1
#> 10:lf1 0.03573499 0.01829770
#> 50:lf1 0.01829770 0.03857539
```
