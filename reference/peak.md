# Peak of hormesis and substrate inhibition curves

The dose at which a curve of the hormesis models braincousens and
cedergreen, or of the substrate inhibition model haldane, reaches its
maximum, and the maximal response, as `MAX` in drc. For the hormesis
models the peak has no closed form: it is found numerically, and its
derivatives, for the standard errors, follow from the implicit function
theorem, as for the effective doses of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md). Curves
without a stimulation (b5 \<= 0) have no peak, and give NA with a
warning; the hormesis models need a decreasing curve (b3 \< 0),
otherwise an error is raised. For haldane, the peak is at the dose
sqrt(b4/b5), the optimal substrate concentration sqrt(Km Ki); without
inhibition (b5 = 0) there is no peak. With a `curvelink`, the peak dose
is that of the curve, as the link is monotone, and the response at the
peak is on the scale of the response, the inverse link of the curve at
the peak.

## Usage

``` r
peak(x, linfct = NULL, newdata = NULL, random = TRUE, vcov. = NULL)
```

## Arguments

- x:

  a drmTMB object of the braincousens, cedergreen or haldane model

- linfct, newdata, random:

  the curves, as in
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)

- vcov.:

  the covariance of the estimated parameters, as in
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md)

## Value

an edTMB object, with entries "dose:curve" and "response:curve"; not
available for fits in the BMD parameterisation

## Details

The result is an edTMB object, as from
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md), with
the dose and the response at the peak of each curve, named "dose" and
"response", and their joint covariance matrix, e.g. for `glht` of
multcomp.
[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md)
gives Wald and bootstrap intervals for both, and profile likelihood
intervals for the dose; for the response, profile intervals are not
available and Wald intervals are returned, with a message.

## See also

[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`confint.edTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md),
[`edTMB-methods`](https://daniel-gerhard.github.io/drcTMB/reference/edTMB-methods.md)

## Examples

``` r
# hormesis: stimulation of the growth of lettuce at low concentrations of isobutyl alcohol
hmod <- drmTMB(weight ~ conc, data = lettuce, model = "braincousens")
(p <- peak(hmod))
#> Peak of the hormesis curves, dose and response:
#> 
#>              Estimate Std.Error
#> dose:lf1         0.72    0.2670
#> response:lf1     1.20    0.0694
confint(p, method = "wald")
#>                  2.5 %   97.5 %
#> dose:lf1     0.1162045 1.324228
#> response:lf1 1.0438171 1.357752
# the dose at which the response returns to the level of the control
ed(hmod, respLev = 0)
#> Effective dose estimates:
#> 
#>       Estimate Std.Error
#> 0:lf1     2.93      1.05
```
