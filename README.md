drcTMB
================

The package drcTMB provides functions to fit dose-response curves using
the Template Model Builder (TMB) framework.

## Installation

``` r
# You can install drcTMB from GitHub by
pak::pak("daniel-gerhard/drcTMB")
```

## The drmTMB function

The main function in the drcTMB package is `drmTMB()`. There is some
similarity to the `drm()` function from the drc package, with a
simplified input of response and dose variables through a formula
interface, and specifying a nonlinear model through the `model`
argument. Available models include “logistic”, “loglogistic”,
“weibull1”, “weibull2”, and “lognormal”.

Each of these models has the same parameter names, b1 for the lower and
b2 for the upper asymptote, b3 for the steepness, b4 for the location of
the inflection point on the dose scale, and b5 for an optional asymmetry
parameter. The steepness b3 is positive for an increasing and negative
for a decreasing curve. The asymmetry parameter b5 is estimated on a
logarithmic scale. The models michaelismenten, asymptoticregression, and
exponentialdecay have no steepness parameter; here b1 is the response at
dose zero and b2 the asymptote at high doses. The `start` argument
provides starting values for the parameters. The `fix` argument
(TRUE/FALSE) allows to fix parameters at these starting values.
Parameters that a model does not use, such as b5 for the Weibull models,
are fixed automatically.

``` r
library(drcTMB)
library(drcData)
data(ryegrass)
# Fitting a 3-parameter log-logistic model
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic", 
              start=c(0, 8, -1, 2, 0), fix=c(TRUE, FALSE, FALSE, FALSE, TRUE))
print(mod)
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, start = c(0, 8, 
#>     -1, 2, 0), fix = c(TRUE, FALSE, FALSE, FALSE, TRUE), model = "loglogistic")
#> 
#> Model:  loglogistic
#> Family: gaussian, identity link for the asymptotes
#> 
#> Coefficients:
#> b2:(Intercept)  b3:(Intercept)  b4:(Intercept)  
#>          7.855          -2.470           3.263  
#> Fixed: b1:(Intercept) = 0, b5:(Intercept) = 0
#> 
#> Residual standard deviation: 0.5253
```

With the `predict()` method, we can visualize the fitted curve, together
with 95% confidence intervals for the mean response
(`interval = "confidence"`) and prediction intervals for a new
observation (`interval = "prediction"`).

``` r
newdata <- data.frame(conc = exp(seq(-0.5, log(30), length.out = 100)))
ci <- predict(mod, newdata = newdata, interval = "confidence")
pi <- predict(mod, newdata = newdata, interval = "prediction")
head(cbind(newdata, ci))
#>        conc      fit    se.fit      lwr      upr
#> 1 0.6065307 7.734345 0.1675733 7.405907 8.062783
#> 2 0.6309088 7.722175 0.1665990 7.395647 8.048703
#> 3 0.6562668 7.708806 0.1657090 7.384022 8.033589
#> 4 0.6826439 7.694123 0.1649320 7.370862 8.017384
#> 5 0.7100813 7.678003 0.1643010 7.355979 8.000027
#> 6 0.7386214 7.660314 0.1638521 7.339169 7.981458
newdata[c("pred", "clwr", "cupr")] <- ci[c("fit", "lwr", "upr")]
newdata[c("plwr", "pupr")] <- pi[c("lwr", "upr")]
library(ggplot2)
ggplot(ryegrass, aes(x = log(conc), y = rootl)) +
  geom_ribbon(data = newdata, aes(x = log(conc), y = pred, ymin = plwr, ymax = pupr), fill = "blue", alpha = 0.1) +
  geom_ribbon(data = newdata, aes(x = log(conc), y = pred, ymin = clwr, ymax = cupr), fill = "blue", alpha = 0.25) +
  geom_point() +
  geom_line(data = newdata, aes(x = log(conc), y = pred), color = "blue") +
  labs(x = "Log Concentration", y = "Root Length") +
  theme_bw()
```

![](README-unnamed-chunk-4-1.png)<!-- -->

Bootstrap intervals are available by passing a parametric bootstrap
sample of the model from `bootstrap()` through the `boot` argument of
`predict()`.

## Effective dose estimation

The `ed()` function can be used to estimate effective doses (e.g., ED50,
ED90) from the fitted model, given a specified response level
(`respLev`). The effective dose ED p is the dose at which the response
has changed by p% of the difference between the asymptote at low doses
(the control) and the asymptote at high doses.

``` r
(effdose <- ed(mod, respLev = c(25, 50, 75)))
#> Effective dose estimates:
#> 
#>        Estimate Std.Error
#> 25:lf1     2.09     0.175
#> 50:lf2     3.26     0.184
#> 75:lf3     5.09     0.394
```

Simultaneous confidence intervals can be constructed through the package
multcomp.

``` r
library(multcomp)
linfct <- diag(3)
rownames(linfct) <- c("ED25", "ED50", "ED75")
effdose_glht <- glht(effdose, linfct = linfct)
confint(effdose_glht)
#> 
#>   Simultaneous Confidence Intervals
#> 
#> Fit: NULL
#> 
#> Quantile = 2.3181
#> 95% family-wise confidence level
#>  
#> 
#> Linear Hypotheses:
#>           Estimate lwr    upr   
#> ED25 == 0 2.0918   1.6870 2.4966
#> ED50 == 0 3.2634   2.8375 3.6892
#> ED75 == 0 5.0910   4.1776 6.0045
```
