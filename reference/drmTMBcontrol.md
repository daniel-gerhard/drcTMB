# Control parameters of drmTMB

Settings of the optimisation of
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md),
passed as its argument control, e.g.
`drmTMB(..., control = drmTMBcontrol(restart = TRUE))` or
`control = list(restart = TRUE)`.

## Usage

``` r
drmTMBcontrol(
  optimiser = c("nlminb", "optim"),
  method = c("L-BFGS-B", "BFGS", "CG", "Nelder-Mead"),
  iter.max = 1000,
  eval.max = 2000,
  optcontrol = list(),
  inner = list(),
  restart = FALSE,
  smallsample = TRUE
)
```

## Arguments

- optimiser:

  the optimiser: "nlminb" (default,
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html)) or "optim"
  ([`optim`](https://rdrr.io/r/stats/optim.html))

- method:

  the method of optim: "L-BFGS-B" (default), "BFGS", "CG" or
  "Nelder-Mead"

- iter.max:

  maximum number of iterations of the optimiser (iter.max of nlminb,
  maxit of optim)

- eval.max:

  maximum number of evaluations of the objective function of nlminb

- optcontrol:

  a list of further control parameters of the optimiser, the argument
  control of [`nlminb`](https://rdrr.io/r/stats/nlminb.html) (e.g.
  rel.tol, x.tol, trace) or
  [`optim`](https://rdrr.io/r/stats/optim.html) (e.g. reltol, factr,
  pgtol, trace, parscale); its iter.max, eval.max or maxit replace those
  above. The tolerances of optim are factr=1e5 (L-BFGS-B) and
  reltol=1e-10 (the other methods) unless given, close to rel.tol=1e-10
  of nlminb.

- inner:

  a list of control parameters of the inner optimisation of the random
  effects for the Laplace approximation, see the argument inner.control
  of [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html), e.g.
  list(maxit=2000, tol10=1e-6)

- restart:

  if TRUE, the fit is restarted once at the estimates and the better fit
  is kept, see Details

- smallsample:

  if TRUE (default), gaussian and lognormal models without random
  effects, priors, bmdpar, correlation structure and variance function
  use the unbiased residual variance with n - p residual degrees of
  freedom and t quantiles, as lm, nls and drc; if FALSE, the maximum
  likelihood estimate and normal quantiles. See
  [`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md).

## Value

a list with control parameters

## Details

The optimiser minimises the negative log-likelihood (with random
effects, its Laplace approximation) over the fixed parameters, with the
gradient from TMB, and for nlminb without random effects also the
Hessian. The same settings are used for the refits of
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md),
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
(test="boot"),
[`bmdma`](https://daniel-gerhard.github.io/drcTMB/reference/bmdma.md)
and
[`vcovrobust`](https://daniel-gerhard.github.io/drcTMB/reference/vcovrobust.md).

optim does not use the Hessian, and its methods other than "L-BFGS-B" do
not allow bounds; models with bounds (the default lower bounds 0 of the
multistage, nec, haldane, exponential and power models, see
[`drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB.md),
user bounds lower and upper, or priors with a bounded support) need
method "L-BFGS-B"; the other methods give an error. Parameter values at
which the negative log-likelihood is not finite (e.g. a negative b4 of
the BMD parameterisation) are given a large finite value for optim.

With restart=TRUE, the model is fitted again starting at the estimates,
which restarts the approximations of the optimiser (e.g. after a false
convergence, a stop at a kink of the likelihood, or in a flat likelihood
of a variance), and the better fit is kept: the restart if its negative
log-likelihood is lower, or as low and it converged. This doubles the
time of the fits at least, also of the bootstrap fits.

## Examples

``` r
drmTMBcontrol(optimiser="optim", optcontrol=list(factr=1e5), restart=TRUE)
#> $optimiser
#> [1] "optim"
#> 
#> $method
#> [1] "L-BFGS-B"
#> 
#> $iter.max
#> [1] 1000
#> 
#> $eval.max
#> [1] 2000
#> 
#> $optcontrol
#> $optcontrol$factr
#> [1] 1e+05
#> 
#> 
#> $inner
#> list()
#> 
#> $restart
#> [1] TRUE
#> 
#> $smallsample
#> [1] TRUE
#> 
fit <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              control = drmTMBcontrol(optimiser = "optim", restart = TRUE))
summary(fit)
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, model = "loglogistic", 
#>     control = drmTMBcontrol(optimiser = "optim", restart = TRUE))
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 24
#> Log-likelihood: -15.88 (df = 6), AIC: 43.76, BIC: 50.82
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b1:(Intercept)   0.3215      0.380
#> b2:(Intercept)   7.7605      0.193
#> b3:(Intercept)  -3.9301      2.194
#> b4:(Intercept)   2.2147      0.920
#> b5:(Intercept)  -0.7594      1.134
#> The asymmetry parameter b5 is on the log scale.
#> 
#> Residual standard deviation: 0.527 on 19 degrees of freedom
#> 
#> Optimiser: optim (L-BFGS-B): CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH, maximum absolute gradient 1.8e-05
#> Restarted at the estimates: no better optimum
```
