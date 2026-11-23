# Fit dose-response models with TMB

Fits a nonlinear dose-response curve f(x) with the parameters b1-b5 by
maximum likelihood with the Template Model Builder (TMB), optionally
with covariates on the parameters (`fform`), random effects (`rform`), a
non-normal response distribution (`family`), correlated or
heteroscedastic residuals (`correlation`, `variance`) and priors
(`prior`). The interface resembles `drm` of the drc package;
[`vignette("drcTMB", package = "drcTMB")`](https://daniel-gerhard.github.io/drcTMB/articles/drcTMB.md)
introduces the package with examples.

## Usage

``` r
drmTMB(
  form,
  fform = NULL,
  rform = NULL,
  data,
  start = NULL,
  fix = NULL,
  weights = NULL,
  lower = -Inf,
  upper = Inf,
  family = "gaussian",
  model = "logistic",
  alpha = NULL,
  smooth = NULL,
  parlink = NULL,
  curvelink = "identity",
  prior = NULL,
  bmdpar = NULL,
  correlation = NULL,
  variance = NULL,
  series = NULL,
  gsl = NULL,
  gsu = NULL,
  gsn = NULL,
  control = list()
)
```

## Arguments

- form:

  a formula response ~ dose. For the binomial and beta-binomial
  families, the response is cbind(successes, failures), or a 0/1
  response for binary data (binomial only); proportions with weights, as
  in drc, are not used. For the event family, count ~ start + end, the
  number of events in each time interval (see the section Time-to-event
  data).

- fform:

  fixed effects for the model parameters: a formula of the form b2 + b4
  ~ Stress, or a list of such formulas, e.g. list(b2 ~ Stress, b4 ~
  Stress + Time). Parameters not listed have a single intercept (~ 1). A
  list of five right-hand side formulas, one for each of b1-b5 (X1 = ~
  1, ..., X5 = ~ 1), is also accepted.

- rform:

  random effects: a formula of the form b1 + b2 + b4 ~ 1 \| cluster for
  random intercepts, or b4 ~ x \| cluster for random intercepts and
  slopes in a covariate x, or a list of such formulas, e.g. list(b1 + b4
  ~ 1 \| exper, b1 ~ 1 \| plate). The random effects within a formula
  are correlated across the listed parameters and their intercepts and
  slopes; different formulas, the levels of nested factors (1 \|
  exper/plate), and intercepts and slopes with \|\| (b4 ~ x \|\|
  cluster) are independent terms. The random slopes of a parameter
  usually go with a fixed effect of the covariate in fform.

- data:

  a data.frame object. Observations with missing values in the variables
  of the model or in the weights are left out.

- start:

  starting values: a vector with one value for each model coefficient,
  in the order of
  [`drmTMBstart`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBstart.md),
  where NA entries are replaced by the starting values of drmTMBstart;
  or a vector named by parameter, e.g. c(b1=0, b5=0), which gives the
  same starting value to all coefficients of each named parameter.

- fix:

  a logical vector with one entry for each model coefficient; TRUE fixes
  a coefficient at its value in start. Alternatively, a vector named by
  parameter, e.g. c(b1=TRUE, b5=TRUE), fixes all coefficients of the
  named parameters, whose values must then be given in start, e.g.
  start=c(b1=0, b5=0). Parameters that the model does not use are always
  fixed, at 0 unless given in start (see Details).

- weights:

  a vector of prior weights, one for each row of data. For the gaussian
  and lognormal families, they are precision weights: the residual
  standard deviation of an observation is sigma/sqrt(w), as in lm and
  drc. For the other families, the log-likelihood contributions are
  multiplied by the weights. Weights cannot be combined with a
  correlation structure.

- lower:

  lower bounds for the model coefficients b1-b5: a single value or a
  vector with one entry for each model coefficient, in the order of
  [`drmTMBstart`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBstart.md).
  Bounds of fixed coefficients are ignored; variance parameters are
  unbounded. By default, some parameters without covariates in fform are
  bounded below by 0: b3-b5 of multistage, the rate b3 and (with the
  identity link) the threshold b4 of nec, the inhibition b5 of haldane,
  and the power b3 of exponential and power; see Details.

- upper:

  upper bounds for the model coefficients b1-b5, as lower.

- family:

  the distribution of the response: "gaussian" (default), "lognormal",
  "binomial", "betabinomial", "beta", "poisson", "nbinom2", "nbinom1",
  or "event" for time-to-event data, e.g. germination times (see the
  section Time-to-event data). For the lognormal family, log(y) is
  normal with mean log f(x) and standard deviation sigma: the curve f is
  the median of a positive response, with a constant coefficient of
  variation, as in PROAST and BMDS; the mean is f(x) exp(sigma^2/2). The
  negative binomial families are for overdispersed counts, with variance
  mu (1 + mu/theta) for nbinom2 and mu (1 + phi) for nbinom1, and the
  beta-binomial family for overdispersed binomial counts, with variance
  n mu (1 - mu) (1 + (n - 1)/(phi + 1)) for n trials, as in glmmTMB.

- model:

  the model function: "logistic" (default), "loglogistic", "weibull1",
  "weibull2", "lognormal", "gompertz", "michaelismenten",
  "asymptoticregression", "exponentialdecay", "gamma", "multistage",
  "nec" (threshold), "braincousens" and "cedergreen" (hormesis),
  "haldane" (substrate inhibition), "exponential" and "power" (no upper
  asymptote), or "constant" (no dose effect); see Details.

- alpha:

  the exponent of the hormesis term of the cedergreen model, a positive
  number, e.g. 1, 0.5 or 0.25 as the models CRS.4a, CRS.4b and CRS.4c of
  drc (CRS.5a-c with b1 estimated); only for the cedergreen model. See
  Details.

- smooth:

  smoothing of the threshold of the nec model: the width of the smoothed
  threshold relative to the smallest interval between doses, a number of
  at least 0; 0 is an exact threshold. The default is 0.1 with random
  effects on b4, and 0 otherwise. Only for the nec model, see Details.

- parlink:

  links of the parameters b1-b5: NULL for the defaults of the model, a
  single link for the asymptotes b1 and b2, or links named by parameter,
  e.g. c(b1="logit", b2="logit", b4="log"), where the parameters not
  named keep their defaults. The available links are:

  - b1, b2: "identity" (default), "log" or "logit" (for the event
    family, b2 has the logit link by default);

  - b4: "identity" (default) or "log", for the models with a positive b4
    (not logistic, gompertz, multistage, exponential, power and
    constant);

  - b3 and b5: identity; the asymmetry b5 of logistic and loglogistic
    and the shape b5 of nec always have the log link.

  The coefficients, and start, fix, lower, upper and linfct of
  [`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md), are
  on the scale of the links, see the section Parameter links and curve
  links.

- curvelink:

  link of the whole curve: "identity" (default), the curve is the mean
  (the median for the lognormal family, the probability for binomial,
  beta-binomial and beta); "log" (gaussian, lognormal, poisson, nbinom2,
  nbinom1) or "logit" and "probit" (binomial, betabinomial, beta), the
  mean (median, probability) is the inverse link of the curve. Unlike
  parlink, which applies to the linear predictors of the parameters,
  curvelink applies to the response; b1 and b2 are then on the scale of
  the curve link and keep the identity parlink. See the section
  Parameter links and curve links.

- prior:

  NULL (default) for maximum likelihood, or a named list of prior
  distributions (see
  [`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md)),
  e.g. list(b2=normal(8, 2), b4=pert(0.5, 3, 10), sd=half_normal(1)).
  The names are:

  - b1-b5, for parameters with a single coefficient, on the scale of the
    parameter (after the inverse link);

  - coefficient names, e.g. "b4:herbicideB", on the scale of the
    coefficient (the link);

  - "sd" for all standard deviations of the random effects, or e.g.
    "sd:b1\|group" for one of them;

  - "sigma" (gaussian, lognormal), "phi" (beta, nbinom1, betabinomial)
    or "theta" (nbinom2) for the residual parameter;

  - "bmd" for the benchmark dose, with bmdpar;

  - the parameters of a correlation structure and variance function:
    "cor:rho" (AR1, compound symmetry), "cor:phi" (CAR1), "cor:range",
    "cor:nugget", "var:delta" (power, exponential) and
    "var:ratio:STRATUM" (e.g. "var:ratio:B", the ratio of the standard
    deviation of stratum B to the first stratum), see
    [`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md).

  See the section Priors.

- bmdpar:

  NULL (default) for the standard parameterisation, or a list with a
  benchmark response bmr and its definition def as in
  [`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md),
  e.g. list(bmr=0.1, def="relative"), for the BMD parameterisation: b4
  is replaced by the log benchmark dose, coefficient log_bmd, and
  computed from it. See the section BMD parameterisation.

- correlation:

  NULL (default) for independent residuals, or a residual correlation
  structure within series, e.g. `cor_car1(~ time | subject)`, for the
  gaussian and lognormal families; see
  [`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md).

- variance:

  NULL (default) for a constant residual standard deviation, or a
  variance function, e.g.
  [`var_power()`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md)
  (a power of the fitted mean) or `var_ident(~ 1 | group)`, for the
  gaussian and lognormal families; see
  [`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md).

- series:

  for the event family: a one-sided formula of the series of the rows,
  e.g. ~ dish for the intervals of each Petri dish. By default, rows
  with the same covariates (fform) and the same levels of the random
  effects form one series. The series are needed by simulate, bootstrap,
  residuals and the dispersion of summary, which need complete series:
  intervals from 0 to Inf without gaps.

- gsl, gsu:

  lower and upper bounds of a grid search for starting values: vectors
  with one entry for each model coefficient, on the scale of the
  coefficients, NA for the default range. If gsl, gsu and gsn are all
  NULL (default), no grid search is performed. See
  [`drmTMBstart`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBstart.md).

- gsn:

  number of grid points for each coefficient, a single number or a
  vector with one entry for each model coefficient; 10 if gsl or gsu is
  given. `gsn=10` alone runs a grid search over the default ranges.

- control:

  a list of control parameters, the result of
  [`drmTMBcontrol`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBcontrol.md)
  or a list of its arguments, e.g. list(restart=TRUE).

## Value

An object of class drmTMB, a list with, among others, the call, the
data, the model and family, the TMB object (obj), the optimiser result
(opt), start, fix, lower and upper. Use the methods (`coef`, `vcov`,
`summary`, `predict`, ..., see
[`drmTMB-methods`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md)),
[`confint.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/confint.drmTMB.md),
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md),
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md),
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
and
[`plot.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/plot.drmTMB.md)
rather than the components.

## Details

All models, except exponential and power, are of the form f(x) = b1 +
(b2 - b1) F(x), with a function F increasing from 0 to 1. For the models
with steepness b3 (logistic, loglogistic, weibull1, weibull2, lognormal,
gompertz), F depends on b3 (x - b4) or b3 (log(x) - log(b4)); b1 is the
lower and b2 the upper asymptote, and b3 is positive for an increasing
and negative for a decreasing curve. The asymmetry parameter b5 of the
logistic and loglogistic models is on the log scale; b5 = 0 is the
symmetric four-parameter model. Compared with the models of drc (e.g.
LL.4 with parameters b, c, d, e, and LL.5 with f), b1 = c, b2 = d, b3 =
-b (a positive b3 is an increasing curve, unlike drc), b4 = e and b5 =
log(f). The formulas of all models are given in the appendix of the
vignette.

For michaelismenten, asymptoticregression and exponentialdecay, b1 is
the response at dose zero and b2 the asymptote at high doses.
asymptoticregression and exponentialdecay are the same function, f(x) =
b1 + (b2 - b1) (1 - exp(-x/b4)), increasing for b2 \> b1 and decreasing
otherwise (AR.3 and EXD.3 of drc).

The parameters that a model does not use are fixed: b5 except for
logistic, loglogistic, braincousens, cedergreen, multistage, nec and
haldane; b3 for michaelismenten, asymptoticregression, exponentialdecay
and haldane; b2 and b5 for exponential and power; b2-b5 for constant.

The constant model, f(x) = b1, has the same mean at all doses. Compared
with a dose-response model by
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md),
it tests for a dose effect, as `noEffect` in drc. Under the constant
model, the parameters of the curve other than an asymptote are not
identifiable, so the chi-square distribution of the likelihood ratio
statistic is only an approximation.

## Models of BMDS and PROAST

The gamma model, F(x) = pgamma(x/b4, shape=b3), is the gamma model of
the benchmark dose software of the US EPA (BMDS) and `gammadr` of drc
(b3 = e, b4 = 1/b). As for michaelismenten, b1 is the response at dose
zero and b2 the asymptote at high doses; for quantal responses, b1 is
the background risk and b2 is fixed at 1, e.g.
`start=c(b2=1), fix=c(b2=TRUE)`. The shape b3 is positive; with b3 = 1,
the gamma model is the asymptoticregression model (the quantal-linear
model of BMDS for quantal responses), so that
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
tests the shape. BMDS restricts the shape to b3 \>= 1 by default, which
can be set by `lower`.

The multistage model of BMDS, F(x) = 1 - exp(-(b3 x + b4 x^2 + b5 x^3)),
has the coefficients beta1-beta3 of BMDS in b3-b5; as for the gamma
model, b1 is the response at dose zero (the background risk) and b2 the
asymptote at high doses, fixed at 1 for quantal responses. All three
coefficients are estimated by default; the model of degree 2 fixes b5 at
0, and of degree 1 also b4, e.g.
`start=c(b2=1, b5=0), fix=c(b2=TRUE, b5=TRUE)`. As in BMDS, b3-b5 are
bounded below by 0 by default (if `lower` is not given, and for
parameters without covariates in fform), so that the curve increases;
estimates at the bound are common, and are listed by `summary`, as their
standard errors and Wald intervals do not account for the bound. Profile
intervals of the effective doses respect the bounds. The log link for b4
is not available.

The exponential and power models of PROAST and BMDS have no upper
asymptote: f(x) = b1 exp(b4 x^b3) and f(x) = b1 + b4 x^b3, with b1 the
response at dose zero, b3 \> 0 the power and b4 the coefficient, whose
sign gives the direction of the curve; b2 and b5 are not used. b4 = 0 is
the constant model, and b3 = 1 gives the exponential model 2 of PROAST
and BMDS (exp2, e.g. `start=c(b3=1), fix=c(b3=TRUE)`), so that
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
tests the power; the exponential model with b3 estimated is exp3
(PROAST: a exp(b x^d); BMDS: a exp(sign \* (b x)^d) with b4 = sign \*
b^b3). The power b3 is bounded below by 0 by default (if `lower` is not
given, without covariates in fform); BMDS restricts the power of its
power model to b3 \>= 1 by default, which can be set by `lower`. Without
an upper asymptote, relative effective doses are not defined:
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md) gives
absolute effective doses (`type="absolute"`), and
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md)
benchmark doses. For the exponential model, the log link for b1,
`parlink=c(b1="log")`, makes the model linear on the log scale, with
log(b1) the intercept.

The linear-quadratic model of radiobiology, S = exp(-alpha x - beta
x^2), is the multistage model of degree 2 with b1 = 1 and b2 = 0: alpha
= b3 and beta = b4, e.g.
`start=c(b1=1, b2=0, b5=0), fix=c(b1=TRUE, b2=TRUE, b5=TRUE)`.

## Threshold model nec

The threshold model nec has no effect up to a threshold b4, the
no-effect concentration (NEC), and F(x) = 1 - exp(-(b3 (x -
b4))^exp(b5)) above it: b1 is the response up to the threshold and b2
the asymptote at high doses, b3 \> 0 the rate and b5 the shape on the
log scale; b3 and b4 are bounded below by 0 by default (if `lower` is
not given, for parameters without covariates in fform, and for b4 with
the identity link). With b5 = 0, it is the NEC model of drc, NEC.4 (b1 =
d, b2 = c, b3 = b, b4 = e); with b2 fixed at 0, NEC.3, e.g.
`start=c(b2=0, b5=0), fix=c(b2=TRUE, b5=TRUE)`, and with b1 also fixed
at 1, NEC.2. With b5 estimated, the curve is sigmoidal above the
threshold for b5 \> 0, as the necsigm model of bayesnec (with the rate
inside the power), and
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
tests the shape. Effective doses are b4 + (-log(1 - p))^exp(-b5) / b3,
above the threshold. The likelihood is smooth in b4 between the doses,
has a kink where b4 equals a dose, and often has several local maxima;
the starting values therefore come from thresholds within each interval
between the doses. A warning is given for a threshold estimated at a
dose, where the standard errors are not valid (profile or bootstrap
intervals remain available), or at or above the highest dose, where it
is not identified. With b5 \> 0, the curve starts flat above the
threshold, and b4 and b5 are often hard to separate.

With random effects on b4, the likelihood of each group with an exact
threshold has kinks at the doses, so that the Laplace approximation is
unreliable and its inner optimisation has several modes. The threshold
is therefore smoothed by default with random effects on b4: the dose
above the threshold, max(x - b4, 0), is replaced by s log(1 + exp((x -
b4)/s)), with a width s of `smooth` (default 0.1) times the smallest
interval between the doses. The curve then leaves b1 within a few widths
s around b4, and effective doses are b4 + h + s log(1 - exp(-h/s)), with
h the term above. In simulations with random thresholds, fits with an
exact threshold mostly failed, and smoothed fits converged with little
bias. `smooth=0` gives an exact threshold also with random effects, and
`smooth` \> 0 a smoothed threshold without them; the width is not
estimated. `summary` reports the width.

## Hormesis and substrate inhibition

The hormesis models braincousens (Brain and Cousens, 1989) and
cedergreen (Cedergreen, Ritz and Streibig, 2005) describe a stimulation
at low doses before a decreasing curve: f(x) = b1 + (b2 - b1 + b5 g(x))
/ (1 + exp(-b3 (log(x) - log(b4)))), with g(x) = x for braincousens and
g(x) = exp(-1/x^alpha) for cedergreen. b2 is the response at dose zero
(the control), b1 the lower asymptote at high doses, and b5 the size of
the hormesis effect on the scale of the response; b5 = 0 gives the
loglogistic model, so that
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
tests for hormesis. They are the models BC.5 and CRS.5a-c of drc, with
b3 = -b; with b1 fixed at 0, BC.4 and CRS.4a-c. The curve is decreasing
at high doses for b3 \< 0, and returns to b1 for the braincousens model
only for b3 \< -1. The models are meant for decreasing curves; a warning
is given otherwise. The exponent alpha of the cedergreen model is not
estimated, as it is hardly identifiable together with b5.

The haldane model of substrate inhibition, f(x) = b1 + (b2 - b1) x /
(b4 + x + b5 x^2), extends michaelismenten: b1 is the response at dose
zero, b2 the maximal response without inhibition (Vmax), b4 the
half-saturation constant (Km), and b5 = 1/Ki the inhibition, bounded
below by 0 by default (if `lower` is not given, without covariates in
fform); b3 is not used. With b5 = 0, it is the michaelismenten model, so
that
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
tests for inhibition (with the estimate on the boundary under the null
hypothesis, the chi-square test is conservative). The curve rises to its
peak at the dose sqrt(b4/b5) and then decreases, so that effective doses
are not defined;
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md)
gives the dose and the response at the peak.

## Parameter links and curve links

With the log link for b4, `parlink=c(b4="log")`, the linear predictor of
b4, including its random effects, is log(b4): the coefficients of b4 in
`coef`, `start`, `fix`, `lower`, `upper` and `linfct` of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md) are on
the log scale. Random effects on b4 then vary b4 by a factor between
groups, and the estimate of b4 is the median over the groups. For the
models on the log-dose scale, this is the same model as the logistic (or
gompertz) model of log(dose) for loglogistic (weibull1), with b4 the log
ED50, but it keeps observations at dose zero and reports effective doses
on the dose scale. Without random effects on b4, both links give the
same fit.

With `curvelink`, the curve f(x) of the model is on the scale of a link
of the response: the mean, the median for the lognormal family, or the
probability is exp(f(x)) for the log link, plogis(f(x)) for logit and
pnorm(f(x)) for probit. The lognormal family with `curvelink="log"` has
the curve on the scale of log(y), log(median) = f(x), the lognormal
models of the Bayesian benchmark dose modelling of EFSA (2022); without
it, as in PROAST and BMDS, the curve is the median. The binomial family
with the power model and `curvelink="probit"` or `"logit"` gives the
quantal probit and logistic models of EFSA (2022), pnorm(b1 + b4 x^b3)
and plogis(b1 + b4 x^b3); with b3 fixed at 1, the probit and logit
models of glm. Effective doses of
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md) refer to
the curve f, so that b4 keeps its meaning, e.g. the ED50 of the curve;
absolute effective doses take responses on the scale of the response.
`predict`, `fitted`, `residuals` and `simulate` are on the scale of the
response, and the benchmark doses of
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md) and
`bmdpar` are defined on the scale of the response, e.g. the relative
change of the median.
[`peak`](https://daniel-gerhard.github.io/drcTMB/reference/peak.md)
gives the peak dose of the curve and the response at the peak on the
scale of the response. `ed(residual=TRUE)` is not available with a curve
link.

## Small-sample inference

Gaussian and lognormal models without random effects, priors, bmdpar,
correlation structure and variance function are estimated by maximum
likelihood, but by default (`smallsample=TRUE` in
[`drmTMBcontrol`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBcontrol.md))
their inference follows lm, nls and drc: the residual standard deviation
is the square root of the unbiased residual variance, the weighted
residual sum of squares divided by n - p residual degrees of freedom,
for n observations with positive weights and p free coefficients b1-b5;
the covariance of the coefficients (`vcov`, `summary`,
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
`predict`) is scaled accordingly, the Wald and prediction intervals use
t quantiles with n - p degrees of freedom, and the profile intervals of
the coefficients and effective doses an F-based cutoff. The
log-likelihood, the estimate log_sigma in `coef`, AIC, `anova`,
`simulate` and
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
remain maximum likelihood. For the lognormal family, this applies to
log(y). All other models use maximum likelihood and normal quantiles.

## Priors

With `prior`, the estimates are the posterior mode: the objective is the
negative log posterior, the negative log-likelihood (with the random
effects integrated out by the Laplace approximation) minus the log
densities of the priors. A prior for a parameter after the inverse link
(b1-b5) or for a standard deviation or residual parameter, which are
optimised on the log scale, is a density on that scale; the log Jacobian
of the transformation is added, so that the prior is the density of the
parameter as given. Priors with a bounded support (pert, and lognormal
and the half distributions for b1-b5) bound the parameter; a prior whose
support extends beyond the range of the parameter, given by its link and
by `lower` and `upper`, is truncated to the range. Parameters without a
prior have a flat (improper) prior on the scale of their coefficients,
the scale of the link (e.g. on log(b4) with the log link for b4), and
the standard deviations of the random effects and the residual parameter
on the log scale; their bounds still apply. The standard errors and Wald
intervals are from the curvature of the log posterior, and profile
intervals profile the log posterior. `logLik` returns the log-likelihood
at the posterior mode, without the priors, so that AIC and
[`anova.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/anova.drmTMB.md)
are approximate for fits with priors. Small-sample inference is not used
with priors. Refits of
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
keep the priors.

## BMD parameterisation

With `bmdpar`, the model is fitted in the BMD parameterisation, as in
the Bayesian benchmark dose modelling of EFSA (2022): the coefficient of
b4 is replaced by log_bmd, the log benchmark dose of the curve of the
intercepts (the conditional BMD of
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md), with
the random effects at zero), and b4 is computed from it and b1-b3 and
b5. Without priors, this is the same model as the standard one, with the
same likelihood, and exp(log_bmd) is the BMD of
`bmd(fit, bmr, def, scale="conditional")` of the standard fit, by
maximum likelihood (`smallsample` is not used). A prior for the BMD is
given as `prior=list(bmd=...)`, on the dose scale. The standard model is
fitted first, for the starting values and the direction of the curve.
The parameterisation needs a model whose response at dose zero does not
depend on b4 (not logistic, gompertz, the hormesis models, haldane,
constant, or the smoothed threshold of nec), b4 with a single
coefficient and without random effects, and b1-b3 and b5 without
covariates; random effects on b1-b3 and b5 are allowed. For multistage
(with b5 estimated) and nec, the computed b4 keeps its default lower
bound 0. The multistage model of degree 2 (b5 fixed at 0) is the
quadratic exponential model of EFSA (2022): the polynomial is b ((1 - c)
x/dmax + c (x/dmax)^2), with the highest dose dmax, and b3 is replaced
by the technical parameter log_d, with c = exp(-exp(log_d)) in (0, 1);
the rate b follows from the BMD, so that both coefficients of the
polynomial stay positive and a linear curve is the limit c -\> 0, not a
bound. `coef`, `vcov`, `summary`, `logLik`, `confint`, `simulate` and
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
are available for these fits; `ed`, `bmd`, `peak` and `predict` need the
standard parameterisation.

## Correlated residuals and variance functions

With `correlation` or `variance` (gaussian and lognormal families), the
residuals of each series are multivariate normal with a correlation
structure, and their standard deviations sigma g depend on a covariate
or the fitted mean, as in gnls and nlme of nlme (see
[`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md));
random effects can be combined with them. The parameters (cor_par:...,
var_par:...) are estimated with the coefficients and follow log_sigma in
`coef` and `vcov`; `summary` and `confint` give them on their natural
scale, and `anova` compares fits with and without them. Small-sample
inference (`smallsample`) is not used. Weights cannot be combined with a
correlation structure. Pearson residuals are divided by sigma g;
quantile residuals are in addition normalized, decorrelated within
series. Prediction intervals use the standard deviation sigma g at
newdata, which then needs the covariate of the variance function;
`ed(residual=TRUE)` and the BMD definition "sd" are not available with a
variance function.

## Time-to-event data

The event family is for time-to-event data, as `type="event"` of drc:
each row is a time interval \[start, end) of a series (e.g. a Petri
dish) with the number of events in it, e.g. germinated seeds, and a row
with end = Inf holds the units without an event by the end of the
experiment. The curve is the cumulative proportion of events by time t,
f(t) = b2 F(t), with b1 fixed at 0 and the proportion b2 of units with
an event, on the logit scale by default (parlink); the models are those
whose F is a distribution function of the time: loglogistic, weibull1,
weibull2, lognormal, gamma, asymptoticregression and michaelismenten.
The log-likelihood is the multinomial sum(count log(f(end) - f(start))),
with f(Inf) = 1, without the multinomial coefficients, as in drc.
Covariates (fform) and random effects (rform), e.g. of the dishes, are
as for the other families. `ed` gives the times at which a percentage of
the units with an event have had it, and with `type="absolute"` the
times at which a proportion of all units has had an event; `predict`
gives the cumulative proportion at the times in the first variable of
the formula (start). `simulate` and
[`bootstrap`](https://daniel-gerhard.github.io/drcTMB/reference/bootstrap.md)
draw multinomial counts with the total of each series; the residuals
compare the counts with these totals times the probabilities of the
intervals, and `summary` gives their Pearson dispersion. Prediction
intervals and benchmark doses are not available.

## References

Brain, P. and Cousens, R. (1989). An equation to describe dose responses
where there is stimulation of growth at low doses. Weed Research 29,
93-96.

Cedergreen, N., Ritz, C. and Streibig, J. C. (2005). Improved empirical
models describing hormesis. Environmental Toxicology and Chemistry 24,
3166-3172.

EFSA Scientific Committee (2022). Guidance on the use of the benchmark
dose approach in risk assessment. EFSA Journal 20(10), 7584.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. M.
(2016). TMB: automatic differentiation and Laplace approximation.
Journal of Statistical Software 70(5), 1-21.

Pinheiro, J. C. and Bates, D. M. (2000). Mixed-Effects Models in S and
S-PLUS. Springer, New York.

Ritz, C., Baty, F., Streibig, J. C. and Gerhard, D. (2015).
Dose-response analysis using R. PLoS ONE 10(12), e0146021.

## See also

[`drmTMBstart`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBstart.md)
for starting values,
[`drmTMBcontrol`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMBcontrol.md)
for the optimiser,
[`priors`](https://daniel-gerhard.github.io/drcTMB/reference/priors.md),
[`correlation`](https://daniel-gerhard.github.io/drcTMB/reference/correlation.md),
[`drmTMB-methods`](https://daniel-gerhard.github.io/drcTMB/reference/drmTMB-methods.md),
[`ed`](https://daniel-gerhard.github.io/drcTMB/reference/ed.md),
[`bmd`](https://daniel-gerhard.github.io/drcTMB/reference/bmd.md),
[`plot.drmTMB`](https://daniel-gerhard.github.io/drcTMB/reference/plot.drmTMB.md)

## Examples

``` r
# three-parameter log-logistic model (LL.3 of drc): lower asymptote b1 = 0, no asymmetry (b5 = 0)
mod <- drmTMB(rootl ~ conc, data = ryegrass, model = "loglogistic",
              start = c(b1 = 0, b5 = 0), fix = c(b1 = TRUE, b5 = TRUE))
summary(mod)
#> Call:
#> drmTMB(form = rootl ~ conc, data = ryegrass, start = c(b1 = 0, 
#>     b5 = 0), fix = c(b1 = TRUE, b5 = TRUE), model = "loglogistic")
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 24
#> Log-likelihood: -18.6 (df = 4), AIC: 45.21, BIC: 49.92
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b2:(Intercept)    7.855      0.204
#> b3:(Intercept)   -2.470      0.342
#> b4:(Intercept)    3.263      0.196
#> Fixed: b1:(Intercept) = 0, b5:(Intercept) = 0
#> 
#> Residual standard deviation: 0.5616 on 21 degrees of freedom
#> 
#> Optimiser: nlminb: both X-convergence and relative convergence (5), maximum absolute gradient 4.3e-07

# binomial response, separate steepness b3 and ED50 b4 for each time
dmod <- drmTMB(cbind(no, total - no) ~ dose, data = daphnids, model = "loglogistic",
               family = "binomial", fform = b3 + b4 ~ 0 + time,
               start = c(b1 = 0, b2 = 1, b4 = 1000, b5 = 0),
               fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE))
dmod
#> Call:
#> drmTMB(form = cbind(no, total - no) ~ dose, fform = b3 + b4 ~ 
#>     0 + time, data = daphnids, start = c(b1 = 0, b2 = 1, b4 = 1000, 
#>     b5 = 0), fix = c(b1 = TRUE, b2 = TRUE, b5 = TRUE), family = "binomial", 
#>     model = "loglogistic")
#> 
#> Model:  loglogistic
#> Family: binomial, identity link for the asymptotes
#> 
#> Coefficients:
#> b3:time24h  b3:time48h  b4:time24h  b4:time48h  
#>      1.174       1.850    5134.036    1509.065  
#> Fixed: b1:(Intercept) = 0, b2:(Intercept) = 1, b5:(Intercept) = 0

# correlated random effects of the experiments on both asymptotes, and a power variance function
vmod <- drmTMB(effect ~ conc, data = vinclozolin, model = "loglogistic",
               rform = b1 + b2 ~ 1 | exper, variance = var_power(),
               start = c(b5 = 0), fix = c(b5 = TRUE),
               lower = c(0, 0, -Inf, 0, -Inf))  # non-negative asymptotes b1, b2 and ED50 b4
summary(vmod)
#> Call:
#> drmTMB(form = effect ~ conc, rform = b1 + b2 ~ 1 | exper, data = vinclozolin, 
#>     start = c(b5 = 0), fix = c(b5 = TRUE), lower = c(0, 0, -Inf, 
#>         0, -Inf), model = "loglogistic", variance = var_power())
#> 
#> Model:  loglogistic, decreasing curve
#> Family: gaussian, identity link for the asymptotes
#> 
#> Number of observations: 53, groups: exper 6
#> Log-likelihood: -335.9 (df = 9), AIC: 689.9, BIC: 707.6
#> 
#> Coefficients:
#>                  Estimate Std. Error
#> b1:(Intercept)  252.72312     78.637
#> b2:(Intercept) 1924.39964    258.759
#> b3:(Intercept)   -0.72717      0.103
#> b4:(Intercept)    0.05837      0.013
#> Fixed: b5:(Intercept) = 0
#> 
#> Random effects:
#>  Groups Name Std.Dev. Corr
#>  exper  b1   156.8        
#>         b2   566.6    0.96
#> 
#> Residual standard deviation: 0.2197
#> 
#> Variance function: power, ~fitted: delta = 0.9453
#> 
#> Optimiser: nlminb: relative convergence (4), maximum absolute gradient 0.00017
```
