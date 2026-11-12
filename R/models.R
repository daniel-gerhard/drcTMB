# starting values of b3 and b4 from F linearised by z: z(F) = b3 (x - b4), or b3 (log x - log b4) with
# the positive doses, by least squares; NA without two distinct doses
linstart <- function(z, logdose){
  function(x, q){
    if (logdose){
      q <- q[x > 0]
      x <- log(x[x > 0])
    }
    if (length(unique(x)) < 2) return(c(NA, NA))
    cf <- coefficients(lm(z(q) ~ x))
    c(cf[2], if (logdose) exp(-cf[1]/cf[2]) else -cf[1]/cf[2])
  }
}

# Properties of the dose-response models. The position of a model in the list is its index (mod) in the
# TMB templates, which define the mean function (drcmean) and the effective dose (drced) of each model;
# F and edh below are the same functions in R. A new model needs an entry here, at the same index as in
# the templates.
#
# b1-b5: role of each parameter, NA if the model does not use it (it is then always fixed)
#   b1, b2: "asymptote", or "mean" for the constant model; the lower and upper asymptote for a steepness
#           b3, otherwise the asymptotes at dose zero and at high doses; b1 "intercept": the response at
#           dose zero of a model without asymptotes
#   b3: "steepness" (its sign gives the direction of the curve), "shape", "rate", "power", or "coefficient" of a
#       polynomial
#   b4: "location" on the dose scale, "scale" (a positive dose, the effective dose is b4 h), "threshold" (a
#       positive dose, the effective dose is b4 + h), or "coefficient" of a polynomial or of b4 x^b3
#   b5: "asymmetry" or "shape" (both positive, with the log link), "hormesis" effect, "inhibition", or
#       "coefficient" of a polynomial
# The functions F and edh take b3-b5 on their original scale, after the inverse links (see modellinks).
# logdose: F depends on the dose through t = b3 (log x - log b4) instead of t = b3 (x - b4)
# F(t, x, b3, b4, b5, smooth): increasing from 0 to 1; NULL for the constant model
# ed: how effective doses are computed: "closed" by edh, "hormesis" numerically above the peak,
#   "polynomial" as the root of the polynomial, "unbounded" absolute only (see mean, level); NA for none
# edh(q, b3, b5, smooth): effective dose at level q of F, b4 h, or b4 + h for a location or threshold b4
# start: starting values in lmstart: "linear" from a linearised F (linear), "asymptotic" from the residual
#   sum of squares profiled over b4, "inhibition" profiled over b4 and b5, "threshold", "polynomial",
#   "constant", "power" (profiled over b3, linear in b1 and b4)
# linear(x, q): starting values of b3 and b4 from the doses x and the values q of F
# lower0: parameters bounded below by 0 by default, if they have a single coefficient (b4 only with the
#   identity link)
# alpha: the model needs the exponent alpha of its hormesis term
# hormesis: hormesis term g(x, alpha) of f1 + (f2 - f1 + f5 g) F, its first and second derivatives dg and
#   d2g w.r.t. the log dose, and the largest b3 of a decreasing curve that the model is meant for
# mean, level: the models without an upper asymptote (exponential, power), f(x) = mean(x, b1, b3, b4) =
#   g^-1(b1, b4 x^b3) with b1 the response at dose zero, b3 > 0 the power and b4 the coefficient (any
#   sign, the direction of the curve); level(y, b1) = g(y, b1) = b4 x^b3 at the dose x with response y,
#   which replaces the level of F, (y - b1)/(b2 - b1), of the other models; ed "unbounded": only absolute
#   effective doses, (level/b4)^(1/b3)
# peak: the dose of the peak of the curve in closed form, dose(b) from the parameters b1-b5, NA without a
#   peak, and its inverse b4(pk, b), the b4 of a curve with its peak at dose pk (as in the TMB template)
drcmodels <- list(
  logistic=list(b1="asymptote", b2="asymptote", b3="steepness", b4="location", b5="asymmetry", logdose=FALSE,
                F=function(t, x, b3, b4, b5, smooth) (1 + exp(-t))^(-b5),
                ed="closed", edh=function(q, b3, b5, smooth) -log(q^(-1/b5) - 1)/b3,
                start="linear", linear=linstart(function(q) qlogis(q), FALSE)),
  loglogistic=list(b1="asymptote", b2="asymptote", b3="steepness", b4="scale", b5="asymmetry", logdose=TRUE,
                   F=function(t, x, b3, b4, b5, smooth) (1 + exp(-t))^(-b5),
                   ed="closed", edh=function(q, b3, b5, smooth) (q^(-1/b5) - 1)^(-1/b3),
                   start="linear", linear=linstart(function(q) qlogis(q), TRUE)),
  weibull1=list(b1="asymptote", b2="asymptote", b3="steepness", b4="scale", b5=NA, logdose=TRUE,
                F=function(t, x, b3, b4, b5, smooth) exp(-exp(-t)),
                ed="closed", edh=function(q, b3, b5, smooth) (-log(q))^(-1/b3),
                start="linear", linear=linstart(function(q) -log(-log(q)), TRUE)),
  weibull2=list(b1="asymptote", b2="asymptote", b3="steepness", b4="scale", b5=NA, logdose=TRUE,
                F=function(t, x, b3, b4, b5, smooth) 1 - exp(-exp(t)),
                ed="closed", edh=function(q, b3, b5, smooth) (-log(1 - q))^(1/b3),
                start="linear", linear=linstart(function(q) log(-log(1 - q)), TRUE)),
  lognormal=list(b1="asymptote", b2="asymptote", b3="steepness", b4="scale", b5=NA, logdose=TRUE,
                 F=function(t, x, b3, b4, b5, smooth) pnorm(t),
                 ed="closed", edh=function(q, b3, b5, smooth) exp(qnorm(q)/b3),
                 start="linear", linear=linstart(function(q) qnorm(q), TRUE)),
  michaelismenten=list(b1="asymptote", b2="asymptote", b3=NA, b4="scale", b5=NA, logdose=FALSE,
                       F=function(t, x, b3, b4, b5, smooth) 1 / (1 + b4 / x),
                       ed="closed", edh=function(q, b3, b5, smooth) q/(1 - q),
                       # 1/F - 1 = b4/x, linear in 1/x through the origin
                       start="linear", linear=function(x, q){
                         q <- q[x > 0]
                         x <- x[x > 0]
                         if (length(unique(x)) < 2) return(c(NA, NA))
                         c(NA, coefficients(lm(I(1/q - 1) ~ 0 + I(1/x)))[1])
                       }),
  asymptoticregression=list(b1="asymptote", b2="asymptote", b3=NA, b4="scale", b5=NA, logdose=FALSE,
                            F=function(t, x, b3, b4, b5, smooth) 1 - exp(-x / b4),
                            ed="closed", edh=function(q, b3, b5, smooth) -log(1 - q),
                            start="asymptotic"),
  exponentialdecay=list(b1="asymptote", b2="asymptote", b3=NA, b4="scale", b5=NA, logdose=FALSE,
                        F=function(t, x, b3, b4, b5, smooth) 1 - exp(-x / b4),
                        ed="closed", edh=function(q, b3, b5, smooth) -log(1 - q),
                        start="asymptotic"),
  gompertz=list(b1="asymptote", b2="asymptote", b3="steepness", b4="location", b5=NA, logdose=FALSE,
                F=function(t, x, b3, b4, b5, smooth) exp(-exp(-t)),
                ed="closed", edh=function(q, b3, b5, smooth) -log(-log(q))/b3,
                start="linear", linear=linstart(function(q) -log(-log(q)), FALSE)),
  constant=list(b1="mean", b2=NA, b3=NA, b4=NA, b5=NA, logdose=FALSE, ed=NA, start="constant"),
  braincousens=list(b1="asymptote", b2="asymptote", b3="steepness", b4="scale", b5="hormesis", logdose=TRUE,
                    F=function(t, x, b3, b4, b5, smooth) plogis(t),
                    ed="hormesis", start="linear", linear=linstart(function(q) qlogis(q), TRUE),
                    hormesis=list(g=function(x, alpha) x, dg=function(x, alpha) x, d2g=function(x, alpha) x,
                                  b3max=-1)),
  cedergreen=list(b1="asymptote", b2="asymptote", b3="steepness", b4="scale", b5="hormesis", logdose=TRUE,
                  F=function(t, x, b3, b4, b5, smooth) plogis(t),
                  ed="hormesis", start="linear", linear=linstart(function(q) qlogis(q), TRUE), alpha=TRUE,
                  hormesis=list(g=function(x, alpha) exp(-x^(-alpha)),
                                dg=function(x, alpha) alpha * x^(-alpha) * exp(-x^(-alpha)),
                                d2g=function(x, alpha){
                                  u <- x^(-alpha)
                                  alpha^2 * u * exp(-u) * (u - 1)
                                },
                                b3max=0)),
  gamma=list(b1="asymptote", b2="asymptote", b3="shape", b4="scale", b5=NA, logdose=FALSE,
             F=function(t, x, b3, b4, b5, smooth) pgamma(x / b4, b3),
             ed="closed", edh=function(q, b3, b5, smooth) qgamma(q, b3),
             start="asymptotic"),
  multistage=list(b1="asymptote", b2="asymptote", b3="coefficient", b4="coefficient", b5="coefficient", logdose=FALSE,
                  F=function(t, x, b3, b4, b5, smooth) 1 - exp(-(b3*x + b4*x^2 + b5*x^3)),
                  ed="polynomial", start="polynomial", lower0=3:5),
  nec=list(b1="asymptote", b2="asymptote", b3="rate", b4="threshold", b5="shape", logdose=FALSE,
           F=function(t, x, b3, b4, b5, smooth) 1 - exp(-(b3*necdose(x - b4, smooth))^b5),
           # with smoothing width smooth > 0, the inverse of the smoothing of the dose above the threshold
           ed="closed", edh=function(q, b3, b5, smooth){
             h <- (-log(1 - q))^(1/b5)/b3
             if (smooth > 0) h + smooth*log(-expm1(-h/smooth)) else h
           },
           start="threshold", lower0=3:4),
  haldane=list(b1="asymptote", b2="asymptote", b3=NA, b4="scale", b5="inhibition", logdose=FALSE,
               F=function(t, x, b3, b4, b5, smooth) x / (b4 + x + b5*x^2),
               ed=NA, start="inhibition", lower0=5,
               peak=list(dose=function(b) if (isTRUE(b[5] > 0 && b[4] > 0)) sqrt(b[4]/b[5]) else NA,
                         b4=function(pk, b) pk^2*b[5])),
  exponential=list(b1="intercept", b2=NA, b3="power", b4="coefficient", b5=NA, logdose=FALSE,
                   mean=function(x, f1, f3, f4) f1 * exp(f4 * x^f3), level=function(y, f1) log(y / f1),
                   ed="unbounded", start="power", lower0=3),
  power=list(b1="intercept", b2=NA, b3="power", b4="coefficient", b5=NA, logdose=FALSE,
             mean=function(x, f1, f3, f4) f1 + f4 * x^f3, level=function(y, f1) y - f1,
             ed="unbounded", start="power", lower0=3)
)

# links of the parameters b1-b5 of a model m (see drcmodels): the possible links of each parameter, the
# first one is the default; a positive b4 may have the log link, and a positive b5 has the log link
modellinks <- function(m){
  list(b1=c("identity", "log", "logit"),
       b2=if (is.na(m$b2)) "identity" else c("identity", "log", "logit"),
       b3="identity",
       b4=if (m$b4 %in% c("scale", "threshold")) c("identity", "log") else "identity",
       b5=if (m$b5 %in% c("asymmetry", "shape")) "log" else "identity")
}

# index of each link, as in the TMB templates: 1 = identity, 2 = log, 3 = logit
linkindex <- function(link){
  match(link, c("identity", "log", "logit"))
}

# parameters b1-b5 from their linear predictors eta (a vector, or a matrix with a column for each
# parameter), by the inverse links lnk (see linkindex)
linkpars <- function(lnk, eta){
  if (is.matrix(eta)) return(do.call(cbind, lapply(1:5, function(k) linkfun(lnk[k], eta[, k], inverse=TRUE))))
  sapply(1:5, function(k) linkfun(lnk[k], eta[k], inverse=TRUE))
}

# derivative of the inverse link lnk at the linear predictor eta
dlinkinv <- function(lnk, eta){
  switch(lnk, rep(1, length(eta)), exp(eta), plogis(eta)*(1 - plogis(eta)))
}

# properties of a model, by name or index
drcmodel <- function(mod){
  drcmodels[[mod]]
}

# index of the model, as in the TMB template
modelindex <- function(model){
  which(model == names(drcmodels))
}

# names of the models with a property, for messages
modelswith <- function(cond){
  nm <- names(Filter(cond, drcmodels))
  if (length(nm) == 1) return(nm)
  paste(paste(nm[-length(nm)], collapse=", "), "and", nm[length(nm)])
}
