// Helper functions shared by the drcTMB objective function.
// Kept in a .h file (not .hpp) so TMBtools does not treat it as a model.

#ifndef DRCFUNCTIONS_H
#define DRCFUNCTIONS_H

// inverse link for the asymptote parameters b1 and b2, and for b4
// lnk: 1 = identity, 2 = log, 3 = logit
template<class Type>
Type linkinv(int lnk, Type eta){
  switch (lnk){
  case 2:
    return exp(eta);
  case 3:
    return 1/(1 + exp(-1*eta));
  default:
    return eta;
  }
}

// dose-response mean function
// mod: 1 = logistic, 2 = loglogistic, 3 = weibull1, 4 = weibull2, 5 = lognormal,
//      6 = michaelismenten, 7 = asymptoticregression, 8 = exponentialdecay, 9 = gompertz
// f(x) = f1 + (f2 - f1) * F(x), with F increasing from 0 to 1 in its argument.
// Models with steepness f3: F(f3 * z), z = x - f4 or log(x) - log(f4); f1 is the lower and
// f2 the upper asymptote, f3 > 0 for an increasing and f3 < 0 for a decreasing curve.
// michaelismenten, asymptoticregression, exponentialdecay: F(x) increases with dose, f1 is
// the asymptote at dose zero and f2 the asymptote at high doses.
template<class Type>
Type drcmean(int mod, Type x, Type f1, Type f2, Type f3, Type f4, Type f5){
  bool logdose = (mod >= 2 && mod <= 5);
  // log-dose models at dose zero: the asymptote approached at low doses
  if (logdose && x == 0){
    if (f3 > 0){
      return f1;
    } else {
      return f2;
    }
  }
  if (mod == 6 && x == 0) return f1;

  Type t = f3 * (logdose ? Type(log(x) - log(f4)) : Type(x - f4));
  Type F = 0;
  switch (mod){
  case 1:
  case 2:
    // (1 + exp(-t))^(-exp(f5)), on the log scale for numerical stability
    F = exp(-exp(f5) * logspace_add(Type(0), Type(-t)));
    break;
  case 3:
  case 9:
    F = exp(-exp(-t));
    break;
  case 4:
    F = 1 - exp(-exp(t));
    break;
  case 5:
    F = pnorm(t);
    break;
  case 6:
    F = 1 / (1 + f4 / x);
    break;
  case 7:
  case 8:
    F = 1 - exp(-x / f4);
    break;
  }
  return f1 + (f2 - f1) * F;
}

// effective dose: inverse of drcmean, the dose at which F equals q
// r: shift on the response scale (residual)
template<class Type>
Type drced(int mod, Type q, Type f1, Type f2, Type f3, Type f4, Type f5, Type r){
  Type qr = q - r / (f2 - f1);
  Type ed = 0;
  switch (mod){
  case 1:
    ed = f4 - log(pow(qr, -1/exp(f5)) - 1) / f3;
    break;
  case 2:
    ed = f4 * pow(pow(qr, -1/exp(f5)) - 1, -1/f3);
    break;
  case 3:
    ed = f4 * pow(-1*log(qr), -1/f3);
    break;
  case 4:
    ed = f4 * pow(-1*log(1 - qr), 1/f3);
    break;
  case 5:
    ed = f4 * exp(qnorm(qr) / f3);
    break;
  case 6:
    ed = f4 * qr / (1 - qr);
    break;
  case 7:
  case 8:
    ed = -1*f4 * log(1 - qr);
    break;
  case 9:
    ed = f4 - log(-1*log(qr)) / f3;
    break;
  }
  return ed;
}

// log-likelihood contribution of a single observation
// fam: 1 = gaussian (phi = sd), 2 = poisson, 3 = binomial (bn = number of trials),
//      4 = beta (phi = precision), 5 = nbinom2 (variance f (1 + f/phi), phi = size),
//      6 = nbinom1 (variance f (1 + phi))
template<class Type>
Type loglik(int fam, Type y, Type f, Type bn, Type phi){
  switch (fam){
  case 1:
    return dnorm(y, f, phi, true);
  case 2:
    return dpois(y, f, true);
  case 3:
    return dbinom(y, bn, f, true);
  case 4:
    return lgamma(phi) - lgamma(f*phi) - lgamma((1-f)*phi) + (f*phi-1)*log(y) + ((1-f)*phi - 1)*log(1 - y);
  case 5:
    // on the log scale of the mean and of the variance minus the mean, for numerical stability
    return dnbinom_robust(y, log(f), Type(2*log(f) - log(phi)), 1);
  case 6:
    return dnbinom_robust(y, log(f), Type(log(f) + log(phi)), 1);
  default:
    return 0;
  }
}

// simulate a single observation with mean f, same parameterisation as loglik
template<class Type>
Type simobs(int fam, Type f, Type bn, Type phi){
  switch (fam){
  case 1:
    return rnorm(f, phi);
  case 2:
    return rpois(f);
  case 3:
    return rbinom(bn, f);
  case 4:
    return rbeta(f*phi, (1-f)*phi);
  case 5:
    return rnbinom2(f, Type(f*(1 + f/phi)));
  case 6:
    return rnbinom2(f, Type(f*(1 + phi)));
  default:
    return 0;
  }
}

#endif
