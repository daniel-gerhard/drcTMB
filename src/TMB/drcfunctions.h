// Helper functions shared by the drcTMB objective function.
// Kept in a .h file (not .hpp) so TMBtools does not treat it as a model.

#ifndef DRCFUNCTIONS_H
#define DRCFUNCTIONS_H

// inverse link for the asymptote parameters b1 and b2
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
template<class Type>
Type drcmean(int mod, Type x, Type f1, Type f2, Type f3, Type f4, Type f5){
  Type fl;
  Type f = 0;
  switch (mod){
  case 1:
    fl = 1 + exp(f3*(x - f4));
    f = f2 + (f1 - f2) / (pow(fl, exp(f5)));
    break;
  case 2:
    if (x == 0){
      if (f3 > 0){
        f = f2 + (f1 - f2);
      } else {
        f = f2;
      }
    } else {
      fl = 1 + exp(f3*(log(x) - log(f4)));
      f = f2 + (f1 - f2) / (pow(fl, exp(f5)));
    }
    break;
  case 3:
    if (x == 0){
      if (f3 > 0){
        f = f2 + (f1 - f2);
      } else {
        f = f2;
      }
    } else {
      fl = exp(-1*exp(f3*(log(x) - log(f4))));
      f = f2 + (f1 - f2) * fl;
    }
    break;
  case 4:
    if (x == 0){
      if (f3 < 0){
        f = f2 + (f1 - f2);
      } else {
        f = f2;
      }
    } else {
      fl = 1 - exp(-1*exp(f3*(log(x) - log(f4))));
      f = f2 + (f1 - f2) * fl;
    }
    break;
  case 5:
    if (x == 0){
      if (f3 < 0){
        f = f2 + (f1 - f2);
      } else {
        f = f2;
      }
    } else {
      fl = f3*(log(x) - log(f4));
      f = f2 + (f1 - f2) * pnorm(fl);
    }
    break;
  case 6:
    fl = 1 + f4 / x;
    f = f2 + (f1 - f2) / fl;
    break;
  case 7:
    f = f2 + (f1 - f2) / (1 - exp(-x / f4));
    break;
  case 8:
    f = f2 + (f1 - f2) / (exp(-x / f4));
    break;
  case 9:
    f = f2 + (f1 - f2) * (exp(-exp(f3*(x - f4))));
    break;
  }
  return f;
}

// effective dose: inverse of drcmean at relative response level p
// r: shift on the response scale (residual), q: relative position between f2 and f1
template<class Type>
Type drced(int mod, Type p, Type f1, Type f2, Type f3, Type f4, Type f5, Type r){
  Type q = p - r / (f1 - f2);
  Type ed = 0;
  switch (mod){
  case 1:
    ed = f4 + 1/f3 * log(pow(q, -1/exp(f5)) - 1);
    break;
  case 2:
    ed = f4 * pow(pow(q, -1/exp(f5)) - 1, 1/f3);
    break;
  case 3:
    ed = f4 * pow(-1*log(q), 1/f3);
    break;
  case 4:
    ed = f4 * pow(-1*log(1 - q), -1/f3);
    break;
  case 5:
    ed = f4 * exp(-1/f3 * qnorm(q));
    break;
  case 6:
    ed = f4 * q / (1 - q);
    break;
  case 7:
    ed = -1*f4 * log(1 - 1/q);
    break;
  case 8:
    ed = f4 * log(q);
    break;
  case 9:
    ed = f4 + 1/f3 * log(-1*log(q));
    break;
  }
  return ed;
}

// log-likelihood contribution of a single observation
// fam: 1 = gaussian (phi = sd), 2 = poisson, 3 = binomial (bn = number of trials),
//      4 = beta (phi = precision)
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
  default:
    return 0;
  }
}

#endif
