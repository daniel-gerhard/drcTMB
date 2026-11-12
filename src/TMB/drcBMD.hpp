#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// response at the benchmark dose: the response mu0 at dose zero changed by the benchmark response bmr in
// the direction inc of the curve (1 increasing, -1 decreasing)
// def: 1 = excess (extra) risk, 2 = additional (added) risk, 3 = relative deviation, 4 = standard deviations
template<class Type>
Type bmdtarget(int def, Type bmr, int inc, Type mu0, Type sigma){
  switch (def){
  case 1:
    return inc > 0 ? Type(mu0 + bmr * (1 - mu0)) : Type(mu0 * (1 - bmr));
  case 2:
    return mu0 + Type(inc) * bmr;
  case 3:
    return mu0 * (1 + Type(inc) * bmr);
  case 4:
    return mu0 + Type(inc) * bmr * sigma;
  case 5:
    // lognormal: bmr standard deviations of log(y)
    return mu0 * exp(Type(inc) * bmr * sigma);
  default:
    return mu0;
  }
}

// marginal mean at dose x: the mean over the random effects, by quadrature with the linear predictors
// eta0 + Z.col(k) of b1-b5 and weights wq(k)
template<class Type>
Type margmean(int mod, vector<int> lnk, Type x, vector<Type> eta0, matrix<Type> Z, vector<Type> wq, Type alpha, Type smooth){
  Type m = 0;
  for (int k = 0; k < wq.size(); k++){
    m += wq(k) * drcmean(mod, x,
                         linkinv(lnk(0), Type(eta0(0) + Z(0, k))),
                         linkinv(lnk(1), Type(eta0(1) + Z(1, k))),
                         linkinv(lnk(2), Type(eta0(2) + Z(2, k))),
                         linkinv(lnk(3), Type(eta0(3) + Z(3, k))),
                         linkinv(lnk(4), Type(eta0(4) + Z(4, k))), alpha, smooth);
  }
  return m;
}

// benchmark dose for use with MakeADFun(..., ADreport=TRUE): fn() returns the BMD, gr() its gradient
// w.r.t. (b1-b5, u1-u5, theta, log_sigma)
template<class Type>
Type drcBMD(objective_function<Type>* obj) {
  DATA_INTEGER(mod);
  DATA_IVECTOR(lnk);    // links of b1-b5: 1 = identity, 2 = log, 3 = logit
  DATA_SCALAR(alpha);
  DATA_SCALAR(smooth);
  DATA_INTEGER(def);    // definition of the benchmark response, see bmdtarget
  DATA_SCALAR(bmr);
  DATA_INTEGER(inc);    // direction of the curve: 1 increasing, -1 decreasing
  DATA_SCALAR(sfac);    // residual standard deviation over exp(log_sigma): 1, or the small-sample factor
  DATA_SCALAR(x0);      // multistage: the root of the polynomial; marginal BMD: the root of the marginal mean, found in R
  DATA_IVECTOR(tnp);    // random effect terms: number of random coefficients per term
  DATA_MATRIX(A);       // loads of the random coefficients (columns) on b1-b5 (rows); zero for terms not integrated
  DATA_MATRIX(E);       // standard normal quadrature nodes, a row for each random coefficient, a column for each node
  DATA_VECTOR(wq);      // quadrature weights; length zero for the conditional BMD (no integration)
  PARAMETER_VECTOR(b);  // b1-b5
  PARAMETER_VECTOR(u);  // random effects of a fitted group on b1-b5
  PARAMETER_VECTOR(theta); // standard deviations and correlations of the random effects, as in drcTMB
  PARAMETER(log_sigma);

  vector<Type> eta0 = b + u;
  Type sigma = sfac * exp(log_sigma);
  Type bmd;

  if (wq.size() == 0){
    // BMD of the curve with linear predictors eta0: the effective dose at the level of F of the target
    Type f1 = linkinv(lnk(0), eta0(0));
    Type f2 = linkinv(lnk(1), eta0(1));
    Type f3 = linkinv(lnk(2), eta0(2));
    Type f4 = linkinv(lnk(3), eta0(3));
    Type f5 = linkinv(lnk(4), eta0(4));
    Type mu0 = drcmean(mod, Type(0), f1, f2, f3, f4, f5, alpha, smooth);
    Type qf = drclevel(mod, bmdtarget(def, bmr, inc, mu0, sigma), f1, f2);
    bmd = mod == 14 ? msed(qf, f1, f2, f3, f4, f5, Type(0), Type(x0)) : drced(mod, qf, f1, f2, f3, f4, f5, Type(0), smooth);
  } else {
    // Cholesky factor C of the block-diagonal covariance of the random coefficients: for each term
    // diag(sd) times L with rows of unit length, where the correlation is L L' as UNSTRUCTURED_CORR
    int nc = A.cols();
    matrix<Type> C(nc, nc);
    C.setZero();
    int toff = 0, coff = 0;
    for (int t = 0; t < tnp.size(); t++){
      int np = tnp(t);
      matrix<Type> L(np, np);
      L.setIdentity();
      int k = 0;
      for (int i = 0; i < np; i++){
        for (int j = 0; j < i; j++){
          L(i, j) = theta(toff + np + k);
          k++;
        }
      }
      for (int i = 0; i < np; i++){
        Type s = 0;
        for (int j = 0; j <= i; j++) s += L(i, j) * L(i, j);
        Type sdi = exp(theta(toff + i));
        for (int j = 0; j <= i; j++) C(coff + i, coff + j) = sdi * L(i, j) / sqrt(s);
      }
      toff += np + np*(np-1)/2;
      coff += np;
    }
    matrix<Type> Z = A * C * E;
    // the marginal BMD: one Newton step on the marginal mean from its root x0, which keeps its value and
    // gives the derivatives of the implicit function; the derivative w.r.t. the dose by central differences
    Type mu0 = margmean(mod, lnk, Type(0), eta0, Z, wq, alpha, smooth);
    Type target = bmdtarget(def, bmr, inc, mu0, sigma);
    Type xs = x0;
    Type h = Type(1e-5) * xs;
    Type G = margmean(mod, lnk, xs, eta0, Z, wq, alpha, smooth) - target;
    Type dG = (margmean(mod, lnk, Type(xs + h), eta0, Z, wq, alpha, smooth) -
               margmean(mod, lnk, Type(xs - h), eta0, Z, wq, alpha, smooth)) / (2 * h);
    bmd = xs - G / dG;
  }
  ADREPORT(bmd);
  return 0;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
