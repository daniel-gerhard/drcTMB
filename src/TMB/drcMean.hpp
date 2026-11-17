#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// sum of the means of the curves at the doses x, for the standard errors of predict: each mean depends
// only on the linear predictors in its own row, so gr() gives the derivatives of every mean w.r.t. its
// linear predictors of b1-b5 in a single sweep
template<class Type>
Type drcMean(objective_function<Type>* obj) {
  DATA_INTEGER(mod);
  DATA_IVECTOR(lnk);     // links of b1-b5: 1 = identity, 2 = log, 3 = logit
  DATA_INTEGER(clink);   // link of the curve (curvelink)
  DATA_SCALAR(alpha);    // exponent of the hormesis term of the cedergreen model
  DATA_SCALAR(smooth);   // smoothing width of the threshold of the nec model
  DATA_VECTOR(x);
  PARAMETER_MATRIX(eta); // linear predictors of b1-b5: a row for each dose

  Type s = 0;
  for (int i = 0; i < x.size(); i++){
    vector<Type> e = eta.row(i);
    s += curveat(mod, clink, Type(x(i)), curvepars(lnk, e, false, Type(0), false, Type(0)), alpha, smooth);
  }
  return s;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
