#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// effective dose for use with MakeADFun(..., ADreport=TRUE):
// fn() returns the ED, gr() its gradient w.r.t. (b1-b5, u1-u5, r)
template<class Type>
Type drcED(objective_function<Type>* obj) {
  DATA_INTEGER(mod);
  DATA_INTEGER(lnk);
  DATA_SCALAR(q);       // level of F in drcmean at the effective dose
  PARAMETER_VECTOR(b);  // b1-b5
  PARAMETER_VECTOR(u);  // random effects u1-u5, evaluated at zero
  PARAMETER(r);         // residual, evaluated at zero

  vector<Type> f = b + u;
  Type ed = drced(mod, q, linkinv(lnk, f(0)), linkinv(lnk, f(1)), f(2), f(3), f(4), r);
  ADREPORT(ed);
  return 0;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
