#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// effective dose for use with MakeADFun(..., ADreport=TRUE):
// fn() returns the ED, gr() its gradient w.r.t. (b1-b5, u1-u5, r)
template<class Type>
Type drcED(objective_function<Type>* obj) {
  DATA_INTEGER(mod);
  DATA_INTEGER(lnk);
  DATA_INTEGER(lnk4);   // link for b4: 1 = identity, 2 = log
  DATA_INTEGER(absolute); // 0: q is the level of F in drcmean, 1: q is the response at the effective dose
  DATA_SCALAR(q);
  PARAMETER_VECTOR(b);  // b1-b5
  PARAMETER_VECTOR(u);  // random effects u1-u5, evaluated at zero
  PARAMETER(r);         // residual, evaluated at zero

  vector<Type> f = b + u;
  Type f1 = linkinv(lnk, f(0));
  Type f2 = linkinv(lnk, f(1));
  // level of F at which the curve reaches the response q
  Type qf = absolute == 1 ? Type((q - f1) / (f2 - f1)) : q;
  Type ed = drced(mod, qf, f1, f2, f(2), linkinv(lnk4, f(3)), f(4), r);
  ADREPORT(ed);
  return 0;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
