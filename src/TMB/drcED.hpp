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
  DATA_SCALAR(alpha);   // exponent of the hormesis term of the cedergreen model
  DATA_SCALAR(smooth);  // smoothing width of the threshold of the nec model
  DATA_SCALAR(s0);      // hormesis and multistage models: log dose of the effective dose or of the peak, found in R
  DATA_INTEGER(target); // 0: effective dose, hormesis models: 1 = dose of the peak, 2 = response at the peak
  PARAMETER_VECTOR(b);  // b1-b5
  PARAMETER_VECTOR(u);  // random effects u1-u5, evaluated at zero
  PARAMETER(r);         // residual, evaluated at zero

  vector<Type> f = b + u;
  Type f1 = linkinv(lnk, f(0));
  Type f2 = linkinv(lnk, f(1));
  // level of F at which the curve reaches the response q
  Type qf = absolute == 1 ? Type((q - f1) / (f2 - f1)) : q;
  Type f4 = linkinv(lnk4, f(3));
  Type ed;
  if (target > 0){
    Type x = exp(hormpeak(mod, f1, f2, f(2), f4, f(4), alpha, s0));
    ed = target == 1 ? x : drcmean(mod, x, f1, f2, f(2), f4, f(4), alpha, smooth);
  } else if (mod == 11 || mod == 12){
    ed = hormed(mod, qf, f1, f2, f(2), f4, f(4), r, alpha, s0);
  } else if (mod == 14){
    ed = msed(qf, f1, f2, f(2), f4, f(4), r, Type(exp(s0)));
  } else {
    ed = drced(mod, qf, f1, f2, f(2), f4, f(4), r, smooth);
  }
  ADREPORT(ed);
  return 0;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
