#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type drcTMB(objective_function<Type>* obj) {
  DATA_INTEGER(fam);
  DATA_INTEGER(mod);
  DATA_INTEGER(lnk);
  DATA_IVECTOR(ind);       // 0/1 for b1-b5: parameter has a random effect
  DATA_INTEGER(simre);     // simulate: 1 = new random effects, 0 = keep the fitted ones
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_VECTOR(w);
  DATA_VECTOR(bn);
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  DATA_MATRIX(X3);
  DATA_MATRIX(X4);
  DATA_MATRIX(X5);
  DATA_SPARSE_MATRIX(Z);   // zero columns without random effects
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(b3);
  PARAMETER_VECTOR(b4);
  PARAMETER_VECTOR(b5);
  PARAMETER_VECTOR(u);     // length zero without random effects
  PARAMETER_VECTOR(theta); // length zero without random effects
  PARAMETER(log_sigma);    // mapped (fixed) for poisson and binomial

  Type phi = exp(log_sigma);
  if (fam == 1){
    Type sigma = phi;
    ADREPORT(sigma);
  }
  if (fam == 4){
    ADREPORT(phi);
  }

  Type nll = 0;
  int n = y.size();

  // random effects: column k of Zu is added to the linear predictor of b(k+1)
  matrix<Type> Zu(n, 5);
  Zu.setZero();
  int rn = ind.sum();
  if (rn > 0){
    int un = u.size()/rn;
    vector<int> dim(2);
    dim << un, rn;
    array<Type> U(u, dim);

    vector<Type> logsd = theta.head(rn);
    vector<Type> corr_transf = theta.tail(theta.size() - rn);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);

    for (int j = 0; j < un; j++){
      nll += scnldens(U.transpose().col(j));
    }

    // new random effects for each cluster
    SIMULATE {
      if (simre == 1){
        vector<Type> uj(rn);
        for (int j = 0; j < un; j++){
          scnldens.simulate(uj);
          for (int k = 0; k < rn; k++) U(j, k) = uj(k);
        }
      }
    }

    matrix<Type> corr = nldens.cov();
    ADREPORT(corr);

    int ucount = 0;
    for (int k = 0; k < 5; k++){
      if (ind(k) > 0){
        vector<Type> Zuk = Z * U.col(ucount).matrix();
        Zu.col(k) = Zuk.matrix();
        ucount = ucount + 1;
      }
    }
  }

  vector<Type> Xb1 = X1 * b1;
  vector<Type> Xb2 = X2 * b2;
  vector<Type> Xb3 = X3 * b3;
  vector<Type> Xb4 = X4 * b4;
  vector<Type> Xb5 = X5 * b5;

  Type f;
  for (int i = 0; i < n; i++){
    f = drcmean(mod, x(i),
                linkinv(lnk, Type(Xb1(i) + Zu(i, 0))),
                linkinv(lnk, Type(Xb2(i) + Zu(i, 1))),
                Type(Xb3(i) + Zu(i, 2)),
                Type(Xb4(i) + Zu(i, 3)),
                Type(Xb5(i) + Zu(i, 4)));
    nll += -w(i) * loglik(fam, y(i), f, bn(i), phi);
    SIMULATE {
      y(i) = simobs(fam, f, bn(i), phi);
    }
  }
  SIMULATE {
    REPORT(y);
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
