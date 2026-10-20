#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type drcTMB(objective_function<Type>* obj) {
  DATA_INTEGER(fam);
  DATA_INTEGER(mod);
  DATA_INTEGER(lnk);
  DATA_IVECTOR(tnp);       // random effect terms: number of parameters per term
  DATA_IVECTOR(tnl);       // number of levels per term
  DATA_IVECTOR(tpar);      // parameter indices (0-4) of all terms, concatenated
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
  DATA_SPARSE_MATRIX(Z);   // design matrices of all terms, side by side; zero columns without random effects
  DATA_VECTOR(xp);         // predictions (predict with se.fit): doses, length zero otherwise
  DATA_MATRIX(Xp1);
  DATA_MATRIX(Xp2);
  DATA_MATRIX(Xp3);
  DATA_MATRIX(Xp4);
  DATA_MATRIX(Xp5);
  DATA_SPARSE_MATRIX(Zp);  // random effects of the fitted levels; zero rows for population predictions
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
  // term t has tnl(t) levels and random intercepts for the tnp(t) parameters tpar,
  // correlated within a term and independent between terms
  matrix<Type> Zu(n, 5);
  Zu.setZero();
  int nt = tnp.size();
  if (nt > 0){
    int ncorr = 0;
    for (int t = 0; t < nt; t++) ncorr += tnp(t)*tnp(t);
    vector<Type> corr(ncorr);

    int uoff = 0, toff = 0, poff = 0, zoff = 0, coff = 0;
    for (int t = 0; t < nt; t++){
      int np = tnp(t);
      int nl = tnl(t);
      vector<Type> ut = u.segment(uoff, nl*np);
      vector<int> dim(2);
      dim << nl, np;
      array<Type> U(ut, dim);

      vector<Type> logsd = theta.segment(toff, np);
      vector<Type> corr_transf = theta.segment(toff + np, np*(np-1)/2);
      vector<Type> sd = exp(logsd);
      density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
      density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);

      for (int j = 0; j < nl; j++){
        nll += scnldens(U.transpose().col(j));
      }

      // new random effects for each level
      SIMULATE {
        if (simre == 1){
          vector<Type> uj(np);
          for (int j = 0; j < nl; j++){
            scnldens.simulate(uj);
            for (int k = 0; k < np; k++) U(j, k) = uj(k);
          }
        }
      }

      matrix<Type> ct = nldens.cov();
      for (int l = 0; l < np; l++){
        for (int k = 0; k < np; k++) corr(coff + k + l*np) = ct(k, l);
      }

      for (int k = 0; k < np; k++){
        vector<Type> Zuk = Z.middleCols(zoff, nl) * U.col(k).matrix();
        Zu.col(tpar(poff + k)) += Zuk.matrix();
      }

      uoff += nl*np;
      toff += np + np*(np-1)/2;
      poff += np;
      zoff += nl;
      coff += np*np;
    }
    ADREPORT(corr);
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

  // mean at the doses xp, for standard errors of the predictions by the delta method
  int npred = xp.size();
  if (npred > 0){
    matrix<Type> Zup(npred, 5);
    Zup.setZero();
    int uoff = 0, poff = 0, zoff = 0;
    for (int t = 0; t < nt; t++){
      int np = tnp(t);
      int nl = tnl(t);
      for (int k = 0; k < np; k++){
        vector<Type> uk = u.segment(uoff + k*nl, nl);
        vector<Type> Zuk = Zp.middleCols(zoff, nl) * uk.matrix();
        Zup.col(tpar(poff + k)) += Zuk.matrix();
      }
      uoff += nl*np;
      poff += np;
      zoff += nl;
    }
    vector<Type> Xpb1 = Xp1 * b1;
    vector<Type> Xpb2 = Xp2 * b2;
    vector<Type> Xpb3 = Xp3 * b3;
    vector<Type> Xpb4 = Xp4 * b4;
    vector<Type> Xpb5 = Xp5 * b5;
    vector<Type> mu(npred);
    for (int i = 0; i < npred; i++){
      mu(i) = drcmean(mod, xp(i),
                      linkinv(lnk, Type(Xpb1(i) + Zup(i, 0))),
                      linkinv(lnk, Type(Xpb2(i) + Zup(i, 1))),
                      Type(Xpb3(i) + Zup(i, 2)),
                      Type(Xpb4(i) + Zup(i, 3)),
                      Type(Xpb5(i) + Zup(i, 4)));
    }
    ADREPORT(mu);
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
