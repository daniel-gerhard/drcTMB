#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type drcTMB(objective_function<Type>* obj) {
  DATA_INTEGER(fam);
  DATA_INTEGER(mod);
  DATA_IVECTOR(lnk);       // links of b1-b5: 1 = identity, 2 = log, 3 = logit
  DATA_SCALAR(alpha);      // exponent of the hormesis term of the cedergreen model
  DATA_SCALAR(smooth);     // smoothing width of the threshold of the nec model, 0 for an exact threshold
  DATA_IVECTOR(tnp);       // random effect terms: number of random coefficients per term
  DATA_IVECTOR(tnl);       // number of levels per term
  DATA_IVECTOR(tpar);      // parameter indices (0-4) of the random coefficients of all terms, concatenated
  DATA_IVECTOR(tzb);       // first column in Z (and Zp) of the design block of each random coefficient
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
  DATA_SPARSE_MATRIX(Z);   // design blocks of all terms, side by side: a column for each level, with the
                           // covariate of a random slope (1 for intercepts); zero columns without random effects
  DATA_VECTOR(xp);         // predictions (predict with se.fit): doses, length zero otherwise
  DATA_MATRIX(Xp1);
  DATA_MATRIX(Xp2);
  DATA_MATRIX(Xp3);
  DATA_MATRIX(Xp4);
  DATA_MATRIX(Xp5);
  DATA_SPARSE_MATRIX(Zp);  // random effects of the fitted levels; zero rows for population predictions
  DATA_IVECTOR(prwhat);    // priors: 0 = coefficient of b1-b5, 1 = theta, 2 = log_sigma; length zero without priors
  DATA_IVECTOR(pridx);     // index (from 0) in the coefficients b1-b5, concatenated, or in theta
  DATA_IVECTOR(prtrans);   // from the parameter to the scale of the prior: 1 identity, 2 exp, 3 inverse logit
  DATA_IVECTOR(prtype);    // distribution, see lpdfprior
  DATA_MATRIX(prhyp);      // hyperparameters (columns 0-3) and the log probability of the range of the parameter (column 4)
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(b3);
  PARAMETER_VECTOR(b4);
  PARAMETER_VECTOR(b5);
  PARAMETER_VECTOR(u);     // length zero without random effects
  PARAMETER_VECTOR(theta); // length zero without random effects
  PARAMETER(log_sigma);    // log sd (gaussian; of log(y) for lognormal), precision (beta, betabinomial), size (nbinom2) or dispersion (nbinom1);
                           // mapped (fixed) for poisson and binomial

  Type phi = exp(log_sigma);
  if (fam == 1 || fam == 8){
    Type sigma = phi;
    ADREPORT(sigma);
  }
  if (fam == 4){
    ADREPORT(phi);
  }

  Type nll = 0;
  int n = y.size();

  // priors: the negative log density of each prior at its parameter, on the scale of the prior, with the
  // log Jacobian of the transformation from the parameter, renormalised to the range of the parameter
  Type nlprior = 0;
  int npr = prwhat.size();
  if (npr > 0){
    vector<Type> ball(b1.size() + b2.size() + b3.size() + b4.size() + b5.size());
    ball << b1, b2, b3, b4, b5;
    for (int k = 0; k < npr; k++){
      Type v = prwhat(k) == 0 ? ball(pridx(k)) : (prwhat(k) == 1 ? theta(pridx(k)) : log_sigma);
      Type z = linkinv(prtrans(k), v);
      Type logjac = prtrans(k) == 2 ? v : (prtrans(k) == 3 ? Type(log(z) + log(1 - z)) : Type(0));
      vector<Type> h = prhyp.row(k);
      nlprior -= lpdfprior(prtype(k), z, h) + logjac - prhyp(k, 4);
    }
  }
  REPORT(nlprior);
  nll += nlprior;

  // random effects: column k of Zu is added to the linear predictor of b(k+1)
  // term t has tnl(t) levels and tnp(t) random coefficients (intercepts or slopes) of the parameters
  // tpar, with design blocks at columns tzb of Z, correlated within a term and independent between terms
  matrix<Type> Zu(n, 5);
  Zu.setZero();
  int nt = tnp.size();
  if (nt > 0){
    int ncorr = 0;
    for (int t = 0; t < nt; t++) ncorr += tnp(t)*tnp(t);
    vector<Type> corr(ncorr);

    int uoff = 0, toff = 0, poff = 0, coff = 0;
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
        vector<Type> Zuk = Z.middleCols(tzb(poff + k), nl) * U.col(k).matrix();
        Zu.col(tpar(poff + k)) += Zuk.matrix();
      }

      uoff += nl*np;
      toff += np + np*(np-1)/2;
      poff += np;
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
                linkinv(lnk(0), Type(Xb1(i) + Zu(i, 0))),
                linkinv(lnk(1), Type(Xb2(i) + Zu(i, 1))),
                linkinv(lnk(2), Type(Xb3(i) + Zu(i, 2))),
                linkinv(lnk(3), Type(Xb4(i) + Zu(i, 3))),
                linkinv(lnk(4), Type(Xb5(i) + Zu(i, 4))), alpha, smooth);
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
    int uoff = 0, poff = 0;
    for (int t = 0; t < nt; t++){
      int np = tnp(t);
      int nl = tnl(t);
      for (int k = 0; k < np; k++){
        vector<Type> uk = u.segment(uoff + k*nl, nl);
        vector<Type> Zuk = Zp.middleCols(tzb(poff + k), nl) * uk.matrix();
        Zup.col(tpar(poff + k)) += Zuk.matrix();
      }
      uoff += nl*np;
      poff += np;
    }
    vector<Type> Xpb1 = Xp1 * b1;
    vector<Type> Xpb2 = Xp2 * b2;
    vector<Type> Xpb3 = Xp3 * b3;
    vector<Type> Xpb4 = Xp4 * b4;
    vector<Type> Xpb5 = Xp5 * b5;
    vector<Type> mu(npred);
    for (int i = 0; i < npred; i++){
      mu(i) = drcmean(mod, xp(i),
                      linkinv(lnk(0), Type(Xpb1(i) + Zup(i, 0))),
                      linkinv(lnk(1), Type(Xpb2(i) + Zup(i, 1))),
                      linkinv(lnk(2), Type(Xpb3(i) + Zup(i, 2))),
                      linkinv(lnk(3), Type(Xpb4(i) + Zup(i, 3))),
                      linkinv(lnk(4), Type(Xpb5(i) + Zup(i, 4))), alpha, smooth);
    }
    ADREPORT(mu);
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
