#include "drcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type drcTMB(objective_function<Type>* obj) {
  DATA_INTEGER(fam);
  DATA_INTEGER(mod);
  DATA_IVECTOR(lnk);       // links of b1-b5: 1 = identity, 2 = log, 3 = logit
  DATA_INTEGER(clink);     // link of the curve: 1 = identity, 2 = log, 3 = logit, 4 = probit; the mean (median, probability) is linkinv(clink, curve)
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
  DATA_INTEGER(bmdpar);    // 1: BMD parameterisation, b4 holds log BMD of the curve of the intercepts (see bmdb4); 0 otherwise
  DATA_INTEGER(b4role);    // BMD parameterisation: role of b4, see bmdb4
  DATA_INTEGER(bdef);      // definition of the benchmark response, see bmdtarget
  DATA_SCALAR(bbmr);       // benchmark response
  DATA_INTEGER(binc);      // direction of the curve: 1 increasing, -1 decreasing
  DATA_INTEGER(b4pos);
  DATA_SCALAR(bdmax);      // BMD parameterisation of the quadratic exponential model (b4role 5): the highest dose     // 1: b4 must not be negative (default lower bound 0 of the standard parameterisation)
  DATA_IVECTOR(prwhat);    // priors: 0 = coefficient of b1-b5, 1 = theta, 2 = log_sigma, 3 = cpar, 4 = vpar; length zero without priors
  DATA_IVECTOR(pridx);     // index (from 0) in the coefficients b1-b5, concatenated, or in theta
  DATA_IVECTOR(prtrans);   // from the parameter to the scale of the prior: 1 identity, 2 exp, 3 inverse logit, 5 tanh, 6 compound symmetry
  DATA_IVECTOR(prtype);    // distribution, see lpdfprior
  DATA_IVECTOR(prpost);    // the prior applies to: 1 the parameter, 2 its negative (negate), 3 its exponential, 4 pnorm, 5 plogis
  DATA_MATRIX(prhyp);      // hyperparameters (columns 0-3) and the log probability of the range of the parameter (column 4)
  // residual correlation structure (gaussian and lognormal): 0 none, 1 AR1, 2 CAR1, 3 compound symmetry,
  // 4 exponential, 5 Gaussian, 6 spherical; series of the observations cord (indices from 0, by series)
  // with the lengths clen, positions cpos (a row per observation), largest integer lag cmaxk (AR1)
  DATA_INTEGER(ctype);
  DATA_INTEGER(cnug);      // 1: nugget (spatial structures)
  DATA_IVECTOR(cord);
  DATA_IVECTOR(clen);
  DATA_MATRIX(cpos);
  DATA_INTEGER(cmaxk);
  // variance function (gaussian and lognormal): 0 none, 1 power, 2 exponential, 3 separate standard
  // deviations of the strata vstr; covariate vcov, or the fitted mean for vfit = 1
  DATA_INTEGER(vtype);
  DATA_INTEGER(vfit);
  DATA_VECTOR(vcov);
  DATA_IVECTOR(vstr);
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(b3);
  PARAMETER_VECTOR(b4);
  PARAMETER_VECTOR(b5);
  PARAMETER_VECTOR(u);     // length zero without random effects
  PARAMETER_VECTOR(theta); // length zero without random effects
  PARAMETER(log_sigma);    // log sd (gaussian; of log(y) for lognormal), precision (beta, betabinomial), size (nbinom2) or dispersion (nbinom1);
                           // mapped (fixed) for poisson and binomial
  PARAMETER_VECTOR(cpar);  // correlation parameters, unconstrained (see corvalue); length zero without
  PARAMETER_VECTOR(vpar);  // variance function parameters: delta (power, exponential) or the log ratios of
                           // the standard deviations of the strata to the first (separate); length zero without

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
      Type v = prwhat(k) == 0 ? ball(pridx(k)) : (prwhat(k) == 1 ? theta(pridx(k)) : (prwhat(k) == 2 ? log_sigma :
               (prwhat(k) == 3 ? cpar(pridx(k)) : vpar(pridx(k)))));
      Type z = linkinv(prtrans(k), v);
      Type logjac = prtrans(k) == 2 ? v : (prtrans(k) == 3 ? Type(log(z) + log(1 - z)) : Type(0));
      // correlation of AR1 (tanh) and of compound symmetry (in (lo, 1) for series of up to m observations)
      if (prtrans(k) == 5){
        z = tanh(v);
        logjac = log(1 - z * z);
      }
      if (prtrans(k) == 6){
        Type lo = clen.maxCoeff() > 1 ? Type(-1.0 / (clen.maxCoeff() - 1)) : Type(-1);
        Type p = invlogit(v);
        z = lo + (1 - lo) * p;
        logjac = log(1 - lo) + log(p) + log(1 - p);
      }
      // a prior of the negative of the parameter, or of its exponential, pnorm or plogis (with the log
      // Jacobian of the inverse link)
      if (prpost(k) == 2) z = -z;
      if (prpost(k) == 3){
        logjac += z;
        z = exp(z);
      }
      if (prpost(k) == 4){
        logjac += dnorm(z, Type(0), Type(1), true);
        z = pnorm(z);
      }
      if (prpost(k) == 5){
        logjac += -logspace_add(Type(0), -z) - logspace_add(Type(0), z);
        z = invlogit(z);
      }
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

  // BMD parameterisation: b4 from the BMD and the other parameters of the curve of the intercepts
  // (without covariates; random effects of the groups on b1-b3 and b5 are added to that curve)
  Type f4bmd = 0;
  Type f3bmd = 0;
  bool qe = bmdpar == 1 && b4role == 5;
  if (qe){
    // quadratic exponential: b3 holds log d, and both coefficients of the polynomial follow from the BMD
    vector<Type> qc = bmdqe(Type(exp(b4(0))), linkinv(lnk(0), b1(0)), linkinv(lnk(1), b2(0)), b3(0), bdef, bbmr, binc, phi, clink, bdmax);
    f3bmd = qc(0);
    f4bmd = qc(1);
    REPORT(f3bmd);
    REPORT(f4bmd);
  } else if (bmdpar == 1){
    f4bmd = bmdb4(mod, b4role, Type(exp(b4(0))), linkinv(lnk(0), b1(0)), linkinv(lnk(1), b2(0)),
                  linkinv(lnk(2), b3(0)), linkinv(lnk(4), b5(0)), bdef, bbmr, binc, phi, smooth, clink);
    REPORT(f4bmd);
  }

  Type f;
  // gaussian and lognormal residuals with a correlation structure or a variance function
  bool resid = (fam == 1 || fam == 8) && (ctype > 0 || vtype > 0);
  vector<Type> fy(n);
  for (int i = 0; i < n; i++){
    f = drcmean(mod, x(i),
                linkinv(lnk(0), Type(Xb1(i) + Zu(i, 0))),
                linkinv(lnk(1), Type(Xb2(i) + Zu(i, 1))),
                qe ? f3bmd : linkinv(lnk(2), Type(Xb3(i) + Zu(i, 2))),
                bmdpar == 1 ? f4bmd : linkinv(lnk(3), Type(Xb4(i) + Zu(i, 3))),
                linkinv(lnk(4), Type(Xb5(i) + Zu(i, 4))), alpha, smooth);
    f = linkinv(clink, f);
    fy(i) = f;
    if (!resid){
      nll +=-w(i) * loglik(fam, y(i), f, bn(i), phi);
      SIMULATE {
        y(i) = simobs(fam, f, bn(i), phi);
      }
    }
  }
  if (resid){
    // residuals (of log(y) for lognormal) and their standard deviations sigma g
    vector<Type> e(n), sd(n);
    for (int i = 0; i < n; i++){
      e(i) = fam == 8 ? Type(log(y(i)) - log(fy(i))) : Type(y(i) - fy(i));
      sd(i) = vtype == 0 ? phi : Type(phi * varfun(vtype, vfit == 1 ? fy(i) : Type(vcov(i)), vstr(i), vpar));
    }
    if (ctype == 0){
      // independent residuals, weighted as the other families
      for (int i = 0; i < n; i++){
        nll -= w(i) * (dnorm(e(i), Type(0), sd(i), true) - (fam == 8 ? Type(log(y(i))) : Type(0)));
        SIMULATE {
          e(i) = rnorm(Type(0), sd(i));
        }
      }
    } else {
      // multivariate normal residuals of each series, with covariance D R D
      vector<Type> rpow(cmaxk + 1);
      rpow(0) = 1;
      if (ctype == 1) for (int k = 1; k <= cmaxk; k++) rpow(k) = rpow(k - 1) * tanh(cpar(0));
      int off = 0;
      for (int s = 0; s < clen.size(); s++){
        int m = clen(s);
        matrix<Type> S(m, m);
        vector<Type> es(m);
        for (int a = 0; a < m; a++){
          int ia = cord(off + a);
          es(a) = e(ia);
          for (int b = 0; b <= a; b++){
            int ib = cord(off + b);
            Type r = a == b ? Type(1) : corvalue(ctype, cnug, cpar, rpow, vector<Type>(cpos.row(ia)), vector<Type>(cpos.row(ib)), clen.maxCoeff());
            S(a, b) = S(b, a) = r * sd(ia) * sd(ib);
          }
        }
        density::MVNORM_t<Type> nld(S);
        nll += nld(es);
        SIMULATE {
          es = nld.simulate();
          for (int a = 0; a < m; a++) e(cord(off + a)) = es(a);
        }
        off += m;
      }
      if (fam == 8) for (int i = 0; i < n; i++) nll += log(y(i));
    }
    SIMULATE {
      for (int i = 0; i < n; i++) y(i) = fam == 8 ? Type(fy(i) * exp(e(i))) : Type(fy(i) + e(i));
    }
  }
  SIMULATE {
    REPORT(y);
  }
  // BMD parameterisation: a computed b4 below its bound 0 is outside the parameter space
  if (bmdpar == 1 && b4pos == 1){
    nll = CppAD::CondExpLt(f4bmd, Type(0), Type(std::numeric_limits<double>::infinity()), nll);
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
                      qe ? f3bmd : linkinv(lnk(2), Type(Xpb3(i) + Zup(i, 2))),
                      bmdpar == 1 ? f4bmd : linkinv(lnk(3), Type(Xpb4(i) + Zup(i, 3))),
                      linkinv(lnk(4), Type(Xpb5(i) + Zup(i, 4))), alpha, smooth);
      mu(i) = linkinv(clink, mu(i));
    }
    ADREPORT(mu);
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
