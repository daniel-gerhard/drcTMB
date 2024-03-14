#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type drcnorm(objective_function<Type>* obj) {
  DATA_INTEGER(mod); 
  DATA_IVECTOR(ind);
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_VECTOR(z0);
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  DATA_MATRIX(X3);
  DATA_MATRIX(X4);
  DATA_MATRIX(X5);
  DATA_SPARSE_MATRIX(Z);
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(b3);
  PARAMETER_VECTOR(b4);
  PARAMETER_VECTOR(b5); 
  PARAMETER_VECTOR(u);
  PARAMETER_VECTOR(theta);  
  PARAMETER(log_sigma);

  Type sigma = exp(log_sigma);
  ADREPORT(sigma);

  // dimensions of u
  int rn = ind.sum();
  int un = u.size()/rn;
  vector<int> dim(2);
  dim << un, rn;
  array<Type> U(u, dim);
  
  Type nll;
  
  vector<Type> logsd = theta.head(rn);
  vector<Type> corr_transf = theta.tail(theta.size() - rn);
  vector<Type> sd = exp(logsd);
  density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
  density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
  
  for(int i = 0; i < un; i++){
    nll += scnldens(U.transpose().col(i));
  }
  
  matrix<Type> corr = nldens.cov(); 
  ADREPORT(corr);
  
  vector<Type> Zu1 = z0;
  vector<Type> Zu2 = z0;
  vector<Type> Zu3 = z0;
  vector<Type> Zu4 = z0;
  vector<Type> Zu5 = z0;
  int ucount = 0;
  if (ind(0) > 0){
    Zu1 = Z * U.col(ucount).matrix();
    ucount = ucount + 1;
  }
  if (ind(1) > 0){
    Zu2 = Z * U.col(ucount).matrix();
    ucount = ucount + 1;   
  }
  if (ind(2) > 0){
    Zu3 = Z * U.col(ucount).matrix();
    ucount = ucount + 1; 
  }
  if (ind(3) > 0){
    Zu4 = Z * U.col(ucount).matrix();
    ucount = ucount + 1; 
  }  
  if (ind(4) > 0){
    Zu5 = Z * U.col(ucount).matrix();
  }
  vector<Type> Xb1 = X1 * b1;
  vector<Type> Xb2 = X2 * b2;
  vector<Type> Xb3 = X3 * b3;
  vector<Type> Xb4 = X4 * b4;
  vector<Type> Xb5 = X5 * b5;
  
  Type f1;
  Type f2;
  Type f3;
  Type f4;
  Type f5; 
  Type fl;
  Type f;
  
  int i;
  for (i=0; i < y.size(); i++){
    f1 = Xb1(i) + Zu1(i);
    f2 = Xb2(i) + Zu2(i);
    f3 = Xb3(i) + Zu3(i);
    f4 = Xb4(i) + Zu4(i);
    f5 = Xb5(i) + Zu5(i);
    switch (mod){
    case 1:
      fl = exp(f3*(x(i) - f4));
      f = f2 + (f1 - f2) / (1 + pow(fl, f5));
      break;
    case 2:
      fl = exp(f3*(log(x(i)) - log(f4)));
      f = f2 + (f1 - f2) / (1 + pow(fl, f5));
      break;
    case 3:
      fl = exp(-1*exp(f3*(log(x(i)) - log(f4))));
      f = f2 + (f1 - f2) * fl;
      break;
    case 4:
      fl = 1 - exp(-1*exp(f3*(log(x(i)) - log(f4))));
      f = f2 + (f1 - f2) * fl;
      break;
    case 5:
      fl = f3*(log(x(i)) - log(f4));
      f = f2 + (f1 - f2) * pnorm(fl);
      break;
    }
    nll += -dnorm(y(i), f, exp(log_sigma), true);
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
