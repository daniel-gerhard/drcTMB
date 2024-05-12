#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type drcbetafix(objective_function<Type>* obj) {
  DATA_INTEGER(mod); 
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  DATA_MATRIX(X3);
  DATA_MATRIX(X4);
  DATA_MATRIX(X5);
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(b3);
  PARAMETER_VECTOR(b4);
  PARAMETER_VECTOR(b5); 
  PARAMETER(log_sigma);
  
  Type phi = exp(log_sigma);
  ADREPORT(phi);
  
  Type nll;
  
  vector<Type> Xb1 = X1 * b1;
  vector<Type> Xb2 = X2 * b2;
  vector<Type> Xb3 = X3 * b3;
  vector<Type> Xb4 = X4 * b4;
  vector<Type> Xb5 = X5 * b5;
  
  Type fl;
  Type f;
  Type logl;

  int i;
  for (i=0; i < y.size(); i++){
    switch (mod){
    case 1:
      fl = exp(Xb3(i)*(x(i) - Xb4(i)));
      f = Xb2(i) + (Xb1(i) - Xb2(i)) / (1 + pow(fl, Xb5(i)));
      break;
    case 2:
      fl = exp(Xb3(i)*(log(x(i)) - log(Xb4(i))));
      f = Xb2(i) + (Xb1(i) - Xb2(i)) / (1 + pow(fl, Xb5(i)));
      break;
    case 3:
      fl = exp(-1*exp(Xb3(i)*(log(x(i)) - log(Xb4(i)))));
      f = Xb2(i) + (Xb1(i) - Xb2(i)) * fl;
      break;
    case 4:
      fl = 1 - exp(-1*exp(Xb3(i)*(log(x(i)) - log(Xb4(i)))));
      f = Xb2(i) + (Xb1(i) - Xb2(i)) * fl;
      break;
    case 5:
      fl = Xb3(i)*(log(x(i)) - log(Xb4(i)));
      f = Xb2(i) + (Xb1(i) - Xb2(i)) * pnorm(fl);
      break;
    }
    logl = lgamma(phi) - lgamma(f*phi) - lgamma((1-f)*phi) + (f*phi-1)*log(y(i)) + ((1-f)*phi - 1)*log(1 - y(i));
    nll += -1*logl;
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this