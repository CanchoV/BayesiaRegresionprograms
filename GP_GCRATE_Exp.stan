
//

functions
{
    real lambertW(real x)
    {
      real y;
      real w;
      y= sqrt( 1 + exp(1) *x);
      w=-1 + 2.036 * log( (1 + 1.14956131 * y)/(1 + 0.45495740*log(1+y)) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      return(w);
    }
}

// The input data is a vector 'y' of length 'n'.
data {
 int<lower=1> n; //number of observation
 int<lower=1> p; //number of regressors
   vector[n] status; //status
   matrix[n, p] X; //predictor
    vector[n] y; // response
}

// The parameters accepted by the model. Our model

parameters {
 vector[p] beta; // regression coefficients
 real<lower=0.000001,upper=0.99990> gama; //parameter of the GP  distribution
 real lambda1; //scale parameter of the Weibull distribution
 //real<lower=0.000001> lambda2; //shape parameter of the Weibull  distribution
}



// The model to be estimated. We model the output

model {
  vector[n] eta = X*beta; // linear predictor
  real mu[n];
  real lambda[n];
  real theta[n]; // Media
  real vf[n];
  real vS[n];
  real ax[n];
  real Spop[n];
  real fpop[n];
 
 for (i in 1:n)
  {
  mu[i]=exp(eta[i]);
  lambda[i]=(1-gama)*mu[i];
  vf[i]=exp(lambda1)*exp(-exp(lambda1)*y[i]);
  vS[i]=exp(-exp(lambda1)*y[i]) ;
  ax[i]=lambertW(-gama*(vS[i])*exp(-gama));
  Spop[i]=exp(-(lambda[i]/gama)*(ax[i]+gama)) ;    
  fpop[i]=-(lambda[i]/gama)*(vf[i]/(vS[i]))*Spop[i]*(ax[i]/(1+ax[i])) ;
  target +=status[i]*log(fpop[i])+(1-status[i])*log(Spop[i]);
   }
  
  gama~beta(2,2);
  lambda1~normal(0, 1);
  for(i in 1:p){
    beta[i]~normal(0, 10); // priori 
  }
}
