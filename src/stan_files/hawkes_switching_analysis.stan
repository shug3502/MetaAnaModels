functions {
  real gaussian(real x, real lengthscale) {
    return(1/(2*pi())*exp(-.5 * x^2/lengthscale^2)/lengthscale^2);
  }
}
data {
  real prior_space[2];
  real prior_time[2];

  int<lower=1> n;
  int<lower=1> nPairs;
  matrix[n,2] Space;
  vector[n] Time; // NB: these must be sorted from smallest to largest!
  int SisterPair[n];
  real time_window;
  real space_window;
}
transformed data {
  matrix[n,n] timeD;
  matrix[n,n] spaceD;
  vector[n] timeD2;
  matrix[nPairs,nPairs] SisterPairD = rep_matrix(1,nPairs,nPairs) - diag_matrix(rep_vector(1,nPairs));
  for(i in 1:n) {
    for(j in 1:n) {
      timeD[i,j] = -(Time[i] - Time[j]);
      spaceD[i,j] = distance(Space[i], Space[j]);
    }
    timeD2[i] = -(time_window - Time[i]); 
  }
print(SisterPairD);
}
parameters {
  real<lower=0> lengthscaleS;
  real<lower=0> lengthscaleT;
  real<lower=0> a;
  real<lower=0> mu0;
}
transformed parameters {
  vector[n] ll = rep_vector(mu0,n);
  real lp;
  for(i in 2:n) { // todo: vectorize these calculations
    for(j in 1:(i-1)) {
      ll[i] = ll[i] + (a * lengthscaleT * exp(timeD[i,j] * lengthscaleT) * gaussian(spaceD[i,j],lengthscaleS)) * SisterPairD[SisterPair[i],SisterPair[j]];
    }
  }
  lp = sum(log(ll)) - mu0 * space_window * time_window +
          a * (sum(exp(timeD2*lengthscaleT)) - n);
}

model {
  lengthscaleS ~ normal(prior_space[1],prior_space[2]);
  lengthscaleT ~ normal(prior_time[1],prior_time[2]);
  //lengthscaleS ~ lognormal(3,1);
  //lengthscaleT ~ lognormal(3,1);
  a ~  normal(0,10);
  mu0 ~ normal(0,10);

  target += lp;    
}

generated quantities {
  vector[n] background;
  real lengthscale;
  lengthscale = 1/lengthscaleT;
  background = (mu0) ./ ll;
}


