//The idea of this model variation is that we can in some sense "subtract"
// the PEF and spring forces such that the remaining force is accounted for by
// microtubule forces which are inferred

data {
  int<lower=1> T;
  int<lower=0> T0;
  int<upper=T> T1;
  int y_missing[T]; //is data misssing at given time point
  matrix[T,2] y; //trajectory data from each sister
  real dt;
  vector[T] cos_phi; //the angle to metaphase plate is phi
}
parameters {
  real<lower=0> alpha;
  real<lower=0> kappa;
  real<lower=0> L;
  real<lower=0> tau;
  matrix[T,2] mt_force;
  real<lower=0> phi;
  real<lower=0> epsilon;
}
model {
  row_vector[2] pred;
  //priors ...
  alpha ~ normal(0.01,0.1) T[0,];
  kappa ~ normal(0.05,0.1) T[0,];
  tau ~ gamma(0.5,1.0/1000);
  L ~ normal(0.790,0.119) T[0,];
  phi ~ gamma(100,100);
  epsilon ~ normal(0,0.1) T[0,];
  mt_force[1] ~ cauchy(0,epsilon);
  for (t in 2:T){
    mt_force[t] ~ cauchy(mt_force[t-1]*phi,epsilon);
  }
  for (t in (T0+2):T1) {
    if (y_missing[t] == 0){
      pred = y[t-1] + dt*[-mt_force[t,1] -kappa*(y[t-1,1]-y[t-1,2]-L*cos_phi[t]) - alpha*y[t-1,1],mt_force[t,2] -kappa*(y[t-1,2]-y[t-1,1]+L*cos_phi[t]) - alpha*y[t-1,2]];
      target += normal_lpdf(y[t] | pred, sqrt(dt/tau));
    }
  }
}
