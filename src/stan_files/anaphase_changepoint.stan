//anaphase changepoint model as armond et al 2019
data {
  int<lower=1> T;
  matrix[T,2] y; //trajectory data from each sister
  vector[T] y_missing; //logical to indicate if data is missing
  real dt;
}
transformed data{
  vector[T] intersister_dist;
  for (t in 1:T){
    intersister_dist[t] = fabs(y[t,1] - y[t,2]);
  }
}
parameters {
  real<lower=0,upper=(T*dt)> t_ana;
  real alpha;
  real<lower=0> beta;
  real<lower=0> tau;
}
transformed parameters {
  real a;
  a = alpha + beta*t_ana;
}
model {
  // priors
  tau ~ gamma(2, 1.0/5);
  alpha ~ normal(0, 100);
  beta ~ normal(1.0/60, 100) T[0,];
  t_ana ~ uniform(0, T*dt);
  // change point model
  for (t in 1:T) {
     if (y_missing[t]<=0){ //check if data is missing
       intersister_dist[t] ~ normal((t*dt < t_ana) ? a : alpha + beta*t*dt,1/sqrt(tau));
     }
  }
}
generated quantities {
  vector[T] D_sim;
  for (t in 1:T) {
    D_sim[t] = normal_rng((t*dt < t_ana) ? a : alpha + beta*t*dt, 1/sqrt(tau));
  }  
}
