//anaphase changepoint model as armond et al 2019
data {
  int<lower=1> T;
  int<lower=1> nTracks;
  matrix[T,2] y[nTracks]; //trajectory data from each sister
  vector[T] y_missing[nTracks]; //logical to indicate if data is missing
  real dt;
}
transformed data{
  vector[T] intersister_dist[nTracks];
  for (i in 1:nTracks){
    for (t in 1:T){
      intersister_dist[i][t] = fabs(y[i][t,1] - y[i][t,2]);
    }
  }
}
parameters {
  real<lower=0,upper=(T*dt)> t_ana[nTracks];
  real alpha[nTracks];
  real<lower=0> beta[nTracks];
  real<lower=0> tau[nTracks];
}
transformed parameters {
  real a[nTracks];
  for (i in 1:nTracks) {
    a[i] = alpha[i] + beta[i]*t_ana[i];
  }
}
model {
  // priors
  tau ~ gamma(2, 1.0/5);
  alpha ~ normal(0, 100);
  t_ana ~ uniform(0, T*dt);
  // change point model
  for (i in 1:nTracks) {
      beta[i] ~ normal(1.0/60, 100) T[0,]; // prior for beta, truncated distns must be univaraiate
    for (t in 1:T) {
      if (y_missing[i][t]<=0){ //check if data is missing
        intersister_dist[i][t] ~ normal((t*dt < t_ana[i]) ? a[i] : alpha[i] + beta[i]*t*dt,1/sqrt(tau[i]));
      }
    }
  }
}
generated quantities {
  vector[T] D_sim[nTracks];
  vector[nTracks] log_lik;
  matrix[nTracks,T] ll;
  for (i in 1:nTracks){
    for (t in 1:T) {
      if (y_missing[i][t]<=0){ //check if data is missing
        ll[i,t] = normal_lpdf(intersister_dist[i][t] | (t*dt < t_ana[i]) ? a[i] : alpha[i] + beta[i]*t*dt,1/sqrt(tau[i]));
      } else {
        ll[i,t] = 0;
      }
      D_sim[i][t] = normal_rng((t*dt < t_ana[i]) ? a[i] : alpha[i] + beta[i]*t*dt, 1/sqrt(tau[i]));
    }
    log_lik[i] = sum(ll[i]);  
  }
}
