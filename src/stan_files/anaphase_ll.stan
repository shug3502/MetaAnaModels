// see https://khakieconomics.github.io/2018/02/24/Regime-switching-models.html
//metaphase model. based on armond et al 2015 plos comp biol

//This version includes updated prior for L based on nocodazole data
// Also includes angle to the metaphase plate (input as data measured via KiT) 
//This adjusts the spring based forces which may not act perpendicular to metaphase plate

functions {
  matrix constructTransitionMatrix(real p_icoh, real p_coh, real p_ana){
    matrix[5,5] P;
    real q_ana;
    real q_coh;
    real q_icoh;
    q_icoh = 1-p_icoh;
    q_coh = 1-p_coh;
    q_ana = 1-p_ana;
    P[1,1] = p_icoh*p_icoh*q_ana;
    P[2,1] = p_coh*q_coh*q_ana;
    P[3,1] = p_coh*q_coh*q_ana;
    P[4,1] = q_icoh*q_icoh*q_ana;
    P[5,1] = 0;
    P[1,2] = p_icoh*q_icoh*q_ana;
    P[2,2] = p_coh*p_coh*q_ana;
    P[3,2] = q_coh*q_coh*q_ana;
    P[4,2] = p_icoh*q_icoh*q_ana;
    P[5,2] = 0;
    P[1,3] = p_icoh*q_icoh*q_ana;
    P[2,3] = q_coh*q_coh*q_ana;
    P[3,3] = p_coh*p_coh*q_ana;
    P[4,3] = p_icoh*q_icoh*q_ana;
    P[5,3] = 0;
    P[1,4] = q_icoh*q_icoh*q_ana;
    P[2,4] = p_coh*q_coh*q_ana;
    P[3,4] = p_coh*q_coh*q_ana;
    P[4,4] = p_icoh*p_icoh*q_ana;
    P[5,4] = 0;
    P[1,5] = p_ana;
    P[2,5] = p_ana;
    P[3,5] = p_ana;
    P[4,5] = p_ana;
    P[5,5] = 1;
    return P;
  }
  matrix odeUpdateMatrix(real[] th) {
    matrix[2,7] M;
    M[1,1] = -th[3] - th[2];
    M[1,2] = th[3];
    M[1,3] = -th[5];
    M[1,4] = -th[5];
    M[1,5] = -th[4];
    M[1,6] = -th[4];
    M[1,7] = 0;
    M[2,1] = th[3];
    M[2,2] = -th[3] - th[2];
    M[2,3] = th[5];
    M[2,4] = th[4];
    M[2,5] = th[5];
    M[2,6] = th[4];
    M[2,7] = 0;
    return M;
  }
  vector construct_transition_column(real p_icoh, real p_coh,
                                     real p_ana, //real p_reversal, real p_revtoana,
                                     int sigma){
    vector[5] P_col;
    real q_icoh = 1-p_icoh;
    real q_coh = 1-p_coh;
    real q_ana = 1-p_ana;
    //real q_reversal = 1-p_reversal;
    //real q_revtoana = 1-p_revtoana;
    if (sigma==1){
    P_col = [p_icoh*p_icoh*q_ana,
             p_coh*q_coh*q_ana,
             p_coh*q_coh*q_ana,
             q_icoh*q_icoh*q_ana,
             0]';
          //0];
    } else if(sigma==2){
    P_col = [p_icoh*q_icoh*q_ana,
              p_coh*p_coh*q_ana,
              q_coh*q_coh*q_ana,
              p_icoh*q_icoh*q_ana,
              0]';
              //0];
    } else if (sigma==3){
    P_col = [p_icoh*q_icoh*q_ana,
              q_coh*q_coh*q_ana,
              p_coh*p_coh*q_ana,
              p_icoh*q_icoh*q_ana,
              0]';
              //0];
    } else if (sigma==4){
    P_col = [q_icoh*q_icoh*q_ana,
              p_coh*q_coh*q_ana,
              p_coh*q_coh*q_ana,
              p_icoh*p_icoh*q_ana,
              0]';
              //0];
    } else if (sigma==5){
    P_col = [p_ana,
              p_ana,
              p_ana,
              p_ana,
              1]';
              //p_revtoana]';
    /*} else if (sigma==6){
    P_col = [0,
              0,
              0,
              0,
              p_reversal,
              q_revtoana]';
              */
    } else {
      P_col=rep_vector(0,5);
    }
    return(P_col);
  }
  vector odeUpdateVector(real[] th, real angleTheta){
    vector[2] mu;
    mu[1] = th[3]*th[8]*angleTheta;
    mu[2] = -th[3]*th[8]*angleTheta;
    return mu;
  }
  real myDiff(vector v){
    real w;
    w = fabs(v[1] - v[2]);
    return w;
  }
}
data {
  int<lower=1> T;
  int<lower=0> T0;
  int<upper=T> T1;
  int nStates;
  matrix[T,2] y; //trajectory data from each sister
  real dt;
  vector[nStates] sigma0;
  real t_ana_input;
  vector[T] cos_phi; //the angle to metaphase plate is phi
}
transformed data {
    vector[2] x0;
    real b; //sharpness of sigmoid switch to anaphase. Fixed not inferred. 
    x0 = y[T0+1, ]';
    b = dt/2;
}
parameters {
  real<lower=0> tau;
  real<lower=0> alpha;
  real<lower=0> kappa;
  real<upper=0> v_minus;
  real<lower=0> v_plus;
  real<lower=0,upper=1> p_icoh;
  real<lower=0,upper=1> p_coh;
  real<lower=0,upper=T*dt> t_ana; //estimate of time of anaphase
  real<lower=0> v_ana;
  real<lower=0> L;
}
transformed parameters{
  real theta[10]; //for metaphase ODE model
  matrix[T, nStates] eta;
  matrix[T, nStates] xi;
  vector[T] f;
  matrix[nStates,nStates] P;
  vector[nStates+2] auxStates;
  vector[T] p_ana;
  
  theta[1]=tau;
  theta[2]=alpha;
  theta[3]=kappa;
  theta[4]=v_minus;
  theta[5]=v_plus;
  theta[6]=p_icoh;
  theta[7]=p_coh;
  theta[8]=L;
  theta[9]=v_ana;
  theta[10]=t_ana;
for (t in 1:T0) {
  for (j in 1:5) {
    eta[t,j] = 0;
    xi[t,j] = 0;
  }
  f[t] = 1;
  p_ana[t] = 1/(1+exp(-(t*dt-t_ana)/b));
}
for (t in (T1+1):T) {
  for (j in 1:5) {
    eta[t,j] = 0;
    xi[t,j] = 0;
  }
  f[t] = 1;
  p_ana[t] = 1/(1+exp(-(t*dt-t_ana)/b));
}  
  // fill in etas
  for (t in (T0+1):T1) {
    if(t==(T0+1)) {
      for (j in 1:nStates){
        eta[t,j] = exp(normal_lpdf(y[t]| x0, sqrt(dt/tau)));
      }
    } else {
    for (j in 1:4){
            //construct vector from appending y[t-1] and hidden state eg [0 1 0 0]
      for (i in 1:(nStates+2)){
        if (i<=2) {
          auxStates[i] = y[t-1,i];
        } else if (i==(j+2)) {
          auxStates[i] = 1.0;
        } else {
          auxStates[i] = 0.0;
        }
      }
        eta[t,j] = exp(normal_lpdf(y[t]' | y[t-1]' + dt*odeUpdateMatrix(theta)*auxStates + dt*odeUpdateVector(theta, cos_phi[t]), sqrt(dt/tau)));
      }
      eta[t,5] = exp(normal_lpdf(y[t] | y[t-1] + dt*[v_ana, -v_ana], sqrt(dt/tau)));
    }
  }

  // work out likelihood contributions
  for(t in (T0+1):T1) {
    p_ana[t] = 1/(1+exp(-(t*dt-t_ana)/b));
    P = constructTransitionMatrix(p_icoh,p_coh,p_ana[t]);
    // for the first observation
    if(t==(T0+1)) {
      //replacing xi[t-1] by sigma0
      f[t] = dot_product((sigma0')*P, eta[t]');
      xi[t] = (((sigma0')*P) .* eta[t]) ./ f[t];
    } else {
    // and for the rest
      f[t] = dot_product(xi[t-1]*P, eta[t]'); 
      xi[t] = ((xi[t-1]*P) .* eta[t]) ./ f[t];
    }
  }
}
model {
  //priors ...
  tau ~ gamma(0.5,1.0/1000);
  alpha ~ normal(0.01,0.1) T[0,];
  kappa ~ normal(0.05,0.1) T[0,];
  v_minus ~ normal(-0.03,0.1) T[,0];
  v_plus ~ normal(0.03,0.1) T[0,];
  p_icoh ~ beta(12,3); //beta(2,1);
  p_coh ~ beta(45,5); //beta(2.5,1);
  v_ana ~ normal(0.03,0.1) T[0,];
  L ~ normal(0.790,0.119) T[0,];
  t_ana ~ normal(t_ana_input,14) T[0,T*dt]; 
  
  // likelihood is really easy here!
  target += sum(log(f));
}
generated quantities {
//sampling from hidden states could go here
  int sigma_sim[T];
  vector[nStates] state_probs = xi[T1]';
  vector[nStates] conditional_state_probs; //on log scale
  vector[nStates] P_col;
  int frame;
  vector[T] log_lik;

  for (t in (T1+1):T) {
    sigma_sim[t] = 0; //for time points with missing data
  }
  for (t in 1:T0) {
    sigma_sim[t] = 0; //for time points with missing data
  }

  //start with the final time pt
  sigma_sim[T1] = categorical_rng(state_probs);
  for (t in 1:(T1-T0-1)){
    frame = T1+1-t; //from T1 to T1+1-(T1-T0-1)=T0+2
    P_col = construct_transition_column(p_icoh,p_coh,p_ana[frame],sigma_sim[frame]);
    conditional_state_probs = log(P_col) + log(xi[frame-1]');
  if (is_inf(sum(conditional_state_probs))){
    sigma_sim[frame-1] = categorical_rng(exp(conditional_state_probs)/sum(exp(conditional_state_probs)));
  } else {
    sigma_sim[frame-1] = categorical_logit_rng(conditional_state_probs);
  }
  }

    log_lik = rep_vector(0,T);
    for (t in (T0+1):T1) {
      log_lik[t] = log(f[t]);
    }
}
