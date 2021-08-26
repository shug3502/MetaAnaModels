// see https://khakieconomics.github.io/2018/02/24/Regime-switching-models.html
//metaphase model. based on armond et al 2015 plos comp biol
//fit multiple trajectories. Assume all the same length.
functions {
  matrix constructTransitionMatrix(real p_icoh, real p_coh, real p_ana, real p_reversal, real p_revtoana){
    matrix[6,6] P;
    real q_ana;
    real q_coh;
    real q_icoh;
    real q_reversal;
    real q_revtoana;
    q_icoh = 1-p_icoh;
    q_coh = 1-p_coh;
    q_ana = 1-p_ana;
    q_reversal = 1-p_reversal;
    q_revtoana = 1-p_revtoana;
    P[1,1] = p_icoh*p_icoh*q_ana;
    P[2,1] = p_coh*q_coh*q_ana;
    P[3,1] = p_coh*q_coh*q_ana;
    P[4,1] = q_icoh*q_icoh*q_ana;
    P[5,1] = 0;
    P[6,1] = 0;
    P[1,2] = p_icoh*q_icoh*q_ana;
    P[2,2] = p_coh*p_coh*q_ana;
    P[3,2] = q_coh*q_coh*q_ana;
    P[4,2] = p_icoh*q_icoh*q_ana;
    P[5,2] = 0;
    P[6,2] = 0;
    P[1,3] = p_icoh*q_icoh*q_ana;
    P[2,3] = q_coh*q_coh*q_ana;
    P[3,3] = p_coh*p_coh*q_ana;
    P[4,3] = p_icoh*q_icoh*q_ana;
    P[5,3] = 0;
    P[6,3] = 0;
    P[1,4] = q_icoh*q_icoh*q_ana;
    P[2,4] = p_coh*q_coh*q_ana;
    P[3,4] = p_coh*q_coh*q_ana;
    P[4,4] = p_icoh*p_icoh*q_ana;
    P[5,4] = 0;
    P[6,4] = 0;
    P[1,5] = p_ana;
    P[2,5] = p_ana;
    P[3,5] = p_ana;
    P[4,5] = p_ana;
    P[5,5] = q_reversal;
    P[6,5] = p_revtoana;
    P[1,6] = 0;
    P[2,6] = 0;
    P[3,6] = 0;
    P[4,6] = 0;
    P[5,6] = p_reversal;
    P[6,6] = q_revtoana;
    return P;
  }
  matrix odeUpdateMatrix(vector th) {
    matrix[2,8] M;
    M[1,1] = -th[3] - th[2];
    M[1,2] = th[3];
    M[1,3] = -th[5];
    M[1,4] = -th[5];
    M[1,5] = -th[4];
    M[1,6] = -th[4];
    M[1,7] = 0;
    M[1,8] = 0;
    M[2,1] = th[3];
    M[2,2] = -th[3] - th[2];
    M[2,3] = th[5];
    M[2,4] = th[4];
    M[2,5] = th[5];
    M[2,6] = th[4];
    M[2,7] = 0;
    M[2,8] = 0;
    return M;
  }
  vector odeUpdateVector(vector th, real angleTheta){
    vector[2] mu;
    mu[1] = th[3]*th[8]*cos(angleTheta);
    mu[2] = -th[3]*th[8]*cos(angleTheta);
    return mu;
  }

  vector construct_transition_column(real p_icoh, real p_coh,
                                     real p_ana, real p_reversal, real p_revtoana,
                                     int sigma){
    vector[6] P_col;
    real q_icoh = 1-p_icoh;
    real q_coh = 1-p_coh;
    real q_ana = 1-p_ana;
    real q_reversal = 1-p_reversal;
    real q_revtoana = 1-p_revtoana;
    if (sigma==1){
    P_col = [p_icoh*p_icoh*q_ana,
             p_coh*q_coh*q_ana,
             p_coh*q_coh*q_ana,
             q_icoh*q_icoh*q_ana,
             0,
             0]';
    } else if(sigma==2){
    P_col = [p_icoh*q_icoh*q_ana,
              p_coh*p_coh*q_ana,
              q_coh*q_coh*q_ana,
              p_icoh*q_icoh*q_ana,
              0,
              0]';
    } else if (sigma==3){
    P_col = [p_icoh*q_icoh*q_ana,
              q_coh*q_coh*q_ana,
              p_coh*p_coh*q_ana,
              p_icoh*q_icoh*q_ana,
              0,
              0]';
    } else if (sigma==4){
    P_col = [q_icoh*q_icoh*q_ana,
              p_coh*q_coh*q_ana,
              p_coh*q_coh*q_ana,
              p_icoh*p_icoh*q_ana,
              0,
              0]';
    } else if (sigma==5){
    P_col = [p_ana,
              p_ana,
              p_ana,
              p_ana,
              1,
              p_revtoana]';
    } else if (sigma==6){
    P_col = [0,
              0,
              0,
              0,
              p_reversal,
              q_revtoana]';              
    } else {
      P_col=rep_vector(0,6);
    }
    return(P_col);
  }

  real myDiff(vector v){
    real w;
    w = fabs(v[1] - v[2]);
    return w;
  }
}
data {
  int<lower=1> nTracks;
  int<lower=1> T;
  int<lower=0> T0[nTracks];
  int<upper=T> T1[nTracks];
  vector[T] cos_phi[nTracks]; //the angle to metaphase plate is phi
  matrix[T,2] y[nTracks]; //trajectory data from each sister
  real t_ana_input[nTracks]; //input for prior based on changept model
  real dt;
}
transformed data {
  int nStates = 6;
  real b; //sharpness of sigmoid switch to anaphase. Fixed not inferred.
  vector[2] x0[nTracks];
  vector[nStates] sigma0[nTracks];
  for (i in 1:nTracks){
    x0[i] = y[i][T0[i]+1, ]';
    sigma0[i] = [0, 0.5, 0.5, 0, 0, 0]';
  }
    b = dt/2;
}
parameters {
  //assume for now that biophysical parameters are individual for all trajectories in a cell
  //add t_ana and v_ana as additional parameters in theta
  //shared parameters are the switching rates between states

  real<lower=0> tau[nTracks];
  real<lower=0> alpha[nTracks];
  real<lower=0> kappa[nTracks];
  real<upper=0> v_minus[nTracks];
  real<lower=0> v_plus[nTracks];
  real<lower=0,upper=T*dt> t_ana[nTracks]; //estimate of time of anaphase
  real<lower=0> v_ana[nTracks];
  real<lower=0> L[nTracks];

  real<lower=0,upper=1> p_icoh;
  real<lower=0,upper=1> p_coh;
  real<lower=0,upper=1> p_reversal; // then consider using a strong prior and or hierarchical model
  real<lower=0,upper=1> p_revtoana;

//  vector[10] mu[nTracks]; //for hierarchical model before transforms
}
transformed parameters{
  vector[10] theta[nTracks]; //for metaphase ODE model
  matrix[T, nStates] eta[nTracks];
  matrix[T, nStates] xi[nTracks];
  vector[T] f[nTracks];
  matrix[nStates,nStates] P;
  vector[nStates+2] auxStates;
  vector[T] p_ana[nTracks];

  // fill in etas
  for (ii in 1:nTracks){
//non centering and transform to appropriate space
  theta[ii][1] = tau[ii];
  theta[ii][2] = alpha[ii];
  theta[ii][3] = kappa[ii];
  theta[ii][4] = v_minus[ii];
  theta[ii][5] = v_plus[ii];
  theta[ii][6] = p_coh;
  theta[ii][7] = p_icoh;
  theta[ii][8] = L[ii];
  theta[ii][9] = v_ana[ii];
  theta[ii][10] = t_ana[ii]; //t_ana \in [0,T*dt]

for (t in 1:T0[ii]) {
  for (j in 1:nStates) {
    eta[ii][t,j] = 0;
    xi[ii][t,j] = 0;
  }
  f[ii][t] = 1;
  p_ana[ii][t] = 1/(1+exp(-(t*dt-t_ana[ii])/b));
}
for (t in (T1[ii]+1):T) {
  for (j in 1:nStates) {
    eta[ii][t,j] = 0;
    xi[ii][t,j] = 0;
  }
  f[ii][t] = 1;
  p_ana[ii][t] = 1/(1+exp(-(t*dt-t_ana[ii])/b));
}

 for (t in (T0[ii]+1):T1[ii]) {
    if(t==(T0[ii]+1)) {
      for (j in 1:nStates){
        eta[ii][t,j] = exp(normal_lpdf(y[ii][t]| x0[ii], sqrt(dt/tau[ii])));
      }
    } else {
    for (j in 1:4){
            //construct vector from appending y[t-1] and hidden state eg [0 1 0 0]
      for (i in 1:(nStates+2)){
        if (i<=2) {
          auxStates[i] = y[ii][t-1,i];
        } else if (i==(j+2)) {
          auxStates[i] = 1.0;
        } else {
          auxStates[i] = 0.0;
        }
      }
        eta[ii][t,j] = exp(normal_lpdf(y[ii][t]' | y[ii][t-1]' + dt*odeUpdateMatrix(theta[ii])*auxStates + dt*odeUpdateVector(theta[ii], cos_phi[ii][t]), sqrt(dt/tau[ii])));
      }
      eta[ii][t,5] = exp(normal_lpdf(y[ii][t] | y[ii][t-1] + dt*[v_ana[ii], -v_ana[ii]], sqrt(dt/tau[ii])));
      eta[ii][t,6] = exp(normal_lpdf(y[ii][t] | y[ii][t-1], sqrt(4*dt/tau[ii]))); //pure diffusion in the reversal state, with additional noise
    }
  }
      // work out likelihood contributions
  for (t in (T0[ii]+1):T1[ii]) {
    p_ana[ii][t] = 1/(1+exp(-(t*dt-t_ana[ii])/b));
    P = constructTransitionMatrix(p_icoh,p_coh,p_ana[ii][t],p_reversal,p_revtoana);
    // for the first observation
    if (t==(T0[ii]+1)) {
      //replacing xi[t-1] by sigma0
      f[ii][t] = dot_product((sigma0[ii]')*P, eta[ii][t]');
      xi[ii][t] = (((sigma0[ii]')*P) .* eta[ii][t]) ./ f[ii][t];
    } else {
    // and for the rest
      f[ii][t] = dot_product(xi[ii][t-1]*P, eta[ii][t]');
      xi[ii][t] = ((xi[ii][t-1]*P) .* eta[ii][t]) ./ f[ii][t];
    }
  }
  }
}
model {
  p_icoh ~ beta(2,1);
  p_coh ~ beta(2.5,1);
  p_reversal ~ beta(1,100);
  p_revtoana ~ beta(1,10);
  for (ii in 1:nTracks){
    //priors for individual parameters based on hierarchical model
    tau[ii] ~ gamma(0.5,1.0/1000);
    alpha[ii] ~ normal(0.01,0.1) T[0,];
    kappa[ii] ~ normal(0.05,0.1) T[0,];
    v_minus[ii] ~ normal(-0.03,0.1) T[,0];
    v_plus[ii] ~ normal(0.03,0.1) T[0,];
    v_ana[ii] ~ normal(0.03,0.1) T[0,];
    L[ii] ~ normal(0.790,0.119) T[0,];
    t_ana[ii] ~ normal(t_ana_input[ii],20*dt) T[0,T*dt];
    target += sum(log(f[ii]));
  }
}
generated quantities {
//sampling from hidden states could go here
  int sigma_sim[nTracks,T];
  vector[nStates] state_probs[nTracks];
  vector[nStates] conditional_state_probs; //on log scale
  vector[nStates] P_col;
  int frame;
  for (ii in 1:nTracks) {
      state_probs[ii] = xi[ii][T1[ii]]';
    for (t in (T1[ii]+1):T) {
      sigma_sim[ii,t] = 0; //for time points with missing data
    }
    for (t in 1:T0[ii]) {
      sigma_sim[ii,t] = 0; //for time points with missing data
    }
    //start with the final time pt
    sigma_sim[ii,T1[ii]] = categorical_rng(state_probs[ii]);
    for (t in 1:(T1[ii]-T0[ii]-1)){
      frame = T1[ii]+1-t; //from T1 to T1+1-(T1-T0-1)=T0+2
      P_col = construct_transition_column(p_icoh,p_coh,p_ana[ii][frame],p_reversal,p_revtoana,sigma_sim[ii,frame]);
      conditional_state_probs = log(P_col) + log(xi[ii][frame-1]');
      if (is_inf(sum(conditional_state_probs))){
        sigma_sim[ii,frame-1] = categorical_rng(exp(conditional_state_probs)/sum(exp(conditional_state_probs))); //get around overflow problems with logit this way
      } else {
        sigma_sim[ii,frame-1] = categorical_logit_rng(conditional_state_probs);
      }
    }
  }
}
