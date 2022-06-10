data{
    int C;  // the number of choices per set
    int K;  // the number of slope parameters per individual
    int N;  // the number of observations
    matrix[C, K] x[N];  // variables in design matrix format
    int y[N]; // index of which alternative was selected
    vector[C] obs; //variable used in test-stat calculation
    vector[C] pos[C];
  }
  
  parameters{
    vector[K] beta;  // rsf coefficients
  }
  
  model{
    
    for(l in 1:K){
      beta[l] ~ normal(0, 10);  // prior distribution for slope parameter mean
    }
  
    for(i in 1:N){
      y[i] ~ categorical_logit(x[i] * beta);
    }
  }
  
  generated quantities{
    simplex[C] expected[N]; //the probabilities of use from dc_model
    vector[N] chis_obs_i;   //chi-square value for each choice set using observed data
    vector[N] chis_sim_i;  // chi-square value from simulated data
    real chis_obs;  //sum of chi square values across all choice sets
    real chis_sim;  //sum of chi square values across all choice sets
    vector[C] rch[N];  //simulated random choice of used alt. within obs. data sets
    vector[N] log_lik; // log likelihood
    int rcat[N];
  
    for(a in 1:N){
      expected[a] = softmax(x[a] * beta); //??Create prob. of use of each used loc 
  
      rcat[a] = categorical_rng(expected[a]); //??simulate 0 or 1 based on covar val*beta at each point
      rch[a] = pos[rcat[a]];
  
      chis_obs_i[a] = sum(((obs - expected[a]) .* (obs - expected[a])) ./
                         expected[a]);
      chis_sim_i[a] = sum(((rch[a] - expected[a]) .* (rch[a] - expected[a])) ./
                         expected[a]);
    
      log_lik[a] = categorical_logit_lpmf(y[a]| x[a] * beta);
    }
    chis_obs = sum(chis_obs_i);
    chis_sim = sum(chis_sim_i);
  
  }

