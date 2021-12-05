################################################################################
# TITLE: Steller Sea Lion Discrete Choice Resource Selection Function 
# - Global fixed effects model code 
# AUTHOR: Kelly Kapsar, CSIS 
# CREATED: 2021-11-15
# LAST UPDATED ON 2021-11-15
################################################################################
# Open libraries
library(rstan)

# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
workdir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/"
# homedir <- "~Documents/SSL/" 
# workdir <- "/mnt/scratch/kapsarke/Documents/SSL/"
resultdir <- "Results/SSL_DummyModel_2021-11-20/"

# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(file.path(workdir, resultdir)), 
       dir.create(file.path(workdir, resultdir)), FALSE)

################################################################################
# SECTION 1: PREPARE & COMPILE STAN MODEL FILE
################################################################################

# Set working directory to scratch 
setwd(workdir)
# Specify path to write model to
modpath <- paste(resultdir, "model.stan", sep = "")

# Write the stan model
write(
  "data{
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
"
, file = modpath)

# Compile the model
model <- stan_model(modpath)

# Save the compiled model object
comp.modpath <- paste(workdir, resultdir, "model.rda", sep = "")
saveRDS(model, file = comp.modpath)


