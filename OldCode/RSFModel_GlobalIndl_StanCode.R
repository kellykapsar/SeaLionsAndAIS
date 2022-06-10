################################################################################
# TITLE: Individual Missouri elk resource selection function analysis: Step 4dev
#   - Compiling individual global stan model 
# PURPOSE: This code is for compiling the Stan model on the HPCC using the dev-
#   intel14 node, so that either intel14 or intel16 nodes can run the job.
#   This model is the bare bones discrete choice model, i.e. no random effects.
#   We generate the chi-square statistic for each observation, to validate the 
#   model using Bayesian p-value later on.
# AUTHOR: Kyle Redilla, RECaP Lab
# CREATED: 2017-02-02
# LAST UPDATED ON 2017-02-02
################################################################################
# Open libraries
library(rstan)

# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
workdir <- paste0(getwd(), "/")
datadir <- paste(homedir, "Data_Processed/", sep = "")
resultdir <- paste(homedir, "Results/", sep = "")
outdir <- paste(resultdir, "SSL_RSF_GlobalIndlResults_",
                datestr, "/", sep = "")

################################################################################
# SECTION 1: PREPARE & COMPILE STAN MODEL FILE
################################################################################

modpath <- paste(resultdir, "SSL_RSF_model.stan", sep = "")

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
  vector[K] beta;  // slope parameters 
  }
  
  model{

  for(k in 1:K){
  beta[k] ~ normal(0, 5);  // prior distribution for slope parameter 
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
  int rcat[N];
  
  for(a in 1:N){
  expected[a] <- softmax(x[a] * beta);
  
  rcat[a] <- categorical_rng(expected[a]);
  rch[a] <- pos[rcat[a]];
  
  chis_obs_i[a] <- sum(((obs - expected[a]) .* (obs - expected[a])) ./
  expected[a]);
  chis_sim_i[a] <- sum(((rch[a] - expected[a]) .* (rch[a] - expected[a])) ./
  expected[a]);
  }
  chis_obs <- sum(chis_obs_i);
  chis_sim <- sum(chis_sim_i);
  
  }"
, file = modpath)

model <- stan_model(modpath)
# path of the compiled model, saved as an .rda
comp.modpath <- paste(resultdir, "SSL_RSF_GlobalIndl_model.rda", sep = "")
# Write the single R stan_model object to a file for later restoration
saveRDS(model, file = comp.modpath)

