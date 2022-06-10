################################################################################
# TITLE: Individual Missouri elk resource selection function analysis: Step 5dev
#   - Code for compiling the individual RSF2 Stan model on the HPCC
# PURPOSE: This code is for compiling the Stan model on the HPCC using the dev-
#   intel14 node, so that either intel14 or intel16 nodes can run the job
#   ERS1_step5_code_2017-01-12.qsub
# AUTHOR: Kyle Redilla, RECaP Lab
# CREATED: 2017-1-9
# LAST UPDATED ON 2017-1-15
################################################################################
# Open libraries
library(rstan)

# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
# homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
# workdir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/"
homedir <- "."
workdir <- "."
datadir <- paste(workdir, "/Data/", sep = "")
resultdir <- paste("/Results/SSL_IndlGlobalMixedEffects_",
                   datestr, "/",sep = "")
# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(paste0(workdir, resultdir)),
       dir.create(paste0(workdir, resultdir)), FALSE)

################################################################################
# SECTION 1: PREPARE & COMPILE STAN MODEL FILE
################################################################################

#  For execution on local, multicore CPU with excess RAM
rstan_options(auto_write = TRUE)

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
  int nWks;  // the number of weeks
  int wk[N];  // indexes the week associated with each observation
  matrix[C, K] x[N];  // variables in design matrix format
  int y[N]; // index of which alternative was selected
  vector[C] obs;
  vector[C] pos[C];

  }
  
  parameters{
  vector[K] beta[nWks];  // slope parameters for each week
  real mu[K];  // mean for each slope parameter
  real<lower = 0> stdev[K];  // standard deviation of each slope parameter
  }
  
  model{
  
//  for(k in 1:K){
//    mu[k] ~ normal(0, 1);  // prior distribution for slope parameter mean
//    for(j in 1:nWks)
//      beta[1,j] ~ normal(mu[k], stdev[k]);
//  }
  
  mu ~ normal(0, 1);  // prior distribution for slope parameter mean
  stdev ~ normal(0, 1); // prior distribution for slope parameter sd

  for (k in 1:K){    
    beta[nWks] ~ normal(mu, stdev);
  }
  
  for(i in 1:N){
    y ~ categorical_logit(x[i] * beta[wk[i]]);
  }
  }
  
  generated quantities{
  simplex[C] expected[N]; //probabilities of use from dc model
  vector[N] chis_obs_i; //chi-square value for each choice based on observations
  vector[N] chis_sim_i; // chi-square value from simulated data
  real chis_obs; //sum of chi-square values
  real chis_sim;
  vector[C] rch[N]; //simulated random choice of used alternative
  int rcat[N];

  
  for(a in 1:N){
    expected[a] = softmax(x[a] * beta[wk[a]]);

    rcat[a] = categorical_rng(expected[a]);
    rch[a] = pos[rcat[a]]; 
    chis_obs_i[a] = sum(((obs - expected[a]) .* (obs - expected[a])) ./
                        expected[a]);
    chis_sim_i[a] = sum(((rch[a] - expected[a]) .* (rch[a] - expected[a])) ./
                        expected[a]);
  }
  chis_obs = sum(chis_obs_i);
  chis_sim = sum(chis_sim_i);
  
  }

"
, file = modpath)

# Compile the model
model <- stan_model(modpath)
print(model)

# Save the compiled model object
comp.modpath <- paste(workdir, resultdir, "model.rda", sep = "")
saveRDS(model, file = comp.modpath)


