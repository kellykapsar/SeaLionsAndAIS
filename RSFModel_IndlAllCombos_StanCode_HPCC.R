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

################################################################################
# SECTION 1: PREPARE & COMPILE STAN MODEL FILE
################################################################################


# Specify path to write model to
modpath <- "../Results/model.stan"

# Write the stan model
write(
  "data{
    int C;  // the number of choices per set
    int K;  // the number of slope parameters per individual
    int N;  // the number of observations
    matrix[C, K] x[N];  // variables in design matrix format
    int y[N]; // index of which alternative was selected
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
    vector[N] log_lik; // log likelihood
  
    for(a in 1:N){
      log_lik[a] = categorical_logit_lpmf(y[a]| x[a] * beta);
    }

  }
"
  , file = modpath)

# Compile the model
model <- stan_model(modpath)

# Save the compiled model object
comp.modpath <- "../Results/model.rda"
saveRDS(model, file = comp.modpath)




