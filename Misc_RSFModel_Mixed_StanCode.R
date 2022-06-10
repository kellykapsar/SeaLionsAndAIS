################################################################################
# TITLE: Steller Sea Lion Discrete Choice Resource Selection Function 
# - Global mixed effects model code 
# AUTHOR: Kelly Kapsar, CSIS 
# CREATED: 2021-11-15
# LAST UPDATED ON 2021-11-15
################################################################################
# Open libraries
library(rstan)

# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")


################################################################################
# SECTION 1: PREPARE & COMPILE STAN MODEL FILE
################################################################################


# Specify path to write model to
modpath <- "../Results/mixedmodel.stan"

write(
  "data{
  int C;  // the number of choices per set
  int K;  // the number of slope parameters per individual
  int N;  // the number of observations
  int nWks;  // the number of individuals
  int Wks[N];  // indexes the indivial associated with each observation
  matrix[C, K] x[N];  // variables in design matrix format
  int y[N]; // index of which alternative was selected
  vector[C] obs;
  vector[C] pos[C];

  }
  
  parameters{
  vector[K] beta[nWks];  // slope parameters for individual elk
  vector[K] mu;  // mean for each slope parameter
  vector<lower = 0>[K] stdev;  // standard deviation of each slope parameter
  }
  
  model{
  // vector[C] psi[N];
  
  //for(l in 1:K){           // Attempting to vectorize, removing loop
  mu ~ normal(0, 10);  // prior distribution for slope parameter mean
  stdev ~ normal(0, 10); // prior distribution for slope parameter sd
  //}
  
  // for(j in 1:nWks){
  for (k in 1:K){    
  beta[nWks] ~ normal(mu, stdev);
  }
  // } 
  
  for(i in 1:N){
  y ~ categorical_logit(x[i] * beta[Wks[i]]);
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
    expected[a] = softmax(x[a] * beta[Wks[a]]);

    rcat[a] = categorical_rng(expected[a]);
    rch[a] = pos[rcat[a]]; 
    chis_obs_i[a] = sum(((obs - expected[a]) .* (obs - expected[a])) ./
                        expected[a]);
    chis_sim_i[a] = sum(((rch[a] - expected[a]) .* (rch[a] - expected[a])) ./
                        expected[a]);
  }
  chis_obs = sum(chis_obs_i);
  chis_sim = sum(chis_sim_i);
  
  }"
, file = modpath)

model <- stan_model(modpath)
comp.modpath <- "../Results/mixedmodel.rds"
saveRDS(model, file = comp.modpath)


