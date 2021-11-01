################################################################################
# TITLE: Test Resource Selection Function model for Steller Sea Lion 
# PURPOSE: 
# AUTHOR: Kelly Kapsar 
# (Built of of previous code by Kyle Redilla, RECaP Lab, MSU)
# CREATED: 2021-10-22
# LAST UPDATED ON 2021-10-22
################################################################################
# Open libraries
library(rstan)

# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
homedir <- paste0(getwd(), "/") # Don't know what the difference between workdir and homedir is
workdir <- paste0(getwd(), "/")
# output directory is not unique to date, is equal to resultdir
outdir <- resultdir <- paste(workdir, "results/", sep = "")

################################################################################
# SECTION 1: PREPARE & COMPILE STAN MODEL FILE
################################################################################

modpath <- paste(resultdir, "SSL_RSF_model.stan", sep = "")

write(
  "data{

  int n; // the number of observations
  
  vector[n] y; //
  
  vector[n] x;
  
  real z[n];


}
parameters{

  real beta;

  real<lower = 0> sigma;

}
model {

  y ~ normal(x * beta, sigma);

}"

  , file = modpath)

model <- stan_model(modpath)
# path of the compiled model, saved as an .rda
comp.modpath <- paste(resultdir, "SSL_Test_model.rda", sep = "")
# Write the single R stan_model object to a file for later restoration
saveRDS(model, file = comp.modpath)
