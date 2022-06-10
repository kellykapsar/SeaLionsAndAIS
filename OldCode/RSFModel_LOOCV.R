################################################################################
# TITLE: Individual Missouri elk resource selection function analysis: Step 5 -
#   Code for fitting individual models, RSF2
# PURPOSE: This code is for fitting the random slopes individual resource
#   selection functions (model RSF2) to the Missouri elk relocation data, i.e. a 
#   model without indicator variable selection. We use Stan 
#   because it is much more efficient than Gibbs sampling for this scenario. 
#   This script is meant to be run on my HPCC ICER account.
# AUTHOR: Kyle Redilla, RECaP Lab
# CREATED: 2017-01-09
# LAST UPDATED ON 2017-01-15
################################################################################
start <- Sys.time()
# OPEN LIBRARIES 
library(ggplot2)
library(rstan)
library(dplyr)
library(foreach)
library(doParallel)
library(loo)


################################################################################
# SECTION 1: PREPARE DATA
################################################################################
# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")


# create output file connection
outpath <- paste("../Results/output.Rout", sep = "")

#Specify individual
ind <- 1

# Read in model rrseults (object name = fitlst)
load(paste0("../Results/SSL_IndlAllCombos_", ind, "_results_2021-11-13.rda"))

# time begin
print("Time begin"); print(start); cat("\n"); 

#  For execution on local, multicore CPU with excess RAM
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Start cluster
# cl <- 28
# clus <- makeCluster(cl)
# registerDoParallel(clus)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

loolist <- foreach(i = 1:5, .packages=c(loo)) %dopar% {
    log_lik1 <- extract_log_lik(fitlst[[i]], merge_chains = FALSE)
    rel_n_eff <- relative_eff(exp(-log_lik1))
    loo(log_lik1, r_eff = rel_n_eff, cores = 4)
}

# Save list of LOO values
save(loolist, file=paste0("../Results/SSL_IndlAllCombos_", ind, "_loo_", datestr, "_TEST.rda"))