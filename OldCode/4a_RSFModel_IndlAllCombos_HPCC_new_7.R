################################################################################
# TITLE: Resource selection function analysis 
# Modified by Kelly Kapsar from script by Kyle Redilla 
# November 2021 
################################################################################
start <- Sys.time()
# OPEN LIBRARIES 
library(ggplot2)
library(rstan)
library(dplyr)
library(foreach)
library(doParallel)
library(loo)
# Define functions
# Create RSF design array
#
# This function creates a design array for resource selection function
#   with no interaction terms  
#
# Args:
#   data: The data.frame object containing rs data
#
#   dims: A vector containing the number of observations, choices per set, and
#           selection coefficients
#
# Returns:
#   X: An N * C * K dimensional design array containing rs data prepared
#              for ERS1_step5dev_code stan model
#
rsf_array <- function(data, dims){
  N <- dims[1]
  C <- dims[2]
  K <- dims[3]
  # Indices for filling design array
  ind1 <- seq(1, N * C, by = 6); ind2 <- ind1 + 1; ind3 <- ind2 + 1
  ind4 <- ind3 + 1; ind5 <- ind4 + 1; ind6 <- ind5 + 1
  rs_data <- data
  # design array
  x <- array(dim = c(N, C, K))
  # 7 = column number at which covariates start
  # "scale" centers the data (mean 0, sd 1)
  x[, , 1] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 7])[, 1]
  x[, , 2] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 8])[, 1]
  x[, , 3] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 9])[, 1]
  x[, , 4] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 10])[, 1]
  x[, , 5] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 11])[, 1]
  x[, , 6] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 12])[, 1]
  x[, , 7] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 13])[, 1]
  x[, , 8] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 14])[, 1]
  x[, , 9] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 15])[, 1]
  x[, , 10] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 16])[, 1]
  # x[, , 11] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 17])[, 1]
  # x[, , 12] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 18])[, 1]
  return(x)
}


################################################################################
# SECTION 1: PREPARE DATA
################################################################################
# Set seed 
set.seed(101557)

# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")


# create output file connection
outpath <- paste("../Results/output.Rout", sep = "")

# Initialize sink - 
# sink(outpath)

# time begin
print("Time begin"); print(start); cat("\n"); 

# resource selection data path
rs.data.path <- paste("../Data/UsedAndAvail_WeeklyKDE_20220311.rds", sep = "")
# read in resource selection data and sort with used point at top of list
rs_data <- readRDS(rs.data.path) %>% group_by(weeklyhr_id) %>% mutate(weeklyhr_id = cur_group_id())

# Create list of covariate names 
covar_names <- c("bathymetry", "dist_land", "dist_500m", "slope", "sst", "wind", "logship", 
                 "logfish", "prox_fish_km", "prox_ship_km")

# number of coefficients
K <- 10
# number of choices in each set
C <- 6 
# number of individuals 
nInd <- 11


# Create grid with all TRUE/FALSE covariate combinations
v1 <- rep(TRUE, K)
v1 <- lapply(v1, append, FALSE)
master <- expand.grid(v1)
# Re-order master based on number of variables
master <- master[order(rowSums(master), decreasing = TRUE),]
master <- master[-which(rowSums(master) == 0),] # Remove no predictor model



################################################################################
# SECTION 2: FIT MODEL
################################################################################

# model compiled on intel14 node
modpath <- "../Results/model.rda"
mod <- readRDS(modpath)

#  For execution on local, multicore CPU with excess RAM
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Start cluster
# cl <- 28
# clus <- makeCluster(cl)
# registerDoParallel(clus)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

# RUN ON ONE INDIVIDUAL
ind <- 7
rs_data_subset <- rs_data[rs_data$ind_id == ind,]

# Create filepath for output
# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(paste0("../Results/SSL_Ind", ind, "/")), 
       dir.create(paste0("../Results/SSL_Ind", ind, "/")), FALSE)

# Calculate number of choices and create data array
N <- length(unique(rs_data_subset$choice_id))
x <- rsf_array(rs_data_subset, c(N, C, K))

# Calculate pairwise correlations among all covars 
corrs <- read.csv(paste0("../Data/Corrs_", ind, ".csv"))

# Adjust all possible combinations matrix to remove rows with correlated variables
for(i in 1:length(corrs$x)){
  truex <- which(master[,corrs$x[i]] == TRUE) # row nums where first correlated var included
  truey <- which(master[,corrs$y[i]] == TRUE) # row nums where second correlated var included
  trueboth <- intersect(truex, truey) # row nums where both correlated covars included
  master <- master[-trueboth,] # Remove rows including both correlated covars
}

master <- master[-which(rowSums(master) == 1),]

write.csv(master, paste0("../Results/ModelVariableList_", ind, ".csv"))

fitlst <- c()

# loop through the models and fit each one
fitlst <- foreach(j = 1:length(master$Var1),.packages = c("rstan", "loo")) %dopar% {
  keep.vars <- unlist(master[j,])
  x.temp <- x[,,keep.vars]
  K <- sum(keep.vars)
  # must enter data into a list
  data <- list(
    C = C, K = K, N = N, 
    x = x.temp,
    y = rep(1, N)
  )
  
  # initial values are best supplied as a function
  # initial values are best supplied as a function
  inits <- function(){
    list(
      beta = runif(K, -10, 10)
    )
  }
  
  # a character vector of parameters to monitor
  #   monitor log_lik only to save space and time
  #   all that is needed to calculate WAIC
  params <- c('beta', 'log_lik')
  # sample from the model object
  fit1 <- sampling(mod, data = data, pars = params, init = inits,
                   chains = 4, iter = 1000, warmup = 200, thin = 1)
  
  # Add leading zero so that model files order correctly 
  ifelse(nchar(j) == 1, jnum <- paste0("00",j), ifelse(nchar(j) == 2, jnum <- paste0("0", j), jnum <- j))
  
  save(fit1, file=paste0("../Results/SSL_Ind", ind, "/Ind_", ind, "_Model_", jnum, ".rda"))
  
  # Check parameter estimates and convergence diagnostics
  paramsums <- data.frame(summary(fit1, pars=c("beta"))$summary)
  
  # Estimate model fit using LOO-CV
  log_lik1 <- extract_log_lik(fit1, merge_chains = FALSE)
  rel_n_eff <- relative_eff(exp(-log_lik1))
  loo1 <- loo(log_lik1, r_eff = rel_n_eff, cores = 4)
  waic1 <- waic(log_lik1)
  
  return(list(paramsums, loo1, waic1))
}

# save model fit
filename <- paste('SSL_IndlAllCombos_Ind', ind, '_results_', datestr,
                  '_new.rda', sep = '')
filepath <- paste("../Results/", filename, sep = '')
save(fitlst, file = filepath)


# browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")

# print finish time
end <- Sys.time()
print("Time finished"); print(end); cat("\n"); 
end - start; cat("\n")

sink()
