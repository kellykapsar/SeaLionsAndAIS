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
# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
workdir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/"
datadir <- paste(homedir, "Data_Processed/", sep = "")
resultdir <- paste("Results/SSL_IndlAllCombos_",
                   datestr, "/", sep = "")
# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(file.path(workdir, resultdir)), 
       dir.create(file.path(workdir, resultdir)), FALSE)

# create output file connection
outpath <- paste(resultdir, "output.Rout", sep = "")

# change to scratch directory
setwd(workdir)

# Initialize sink - 
# sink(outpath)
# time begin
print("Time begin"); print(start); cat("\n"); 

# resource selection data path
rs.data.path <- paste(datadir, "/Telemetry/UsedAndAvail_WeeklyKDE_20211104.rds", sep = "")
# read in resource selection data and sort with used point at top of list
rs_data <- readRDS(rs.data.path) %>% group_by(weeklyhr_id) %>% mutate(weeklyhr_id = cur_group_id())



# number of coefficients
K <- 10
# number of choices in each set
C <- 6 
# number of individuals 
nInd <- 11

################################################################################
# SECTION 2: FIT MODEL
################################################################################

# model compiled on intel14 node
modpath <- paste(resultdir, "model.rda", sep = "")
mod <- readRDS(modpath)

#  For execution on local, multicore CPU with excess RAM
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Start cluster
# cl <- 11
# clus <- makeCluster(cl)
# registerDoParallel(clus)

# Create grid with all TRUE/FALSE covariate combinations
v1 <- rep(TRUE, K)
v1 <- lapply(v1, append, FALSE)
master <- expand.grid(v1)
# Re-order master based on number of variables
master <- master[order(rowSums(master), decreasing = TRUE),]
# loop through the models and fit each one
#   (preliminary analysis suggested that models with fewer variables
#    performed worse according to WAIC, so only fitting first 512)

# TEST RUN ON ONE INDIVIDUAL
ind <- 1
rs_data_subset <- rs_data[rs_data$ind_id == ind,]
N <- length(unique(rs_data_subset$choice_id))
x <- rsf_array(rs_data_subset, c(N, C, K))

for(i in 1:5){
  print(i)
  keep.vars <- unlist(master[i,])
  x.temp <- x[,,keep.vars]
  K <- dim(x.temp)[3]
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
  params <- c('log_lik')
  # sample from the model object
  fit1 <- sampling(mod, data = data, pars = params, init = inits,
                   chains = 4, iter = 1000, warmup = 200, thin = 1)
  fit <- list(fit1, keep.vars)
  # save model fit
  filename <- paste('SSL_IndlAllCombos_', ind, "_", i, '_results_', datestr,
                    '.rda', sep = '')
  filepath <- paste(resultdir, filename, sep = '')
  save(fit, file = filepath)
}

browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")

# print finish time
end <- Sys.time()
print("Time finished"); print(end); cat("\n"); 
end - start; cat("\n"); 
