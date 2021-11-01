################################################################################
# TITLE: Test Resource Selection Function model for Steller Sea Lion 
# PURPOSE: 
# AUTHOR: Kelly Kapsar 
# (Built from previous code by Kyle Redilla, RECaP Lab, MSU)
# CREATED: 2021-10-22
# LAST UPDATED ON 2021-10-22
################################################################################
# Open libraries
start <- Sys.time()
# OPEN LIBRARIES 
library(rstan)
library(dplyr)
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
  x[, , 11] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 17])[, 1]
  x[, , 12] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 18])[, 1]
  return(x)
}

################################################################################
# SECTION 1: PREPARE DATA
################################################################################
# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
workdir <- paste0(getwd(), "/")
datadir <- paste(homedir, "Data_Processed/", sep = "")
resultdir <- paste(homedir, "Results/", sep = "")
outdir <- paste(resultdir, "SSL_RSF_Results_",
                datestr, "/", sep = "")

# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(file.path(outdir)), 
       dir.create(file.path(outdir)), FALSE)

# create output file connection
outpath <- paste(outdir, "output.Rout", sep = "")

# change to scratch directory
setwd(workdir)

# Initialize sink - 
# sink(outpath)

# time begin
print("Time begin"); print(start); cat("\n"); 

# resource selection data path
rs.data.path <- paste(datadir, "/Telemetry/UsedAndAvail_WeeklyKDE_20211029.rds", sep = "")
# read in resource selection data and sort with used point at top of list
rs_data <- readRDS(rs.data.path)

################################################################################
# SECTION 2: FIT MODEL
################################################################################
#  For execution on local, multicore CPU with excess RAM
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# model compiled on intel14 node, saved in scratch results
comp.modpath <- paste(resultdir, "SSL_RSF_model.rda", sep = "")
modpath <- paste(resultdir, "SSL_RSF_model.stan", sep="")

# read in compiled model object
mod <- readRDS(comp.modpath)


# Start cluster
cl <- 2
clus <- makeCluster(cl)
registerDoParallel(clus)

# number of coefficients
K <- 12
# number of choices in each set
C <- 6 
# number of individuals
nInd <- 11

# For each elk, sample 4 chains of global model
fitlst <- foreach(j = 1:nInd,.packages = c("rstan")) %dopar% {
  # subset data array for jth model
  ind_rs_data <- rs_data[rs_data$ind_id == j, ]
  #extract id's for each choice set
  cs <- unique(ind_rs_data[, 'choice_id']) 
  # number of choice sets
  N <- length(cs$choice_id) 
  # create design array
  x <- rsf_array(ind_rs_data, c(N, C, K))
  
  # must enter data into a list
  data <- list(
    C = C, K = K, N = N,
    x = x, y = rep(1, N),
    obs=c(1,0,0,0,0,0),
    pos = diag(1, 6) # kk - NOT SURE WHAT Y, OBS, AND POS ARE... 
  )
  
  # initial values are best supplied as a function
  inits <- function(){
    list(
      beta = runif(K, -5, 5)
    )
  }
  # a character vector of parameters to monitor
  params <- c('beta', 'chis_obs', 'chis_sim')
  
  fit <- sampling(mod, data = data, pars = params, init = inits,
                  chains = 4, iter = 1000, warmup = 200, thin = 1)
}
stopCluster(clus)

# save model fit
fitpath <- paste(outdir, "SSL_GlobalIndlRSF_Results_", datestr,
                 ".rda", sep = "")
save(fitlst, file = fitpath)

# print finish time
end <- Sys.time()
print("Time finished"); print(end); cat("\n"); 
end - start; cat("\n")
browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")