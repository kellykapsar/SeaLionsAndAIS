################################################################################
# TITLE: Test Resource Selection Function model for Steller Sea Lion 
# PURPOSE: 
# AUTHOR: Kelly Kapsar 
# (Built of of previous code by Kyle Redilla, RECaP Lab, MSU)
# CREATED: 2021-10-22
# LAST UPDATED ON 2021-10-22
################################################################################
# Open libraries
start <- Sys.time()
# OPEN LIBRARIES 
library(rstan)
library(dplyr)


################################################################################
# SECTION 1: PREPARE DATA
################################################################################
# Specify directories
datestr <- format(Sys.time(), "%Y-%m-%d")
homedir <- paste0(getwd(), "/") # kk - Don't know what the difference between workdir and homedir is
workdir <- paste0(getwd(), "/")
datadir <- paste(workdir, "data/", sep = "")
resultdir <- paste(workdir, "results/", sep = "")
outdir <- paste(resultdir, "ERS1_step4_results_",
                datestr, "/", sep = "")

# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(file.path(outdir)), 
       dir.create(file.path(outdir)), FALSE)

# create output file connection
outpath <- paste(outdir, "output.Rout", sep = "")

# change to scratch directory
setwd(workdir)


################################################################################
# SECTION 2: FIT MODEL
################################################################################
#  For execution on local, multicore CPU with excess RAM
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# model compiled on intel14 node, saved in scratch results
comp.modpath <- paste(resultdir, "SSL_Test_model.rda", sep = "")
modpath <- paste(resultdir, "SSL_RSF_model.stan", sep="")

# read in compiled model object
mod <- readRDS(comp.modpath)


beta <- 0.2

sigma <- 0.5

x <- -1:200

data <- list(y = x * beta + rnorm(length(x), 0, sigma), n = length(x), x = x,
             z = x, g = x)

plot(x,data$y)



scratch <- stan(file = modpath,
                iter = 2000, 
                warmup = 1000,
                data = data)

plot(scratch)


# print finish time
end <- Sys.time()
print("Time finished"); print(end); cat("\n"); 
end - start; cat("\n")
