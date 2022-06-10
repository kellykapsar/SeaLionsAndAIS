################################################################################
# TITLE: Steller Sea Lion Discrete Choice Resource Selection Function 
# - Global fixed effects code 
# AUTHOR: Kelly Kapsar, CSIS 
# CREATED: 2021-11-15
# LAST UPDATED ON 2021-11-15
################################################################################
start <- Sys.time()
# OPEN LIBRARIES 
library(ggplot2)
library(rstan)
library(dplyr)
library(foreach)
library(doParallel)
library(rethinking)
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
  x[, , 1] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 3])[, 1]
  x[, , 2] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 4])[, 1]
  x[, , 3] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 5])[, 1]
  # x[, , 4] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 10])[, 1]
  # x[, , 5] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 11])[, 1]
  # x[, , 6] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 12])[, 1]
  # x[, , 7] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 13])[, 1]
  # x[, , 8] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 14])[, 1]
  # x[, , 9] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 15])[, 1]
  # x[, , 10] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 16])[, 1]
  # x[, , 11] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 17])[, 1]
  # x[, , 12] <- scale(rs_data[c(ind1, ind2, ind3, ind4, ind5, ind6), 18])[, 1]
  return(x)
}


################################################################################
# SECTION 1: PREPARE DATA
################################################################################

# Initialize sink - 
# sink(outpath)

# time begin
print("Time begin"); print(start); cat("\n"); 

# Set seed
set.seed(101557)

# number of coefficients
K <- 3
# number of choices in each set
C <- 6
# choice sets
N <- 500

# covariate matrix
mX <- matrix(runif(N*C*K), C*N, K) %>% as.data.frame()
mX$choice_id <- rep(1:N, each=C)
mX$used <- 0

# coefficients for each choice
b1 = c(2,0,0,0,0,0)
b2 = c(4,0,0,0,0,0)
b3 = c(-3,0,0,0,0,0)


newchoice <- data.frame()
# vector of probabilities
choices <- list()
for(i in 1:N){
  choiceset <- mX[mX$choice_id == i,]
  score <- b1*choiceset[,1] + b2*choiceset[,2] + b3*choiceset[,3]
  p <- softmax(score)
  choice <- sample(1:C, size=1, prob=p)
  choiceset$used[choice] <- 1
  newchoice <- rbind(newchoice, choiceset)
  choices <- c(choices, choice)
}


# Order data appropriately 
newchoice <- newchoice %>% arrange(choice_id, desc(used))
# Arrange columns
test_data <- newchoice[,c("choice_id", "used", "V1", "V2", "V3")]

################################################################################
# SECTION 2: FIT MODEL
################################################################################

# read in compiled model object
mod <- readRDS("model.rds")


#extract id's for each choice set
N <- length(unique(test_data$choice_id))
# create design array
x <- rsf_array(test_data, c(N, C, K))
# must enter data into a list
data <- list(
  C = C, K = K, N = N,
  x = x, y = rep(1, N),
  obs=c(1,0,0,0,0,0),
  pos = diag(1, 6)
)

# initial values are best supplied as a function
inits <- function(){
  list(
    beta = runif(K, -10, 10)
  )
}
# a character vector of parameters to monitor
params <- c('beta', 'chis_obs', 'chis_sim', "rcat", "rch", "expected")

fit <- sampling(mod, data = data, pars = params, init = inits,
                chains =4, iter = 1000, warmup = 200, thin = 1)


# save model fit
saveRDS(fit, "modelfit.rds")


# Figure 1: plot results of Posterior 
# extract observed discrepancies
  obs.disc <- rstan::extract(fit, 'chis_obs', F)[,,1]
  # extract simulated discrepancies
  sim.disc <- rstan::extract(fit, 'chis_sim', F)[,,1]
  # number of post-warmup draws
  pdraw <- dim(obs.disc)[1] * dim(obs.disc)[2]
  # Calculate Bayesian P-value
  pval <- sum(rstan::extract(fit, 'chis_sim', F)[, , 1] >
                rstan::extract(fit, 'chis_obs', F)[, , 1]) / pdraw
  print(pval)
  
  #  NEW VERION I'M TRYING 
  obs.disc <- rstan::extract(fit, 'chis_obs', T)
  # extract simulated discrepancies
  sim.disc <- rstan::extract(fit, 'chis_sim',T)
  # number of post-warmup draws
  pdraw <- length(obs.disc$chis_obs)
  # Calculate Bayesian P-value
  pval <- sum(obs.disc$chis_obs > sim.disc$chis_sim) / pdraw
  print(pval)

# Messing around trying to figure out how the chi-square test is working 
# See one note from 12/10/21 for more details 
obs <- matrix(c(1,0,0,0,0,0), nrow=500, ncol=6, byrow=TRUE)
pos <- diag(1,6)
expected <- rstan::extract(fit, "expected")
rcat <- rstan::extract(fit, "rcat")
rch <- rstan::extract(fit, "rch")
chio <- rstan::extract(fit, 'chis_obs')
chisim <- rstan::extract(fit, 'chis_sim')


## NOT SURE WHY MY MANUALLY CALCULATED CHI-SQUARE IS NOT FITTING THE DATA 

# [1] 1 0 0 0 0 0
expected[[1]][1:50,1,1:6]
#0.2262241 0.2082138 0.1080995 0.1140147 0.1593858 0.1840622
obs - expected[[1]][1,1,1:6]
#0.7737759 -0.2082138 -0.1080995 -0.1140147 -0.1593858 -0.1840622
(obs - expected[[1]][1,1,1:6])**2
#0.59872911 0.04335298 0.01168550 0.01299934 0.02540383 0.03387888
(obs - expected[[1]][1,1,1:6])**2/expected[[1]][1,1,1:6]
#2.6466192 0.2082138 0.1080995 0.1140147 0.1593858 0.1840622
sum((obs[1,] - expected[[1]][1,1,1:6])**2/expected[[1]][1,1,1:6])

sum((obs[1:500,] - expected[[1]][1,1:500,1:6])**2/expected[[1]][1,1:500,1:6])

obs.disc[1]
chio_permT <- rstan::extract(fit, 'chis_obs', T)
chio_permT$chis_obs[1]
chio[[1]][1]


#3.420395
rch[[1]][1,1,1:6]
#0 0 0 0 0 1
sum((rch[[1]][1,600,1:6] - expected[[1]][1,600,1:6])**2/expected[[1]][1,600,1:6])

sum((rch[[1]][1,1:500,1:6] - expected[[1]][1,1:500,1:6])**2/expected[[1]][1,1:500,1:6])
chisim[[1]][1]

  


hist(expected[[1]][,2,2], col="blue", xmin=0, xmax=0.5)
hist(expected[[1]][,2,3], add=T, col="red")
hist(expected[[1]][,2,4], add=T, col="green")
hist(expected[[1]][,2,5], add=T, col="yellow")
hist(expected[[1]][,2,6], add=T, col="purple")
hist(expected[[1]][,2,1], add=T)


# print finish time
end <- Sys.time()
print("Time finished"); print(end); cat("\n"); 
end - start; cat("\n")



