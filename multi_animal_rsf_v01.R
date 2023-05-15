setwd("C:/Users/12487/OneDrive/sea_lion")

d <- readRDS("UsedAndAvail_WeeklyKDE_20220706.rds") 

all <- d %>%
  as_tibble() %>% 
  group_by(ind_id, choice_id) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() %>% 
  group_by(ind_id) %>% 
  mutate(smallest_set = min(set)) %>% 
  mutate(set = ( set - smallest_set) + 1 ) %>% 
  group_by(ind_id, set) %>% 
  mutate( choice = row_number()) %>%
  mutate( ind_set = cur_group_id()) %>% 
  ungroup() %>% 
  mutate( row_id = row_number(),
          used = as.numeric(used) - 1,
          prox_fish_km = as.numeric(scale(prox_fish_km))) %>% 
  group_by(ind_set) %>% 
  mutate( start = min(row_id), 
          end = max(row_id)) %>% 
  dplyr::select( ind = ind_id,
                 set, 
                 ind_set,
                 start, 
                 end, 
                 choice,
                 used,
                 prox_fish_km ) 

only_used <- all %>% 
  filter(used == 1) %>% 
  dplyr::select( ind, set, ind_set, start, end, used)

data <- list(
  used = only_used$used, 
  prox_fish = all$prox_fish_km
)

constants <- list(
  n_animals = length(unique(all$ind)),
  n_used = length(data$used), 
  n_all = length(data$prox_fish),
  start = only_used$start, 
  end = only_used$end, 
  individual = all$ind, 
  ind_set = all$ind_set
)

code <- nimble::nimbleCode({
  
  # Priors for population-level hyperparameters
  mu_beta0 ~ dnorm(0, sd = 2)
  mu_beta1 ~ dnorm(0, sd = 2)
  sd_beta0 ~ dexp(1)
  sd_beta1 ~ dexp(1)
  
  # Priors for individual-level parameters (intercept & slope)
  for( i in 1 : n_animals ) {
    beta0[i] ~ dnorm( mu_beta0, sd = sd_beta0 )
    beta1[i] ~ dnorm( mu_beta1, sd = sd_beta1 )
  }
  
  for( i in 1:n_used ) {
    used[i] ~ dcat( p[ start[i]:end[i]] )
    sum_ep[i] <- sum( ep[ start[i]:end[i] ] )
  }
  
  for( i in 1:n_all){
    ep[i] <- exp( beta0[individual[i]] + beta1[individual[i]] * prox_fish[i] )
    p[i] <- ep[i] / sum_ep[ ind_set[i] ]
  }
  
})

inits <- function(){
  list(
    mu_beta0 = rnorm(1, 0, 1), 
    mu_beta1 = rnorm(1, 0, 1), 
    sd_beta0 = rexp(1, 1), 
    sd_beta1 = rexp(1, 1), 
    beta0 = rnorm(constants$n_animals, 0, 1), 
    beta1 = rnorm(constants$n_animals, 0, 1)
  )}

params <- c( "beta0", "beta1", "mu_beta0", "mu_beta1", "sd_beta0", "sd_beta1" )

library(nimble)

begin <- Sys.time()
out <- nimble::nimbleMCMC(
  code = code, 
  constants = constants, 
  data = data, 
  inits = inits(), 
  monitors = params, 
  thin = 1, 
  niter = 1000, 
  nburnin = 500, 
  nchains = 1
)
done <- Sys.time()
done - begin

MCMCvis::MCMCsummary( out )
