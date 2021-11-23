################################################################################
# TITLE: Individual Missouri elk resource selection function analysis: Step 6 -
#   code for model checking
# PURPOSE: Checking the convergence of the models, validating the models, and 
#   comparing models, creating output for determining which models need to be
#   run longer, which are bogus (e.g., if there is no RX fire in any of used or
#   available, shouldn't even be in model)
# AUTHOR: Kyle Redilla, RECaP Lab
# CREATED: 2017-01-11
# LAST UPDATED ON 2017-01-11
################################################################################
# OPEN LIBRARIES 
library(rstan)
library(dplyr)
library(ggmcmc)
library(loo)
library(bayesplot)
library(fmsb)

begin <- Sys.time()
set.seed(011392)

# specify directories
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
resultdir <- paste(homedir, "Results/", sep = "")
datestr <- format(Sys.time(), "%Y-%m-%d")

# Number of individuals to process
nInd <- 11

# Specify custom functions 
# Radar chart customization function 
# From: https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 1,
                                        caxislabels = NULL, title = NULL, ...){
  fmsb::radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


################################################################################
# SECTION 1: POSTERIOR PREDICTIVE CHECK - Individual All Combos Models 
################################################################################

# Set result  directory for individual all combination models
comboresults <- paste(resultdir, "SSL_IndlAllCombos_2021-11-16/", sep = "")
# set it
setwd(comboresults)
# Get list of all files in results directory
files <- list.files()

# Load in rs_data
rs_data <- readRDS("../../Data_Processed/Telemetry/UsedAndAvail_WeeklyKDE_20211104.rds")

# Covariate names (in order)
covar_names <- c("bathymetry", "dist_land", "dist_500m", "slope", "sst", "wind", "logship", 
                 "logfish", "prox_fish_km", "prox_ship_km")
covar_labels <- c("Bathymetry", "Distance to\nland", "Distance to\nshelf", 
                  "Slope ", "SST", "Wind speed", "Shipping Intensity", "Fishing\nIntensity", 
                  "Distance to\nfishing", "Distance to\nshipping")

megasumms <- list()
megaloos <- list()
megawaics <- list()

# Check Rhat values & neff
for(i in 1:nInd){
  # load stan model fit
  load(files[grepl(paste0("Ind", i, "_"), files)])
  
  # Extract summary info 
  summs <- lapply(fitlst, "[[", 1)
  megasumms[[i]] <- summs
  
  loos <- lapply(fitlst, "[[", 2)
  megaloos[[i]] <- loos
  
  waics <- lapply(fitlst, "[[", 3)
  megawaics[[i]] <- waics
  
  # Calculate Rhat values > 1.1
  rhats <- lapply(summs, function(x) x %>% dplyr::select(Rhat)) %>% bind_rows(., .id="column_label")
  print(paste0("SSL ", i, " has ", length(which(rhats$Rhat > 1.1)), " Rhat values greater than 1.1"))
  # Calculate effective sample sizes < 100
  neffs <- lapply(summs, function(x) x %>% dplyr::select(n_eff)) %>% bind_rows(., .id="column_label")
  print(paste0("SSL ", i, " has ", length(which(neffs$n_eff < 100)), " Neff values less than 100"))
}

############################### 
######### RADAR PLOTS ######### 
############################### 

for(i in 1:nInd){
  
  mods <- read.csv(paste0(comboresults, "ModelVariableList_", i, ".csv")) %>% select(-X)
  # Adjust column names 
  colnames(mods) <- covar_names
  
  looranks <- as.data.frame(loo_compare(megaloos[[i]])) %>% tibble::rownames_to_column()
  looranks$modnum <- substr(looranks$rowname, 6,10)
  
  M <- length(mods$bathymetry)
  
  # what is the 95% mark
  top95 <- M - floor(M*0.95)
  # order and keep
  keep <- looranks[1:top95,]
  # Keep only the top 5% of models
  topmods <- mods[keep$modnum,]
  # Calculate percentage of top models containing each covariate
  pctcovars <- data.frame(t(colSums(topmods)/top95))
  
  pctcovars <- rbind(pctcovars, rep(0, 10))
  pctcovars <- rbind(pctcovars, rep(1, 10))
  rownames(pctcovars) <- c(i, "Min", "Max")
  pctcovars <- pctcovars[c("Max", "Min", i),]
  
  filename <- paste0(comboresults, "RadarPlot_",i,".png")
  png(filename = filename)
  create_beautiful_radarchart(pctcovars, vlabels=covar_labels, caxislabels = c(0, 0.25, 0.50, 0.75, 1))
  dev.off()
}






# Extract loo values and compare
looranks <- lapply(megaloos, function(x){as.data.frame(loo_compare(x)) %>% tibble::rownames_to_column()})



# Calculate probabilities? 
logOdds <- as.data.frame(summary(fitlst[[1]], pars = "beta", probs = c(0.025, 0.975))$summary)$mean
logitToProb <- function(lo) exp(lo)/(1+exp(lo))
groupProbs <- sapply(1:length(logOdds), function (i) round(logitToProb(logOdds[i])*100,1))
groupProbs

# LOO-PIT plot
## https://github.com/jgabry/bayes-vis-paper/blob/master/bayes-vis.R
y <- as.numeric(as.character(rs_data$used[which(rs_data$used == 1 & rs_data$ind_id == ind)]))
yrep1 <- as.matrix(fitlst[[1]], pars="log_lik")
loopsis <- lapply(fitlst, function(x)loo(x, save_psis=TRUE))
color_scheme_set("blue")
ppc_loo_pit_overlay(y, yrep1, lw = weights(loopsis[[1]]$psis_object)) + legend_none()
ggsave(filename = "plots/ppc_loo_pit_overlay1-new.png", width = 4.5, height = 3.75)



# ppc_dens_overlay --- NOT CORRECT
# y <- as.numeric(as.character(rs_data$used[which(rs_data$used == 1 & rs_data$ind_id == 1)]))
# yrep1 <- as.matrix(fitlst[[1]], pars="log_lik")
# 
# samp100 <- sample(nrow(yrep1), 100)
# 
# # overlaid densities
# color_scheme_set("blue")
# ppc_dens_overlay(y, yrep1[samp100, ]) + 
#   coord_cartesian(ylim = c(0, 0.7), xlim = c(0, 6)) +
#   legend_none()
# ggsave(filename = "plots/ppc_dens1.png", width = 4.5, height = 3.75)

# Extract WAIC values and compare

waicranks <- as.data.frame(loo_compare(waics))


# Simulate data? 
logOdds <- as.data.frame(summary(fitlst[[1]], pars = "beta", probs = c(0.025, 0.975))$summary)$mean

modparams <- read.csv(paste0("ModelVariableList_", ind, ".csv"))
################################################################################
# SECTION 2: POSTERIOR PREDICTIVE CHECK - Individual Global Models 
################################################################################

# Figure 1: plot results of Posterior 
# extract observed discrepancies
for(i in 1:11){
obs.disc <- rstan::extract(fitlst[[i]], 'chis_obs', F)[, , 1]
# extract simulated discrepancies
sim.disc <- rstan::extract(fitlst[[i]], 'chis_sim', F)[, , 1]
# number of post-warmup draws
pdraw <- dim(obs.disc)[1] * dim(obs.disc)[2]
# Calculate Bayesian P-value
pval <- sum(rstan::extract(fitlst[[i]], 'chis_sim', F)[, , 1] >
            rstan::extract(fitlst[[i]], 'chis_obs', F)[, , 1]) / pdraw
print(pval)
}

# initiate figure
# figure 1 name
fig1name <- paste("ObsVsSimDiscrepancy", datestr, ".pdf", sep = "")
pdf(fig1name, family = "Times", width = 6, height = 6)
# pick limits 
hi <- 2 * (range(obs.disc)[2]-range(obs.disc)[1]) + max(obs.disc, sim.disc)
lo <- min(obs.disc, sim.disc) - (range(obs.disc)[2]-range(obs.disc)[1])/10
plot(obs.disc, sim.disc, 
     xlab = 'Observed discrepancy', 
     ylab = 'Simulated discrepancy',
     xlim = c(lo, hi), ylim = c(lo, hi),
     type = "n")
points(obs.disc, sim.disc, pch = 19, cex = 0.3)
# add 1 to 1 line
lines(c(0, 800000), c(0, 800000))
# 
text(hi - (0.6 * (hi - lo)), 
     hi - (0.3 * (hi - lo)), 
     labels = paste("p-value = ", 
                    as.character(round(pval, digits = 3)), sep = ""))
dev.off()




P_labels <- data.frame(
  Parameter=c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", 
              "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]"),
  Label=c("Bathymetry (m)", "Distance to land (m)", "Distance to shelf break (m)", 
          "Slope (degrees)", "Avg. sea surface\ntemperature (C)", 
          "Avg. wind speed (m/s)", "Log(Shipping traffic (km))", "Log(Fishing traffic (km))", 
          "Distance to fishing (km)", "Distance to shipping (km)"))
M_labels = c("774PWS", "775PWS", "776PWS", "777PWS", 
             "781KOD", "782KOD", "783KOD", "784KOD", "785KOD", "786KOD", "788KOD")

# figure 2 - diagnostic plot report
# Using package ggmcmc to plot variety of diagnostics plots
for(i in 1:11){
  # create coda mcmc list
  s <- As.mcmc.list(fitlst[[i]], pars = c("beta")) 
  # prepare for plotting functions
  S <- ggs(s, par_labels=P_labels) 
  # number of parameters monitored
  pK <- length(levels(S$Parameter))
  # initialize traceplots figure
  fig2name <- paste("Traceplot_",M_labels[i],"_", datestr, ".pdf", sep = "")
  ggmcmc(S, file = fig2name, param_page = 6, 
         plot = c("traceplot", "autocorrelation", "density", 
                  "running", "caterpillar"))
  
  # figure 3 - Potential scale reduction factor and Geweke diagnostics
  fig3name <- paste("ScaleAndGeweke_",M_labels[i],"_", datestr, ".pdf", sep = "")
  # how many "pages worth" at 
  # pgs <- pK %/% 25
  pdf(fig3name)
  ggs_Rhat(S) + xlab("R_hat")
  ggs_geweke(S)
  dev.off()
}

# Isolate out just caterpillar plots for presentation
for(i in 1:11){
  # create coda mcmc list
  s <- As.mcmc.list(fitlst[[i]], pars = c("beta")) 
  # prepare for plotting functions
  S <- ggs(s, par_labels=P_labels) 
  ggs_caterpillar(S, model_labels=M_labels[i]) +
    theme(axis.text.x = element_text(size=25), 
          axis.text.y = element_text(size=25)) +
    ylab("")
  ggsave(paste0(resultdir, "SSL_IndlGlobalFixedEffects_2021-11-09/CaterpillarPlot_", M_labels[i], ".png"), width=8, height=9, units="in")
}

# Leave one out cross validation (Loo-CV)
# Code from: http://blackwell.math.yorku.ca/MATH6635/files/Stan_first_examples.html#step-5-check-whether-hmc-worked--
fitlst  %>% 
  lapply(function(fit) {
    log_lik1 <- extract_log_lik(fit, merge_chains = FALSE)
    rel_n_eff <- relative_eff(exp(-log_lik1))
    loo(log_lik1, r_eff = rel_n_eff, cores = 4)
  }) -> loolist

loo_compare(loolist)

# END OF SECTION 1
################################################################################

## PARKING LOT

## Bayesian p-value for each "top model" 
#   load step 7 results (top models for individuals)
step7.results.path <- paste(resultdir, "ERS1_step7_results_2017-04-09.rda", 
                            sep = '')
load(step7.results.path)
# fitlst is list object containing fits with chi-square simulated at each
#   data point
for(i in 1:88){
  # extract observed discrepancies
  obs.disc <- rstan::extract(fitlst[[i]], 'chis_obs', F)[, , 1]
  # extract simulated discrepancies
  sim.disc <- rstan::extract(fitlst[[i]], 'chis_sim', F)[, , 1]
  # number of post-warmup draws
  pdraw <- dim(obs.disc)[1] * dim(obs.disc)[2]
  # Calculate Bayesian P-value
  pval <- sum(rstan::extract(fitlst[[i]], 'chis_sim', F)[, , 1] >
                rstan::extract(fitlst[[i]], 'chis_obs', F)[, , 1]) / pdraw
  
  # initiate figure
  # figure 1 name
  figname <- paste("ERS1_step6_fig1_elk", i, "_2017-01-11.pdf", sep = "")
  # output directory
  outdir <- paste(resultdir, "ERS1_step6_figure1_top_models/", sep = '')
  # output path for elk i
  outpath <- paste(outdir, figname, sep = '')
  # initialize figure
  pdf(outpath, family = "Times", width = 6, height = 6)
  # pick limits 
  hi <- 2 * (range(obs.disc)[2]-range(obs.disc)[1]) + max(obs.disc, sim.disc)
  lo <- min(obs.disc, sim.disc) - (range(obs.disc)[2]-range(obs.disc)[1])/10
  # paste string for identifying title
  title.str <- paste("Elk Number: ", i, sep = '')
  plot(obs.disc, sim.disc, 
       xlab = 'Observed discrepancy', 
       ylab = 'Simulated discrepancy',
       main = title.str,
       xlim = c(lo, hi), ylim = c(lo, hi),
       type = "n")
  points(obs.disc, sim.disc, pch = 19, cex = 0.3)
  # add 1 to 1 line
  lines(c(0, 800000), c(0, 800000))
  # 
  text(hi - (0.6 * (hi - lo)), 
       hi - (0.3 * (hi - lo)), 
       labels = paste("p-value = ", 
                      as.character(round(pval, digits = 3)), sep = ""))
  dev.off()
}
