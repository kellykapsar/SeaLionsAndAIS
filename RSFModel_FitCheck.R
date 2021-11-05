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
begin <- Sys.time()

# specify directories
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
resultdir <- paste(homedir, "Results/", sep = "")
datestr <- format(Sys.time(), "%Y-%m-%d")
# now change over to step 6 result dir
step6.resultdir <- paste("SSL_ModelCheck_Results_", 
                         datestr,"/", sep = "")
# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(file.path(resultdir, step6.resultdir)), 
       dir.create(file.path(resultdir, step6.resultdir)), FALSE)

################################################################################
# SECTION 1: POSTERIOR PREDICTIVE CHECK STEP 3 # KK - CHANGED TO SEPT 4... 
################################################################################

# result directory for population-level analysis
step4.resultdir <- paste(resultdir, "SSL_RSF_Results_2021-11-01/", sep = "")
# set it
setwd(step4.resultdir)
# load stan model fit
load("SSL_RSF_results_2021-11-01.rda")

# result directory for step 6
step6.resultdir <- paste(resultdir, step6.resultdir, sep = "")
# set it
setwd(step6.resultdir)

# Figure 1: plot results of Posterior 
# extract observed discrepancies
obs.disc <- rstan::extract(fitlst[[1]], 'chis_obs', F)[, , 1]
# extract simulated discrepancies
sim.disc <- rstan::extract(fitlst[[1]], 'chis_sim', F)[, , 1]
# number of post-warmup draws
pdraw <- dim(obs.disc)[1] * dim(obs.disc)[2]
# Calculate Bayesian P-value
pval <- sum(rstan::extract(fitlst[[1]], 'chis_sim', F)[, , 1] >
            rstan::extract(fitlst[[1]], 'chis_obs', F)[, , 1]) / pdraw

# initiate figure
# figure 1 name
fig1name <- paste("SSL_RSF_results_Figure1_", datestr, ".pdf", sep = "")
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


# figure 2 - diagnostic plot report
# Using package ggmcmc to plot variety of diagnostics plots
# create coda mcmc list
s <- As.mcmc.list(fit1, pars = c("mu", "stdev", "lp__")) 
# prepare for plotting functions
S <- ggs(s) 
# number of parameters monitored
pK <- length(levels(S$Parameter))
# initialize traceplots figure
fig2name <- paste("ERS1_step6_figure2_", datestr, ".pdf", sep = "")
ggmcmc(S, file = fig2name, param_page = 6, 
       plot = c("traceplot", "autocorrelation", "density", 
                "running", "caterpillar"))


# figure 3 - Potential scale reduction factor and Geweke diagnostics
fig3name <- paste("ERS1_step6_figure3_", datestr, ".pdf", sep = "")
# how many "pages worth" at 
pgs <- pK %/% 25
pdf(fig3name, height = 9 * pgs)
ggs_Rhat(S) + xlab("R_hat")
ggs_geweke(S)
dev.off()

# messing
pgs <- pK %/% 6

f1 <- ggs_traceplot(S) + 
  theme_fivethirtyeight()
f2 <- ggs_density(S) + 
  theme_solarized(light=TRUE)
pdf("test.pdf", height = 9 * pgs)
grid.arrange(f1, f2, ncol=2, nrow=1)
dev.off()

# END OF SECTION 1
################################################################################


