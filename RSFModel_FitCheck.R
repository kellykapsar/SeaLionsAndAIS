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
library(RColorBrewer)

begin <- Sys.time()
set.seed(011392)

# specify directories
homedir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/" # kk - Don't know what the difference between workdir and homedir is
resultdir <- paste(homedir, "Results/", sep = "")
datestr <- format(Sys.time(), "%Y-%m-%d")

# Number of individuals to process
nInd <- 11

#################################### 
######### IMPORT FUNCTIONS ######### 
####################################


# Specify custom functions 
# Radar chart customization function 
# From: https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data),
                                        caxislabels = NULL, title = NULL, ...){
  fmsb::radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, plwd = 3, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlabels = vlabels,
    cex.main = 2.5,
    caxislabels = caxislabels, 
    title = title, 
    calcex=2, palcex=2, vlcex=2,...
  )
}


# Multiplot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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

M_labels = c("774PWS", "775PWS", "776PWS", "777PWS", 
             "781KOD", "782KOD", "783KOD", "784KOD", "785KOD", "786KOD", "788KOD")

# Get colors for individuals
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

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
topmodnums <- list()

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
  
  # Calculate top model numbers for each individual
  topmod <- as.data.frame(loo_compare(megaloos[[i]])) %>% tibble::rownames_to_column()
  topmod$modnum <- as.numeric(substr(topmod$rowname, 6,10))
  topmodnums[[i]] <- topmod$modnum[1]
  
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

############# Big central radar plot with labels

mods <- read.csv(paste0(comboresults, "ModelVariableList_1.csv")) %>% select(-X)
# Adjust column names 
colnames(mods) <- covar_names

looranks <- as.data.frame(loo_compare(megaloos[[1]])) %>% tibble::rownames_to_column()
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
filename <- paste0(comboresults, "RadarPlot_1BIG.png")
png(filename = filename, width = 15, height=15, units="in", res=200)
create_beautiful_radarchart(pctcovars, vlabels=covar_labels, 
                            caxislabels = c(0, 0.25, 0.50, 0.75, 1), 
                            title=M_labels[[1]], color=getPalette(11)[1])
dev.off()


pctcovarsall <- data.frame()
# All other covar plots without labels
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
  
  # Add to overall data frame 
  pctcovarsall <- rbind(pctcovarsall, pctcovars)
  
  # Prep min and max values for radar 
  pctcovars <- rbind(pctcovars, rep(0, 10))
  pctcovars <- rbind(pctcovars, rep(1, 10))
  rownames(pctcovars) <- c(i, "Min", "Max")
  pctcovars <- pctcovars[c("Max", "Min", i),]
  
  filename <- paste0(comboresults, "RadarPlot_",i,".png")
  png(filename = filename)
  create_beautiful_radarchart(pctcovars, vlabels="", 
                              caxislabels = c(0, 0.25, 0.50, 0.75, 1), 
                              title=M_labels[[i]], color=getPalette(11)[i])
  dev.off()
}

pctcovarsall <- round(pctcovarsall, 2)

############################################ 
######### Chi-square on Top Models ######### 
############################################

# model path
comp.modpath <- paste(resultdir, "SSL_IndlGlobalFixedEffects_2021-11-09/model.rda", sep = "")

# read in compiled model object
mod <- readRDS(comp.modpath)

topfitlst <- list()

for(i in 1:nInd){
  
  print(paste0("Processing individual: ", i))
  
  # Import variable combinations for this individual 
  mods <- read.csv(paste0(comboresults, "ModelVariableList_", i, ".csv")) %>% select(-X)
  
  # Adjust column names 
  colnames(mods) <- covar_names
  
  # Identify covariates included in the top model 
  betatruefalse <- mods[topmodnums[[i]],]
  betanums[[i]] <- which(betatruefalse[1,] == TRUE)
  
  # Calculate total number of covars in top model
  K <- sum(betatruefalse[1,])
  
  # Specify number of choices per choice set
  C <- 6 
  
  # Isolate out individual data 
  rs_data_subset <- rs_data[rs_data$ind_id == i, ]
  
  # Specify number of choice sets for this individual
  N <- length(unique(rs_data_subset$choice_id))
  
  # create design array with all covariates 
  x <- rsf_array(rs_data, c(N, C, 10)) 
  
  # Remove unused covariates from data array
  x.temp <- x[,,betanums[[i]]]

  # must enter data into a list
  data <- list(
    C = C, K = K, N = N, 
    x = x.temp,
    y = rep(1, N),
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
  params <- c('beta', 'chis_obs', 'chis_sim')
  
  fit <- sampling(mod, data = data, pars = params, init = inits,
                  chains =4, iter = 1000, warmup = 200, thin = 1)
  
  topfitlst <- c(topfitlst, fit)
  
}

saveRDS(topfitlst, paste0(resultdir, "SSL_IndlAllCombos_2021-11-16/TopModelFits_ChiSquare.rds"))
topfitlst <- readRDS(paste0(resultdir, "SSL_IndlAllCombos_2021-11-16/TopModelFits_ChiSquare.rds"))

# Figure 1: plot results of Posterior 
# extract observed discrepancies
for(i in 1:nInd){
  obs.disc <- rstan::extract(topfitlst[[i]], 'chis_obs', F)[, , 1]
  # extract simulated discrepancies
  sim.disc <- rstan::extract(topfitlst[[i]], 'chis_sim', F)[, , 1]
  # number of post-warmup draws
  pdraw <- dim(obs.disc)[1] * dim(obs.disc)[2]
  # Calculate Bayesian P-value
  pval <- sum(rstan::extract(topfitlst[[i]], 'chis_sim', F)[, , 1] >
                rstan::extract(topfitlst[[i]], 'chis_obs', F)[, , 1]) / pdraw
  print(pval)
}


############################################ 
######### Caterpillar plots for Top Models ######### 
############################################


covar_labels <- data.frame(
  Name = covar_names, 
  Label=c("Bathymetry (m)", "Distance to land (m)", "Distance to shelf break (m)", 
          "Slope (degrees)", "Avg. sea surface\ntemperature (C)", 
          "Avg. wind speed (m/s)", "Log(Shipping traffic (km))", "Log(Fishing traffic (km))", 
          "Distance to fishing (km)", "Distance to shipping (km)"))

Parameters <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", 
                            "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]")


P_labellist <- list() # Names of covariates in each top model (in order)
betanums <- list() # Order number for covariates in each top model (1-10)

# Identify covariates in each  top model 
for(i in 1:nInd){
  
  # Import variable combinations for this individual 
  mods <- read.csv(paste0(comboresults, "ModelVariableList_", i, ".csv")) %>% select(-X)
  
  # Adjust column names 
  colnames(mods) <- covar_names
  
  # Identify covariates included in the top model 
  betatruefalse <- mods[topmodnums[[i]],]
  betanums[[i]] <- which(betatruefalse[1,] == TRUE)
  
  covar_l <- covar_labels[betanums[[i]],]
  covar_l$Parameter <- Parameters[1:length(betanums[[i]])]
  
  P_labellist[[i]] <- covar_l
}


# Isolate out just caterpillar plots for presentation
for(i in 1:11){
  # create coda mcmc list
  s <- As.mcmc.list(topfitlst[[i]], pars = c("beta")) 
  # prepare for plotting functions
  S <- ggs(s, par_labels=P_labellist[[i]]) 
  ggs_caterpillar(S, model_labels=M_labels[i]) +
    theme(axis.text.x = element_text(size=25), 
          axis.text.y = element_text(size=25)) +
    ylab("")
  ggsave(paste0(resultdir, "SSL_IndlAllCombos_2021-11-16/CaterpillarPlot_", M_labels[i], ".png"), width=8, height=9, units="in")
}

# Other miscellaneous plots 
fitlstmcmc <- lapply(topfitlst, As.mcmc.list, pars=c("beta"))

mcmc_areas(fitlstmcmc[[1]], prob=0.95)
mcmc_scatter(fitlstmcmc[[1]], pars=c("beta[2]", "beta[8]"))

rethinking::precis(topfitlst[[1]])

####################
# Significance counts for variables in best fitting models

topmodbetas <- data.frame()

for(i in 1:11){
temp <- as.data.frame(rstan::summary(topfitlst[[i]], pars=c("beta"))$summary)
temp$Parameter <- rownames(temp)
temp$ind_id <- i
temp <- left_join(P_labellist[[i]], temp)
topmodbetas <- rbind(topmodbetas, temp)
}

topmodbetas$sigpos <- ifelse(topmodbetas$`2.5%` > 0 & topmodbetas$`97.5%` > 0, 1, 0)
topmodbetas$signeg <- ifelse(topmodbetas$`2.5%` < 0 & topmodbetas$`97.5%` < 0, 1, 0)
topmodbetas$signonsig <- ifelse(topmodbetas$`2.5%` < 0 & topmodbetas$`97.5%` > 0, 1, 0)


test <- topmodbetas %>% group_by(Label) %>% select(Label, sigpos, signeg, signonsig) %>% 
          gather(key="Significance", value="Count", -Label) %>% ungroup() %>% 
          group_by(Label, Significance) %>% summarize(Count=sum(Count))

## set the levels in order we want


p1 <- ggplot(test, aes(fill=factor(Significance, levels=c("sigpos", "signonsig", "signeg")), y=Count, x= reorder(Label, -Count))) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "PRGn", labels=c("Sig. Positive", "Non-significant","Sig. Negative")) +
  # theme(legend.text = element_text(c("Sig. Positive", "Non-significant","Sig. Negative"))) +
  labs(fill='Effect') +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10, 12)) +
  ylab("Number of Models") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust=1), text = element_text(size = 20))

ggsave(plot=p1, filename=paste0(resultdir, "SSL_IndlAllCombos_2021-11-16/SignificanceCountsBarPlot.png"), 
       width=8, height=9, units="in")


#################################################################
######### CATERPILLAR PLOT FOR ALL INDIVIDUALS COMBINED #########
################### EXAMPLE CODE NOT MODIFIED ###################
#################################################################

### figure 2: Population & individual global model point estimates 

## Create figure 2 with credible intervals on individual estimates
# figure 2 filepath
fig1path <- paste(imagedir, "ERS1_step7_figure1_",
                  datestr, ".jpg", sep = "")
# initialize ggplot
pop.p <- ggplot(data = mean.df_long, aes(x = variable, y = mean)) +
  # add horizontal line for zero
  geom_hline(yintercept = 0, color = "navy", lwd = 1.5) + 
  geom_pointrange(aes(ymin = CI_L, ymax = CI_H, x = variable, y = mean),
                  size = 0.2, position = position_jitter(0.3),
                  alpha = 0.3, color = "darkcyan") + 
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black")) + 
  ylab('Mean coefficient estimate') + 
  xlab('Resource covariate') + 
  scale_y_continuous(limits = c(-6, 6)) + 
  geom_crossbar(data = pop.results, 
                aes(ymin=pop.results$'5%.m' - pop.results$'95%.s', 
                    ymax=pop.results$'95%.m' + pop.results$'95%.s',
                    x = vars, y = mean.m), color = "gray50",
                width = 0.7, size = 1, fatten = 0.5) + 
  geom_errorbar(data = pop.results, 
                aes(ymin=mean.m-mean.s, ymax=mean.m+mean.s,
                    x = vars, y = mean.m), color = "gray75",
                width = 0.7, size = 1) 
# Initialize plotting device
jpeg(fig1path, width = 12, height = 7, units = "in", res = 300)
# plot the figure and close
pop.p
dev.off()


######################################################################################

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
