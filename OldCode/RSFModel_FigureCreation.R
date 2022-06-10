################################################################################
# TITLE: Individual elk RSF analysis: Step 7 - Code for conducting MV analyses
#   and generating figures
# PURPOSE: This code is for multivariate analysis on the discerete choice
#   model coefficients, and generating figures from results of 
#   population-level global RSF fits and individual-level global 
#   RSF fits, and the respective comparison/contrasting between the two
# AUTHOR: Kyle Redilla, RECaP Lab
# CREATED: 2016-12-16 
# LAST UPDATED ON: 2017-04-19
################################################################################
# OPEN LIBRARIES 
library(rstan)
library(ggplot2)
library(dplyr)
# library(MASS)
# library(reshape2)
# library(cluster)
# library(vegan)
# library(rgdal)
library(RColorBrewer)
# library(raster)
# library(viridis)
# library(rasterVis)
# library(qdap)
library(foreach)
# library(gridExtra)
library(fmsb)
library(loo)
################################################################################
# Start with clean environment
rm(list = ls())
# Create relevant working directories
datestr <- format(Sys.time(), "%Y-%m-%d")
# Work directory
workdir <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/"
# Data directory
datadir <- paste(workdir, "Data_Processed/", sep = "")
# Results directory (store figures)
resultdir <- paste(workdir, "Results/", sep = "")
# output directory
outdir <- resultdir
# imagedir is outdir now
imagedir <- outdir
# individual model selection results path
step5dir <- paste(resultdir, "SSL_IndlAllCombos_2021-11-16/", sep = '')

# scaled raster data directory
# rasdatadir <- paste(datadir, "scaled_rasters/", sep = '')
# erz layer directory
# erzlayerdir <- paste(datadir, "created_layers/", sep = '')
# Start timing the analysis.
print(Timing <- data.frame(Start = Sys.time(), End = NA, Duration = NA))

# Covariate names (in order)
covar_names <- c("bathymetry", "dist_land", "dist_500m", "slope", "sst", "wind", "logship", 
                 "logfish", "prox_fish_km", "prox_ship_km")
covar_labels <- c("Bathymetry", "Distance to\nland", "Distance to\nshelf", 
                  "Slope ", "SST", "Wind speed", "Shipping Intensity", "Fishing\nIntensity", 
                  "Distance to\nfishing", "Distance to\nshipping")

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
# SECTION 1: COEFFICIENT ESITMATES CATERPILLAR PLOTS
################################################################################
# Specify the number of individuals with data
nInd <- 11


### Read in data

## Data 
# read in uncopmressed dataset for elk rs
# resource selection data path
all.rs.data.path <-  paste(datadir, "Telemetry/UsedAndAvail_WeeklyKDE_20211104.rds", sep = "")
all_rs_data <- readRDS(all.rs.data.path)

# Read in capture individual info
# capture.data.path <- paste(datadir, "elk_age_capture_info.csv", sep = "")
# capture_data <- read.csv(capture.data.path)

# read in loc data
# loc.data.path <- paste(datadir, "all_elk_locations.csv", sep = "")
# loc_data <- read.csv(loc.data.path)

## Global model results
# Load individual model results
# indiv.stan.results.path <- paste(resultdir, 
#                                  "SSL_IndlGlobalFixedEffects_2021-11-05/SSL_IndlGlobalRSF_Results2021-11-05.rda", sep = "")
# load(indiv.stan.results.path)
# load population model results
# pop.stan.results.path <- paste(resultdir, "ERS1_step2_results_2017-02-16.rda", 
#                                sep = "")
# load(pop.stan.results.path)

## Model selection results
setwd(step5dir)
# Get list of all files in results directory
files <- list.files()

# loop through 
for(i in 1:nInd){
  # load stan model fits
  load(files[grepl(paste0("Ind", i, "_"), files)])
  # Models 
  mods <- read.csv(paste0(step5dir, "ModelVariableList_", i, ".csv")) %>% select(-X)
  # Adjust column names 
  colnames(mods) <- covar_names
  
  # extract model data frame from first part of each list
  mod.loos <- lapply(fitlst, '[[', 2)
  looranks <- as.data.frame(loo_compare(mod.loos)) %>% tibble::rownames_to_column()
  looranks$modnum <- substr(looranks$rowname, 6,10)
  
  M <- length(mod.loos)

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
  
  filename <- paste0(step5dir, "RadarPlot_",i,".png")
  png(filename = filename)
  create_beautiful_radarchart(pctcovars, vlabels=covar_labels, caxislabels = c(0, 0.25, 0.50, 0.75, 1))
  dev.off()
}


############################################################################




top.mods <- foreach(i = 1:88) %do% {
  # filepath
  path <- paste(step5dir, 'job4_results_elk_', i, '.rdata', sep = '')
  # load it
  load(path)
  # extract model data frame from first part of each list
  mod.df.lst <- lapply(fitlst, '[[', 1)
  # Remove timer element
  mod.df.lst[length(mod.df.lst)] <- NULL
  # number of models in model set (minus 1 for time calculation)
  M <- length(mod.df.lst)
  # unlist the model df list and populate matrix
  mod.comp.mat <- matrix(unlist(mod.df.lst), 
                         nrow = M, byrow = TRUE)
  # change to numeric
  mod.comp.mat <- apply(mod.comp.mat, 2, as.numeric)
  # discard all rows containing dist_2_gravel parameter
  mod.comp.mat <- mod.comp.mat[mod.comp.mat[, 5] == 1,]
  mod.comp.mat <- mod.comp.mat[,-5]
  # order and keep
  keep <- mod.comp.mat[order(mod.comp.mat[, 21]),][1,]
  # turn all '1' to 'FALSE'
  keep[keep == 1] = 0
  # turn all NA to TRUE
  keep[is.na(keep)] = TRUE
  keep
}



### Prepare results/data for figures

## Assemble individual information 
require(rstan)
require(reshape2)
require(dplyr)
# variable names (distance to gravel removed)
vars <- c("Distance to 2-track", "Distance to paved road", "Slope", 
          "Aspect", "Road Density", 
          "Canopy cover", "Distance to edge", "IJI", 
          "Years since Rx burn", "Forest", "Woodland", 
          "Savanna", "Shrubland", "Glade", "Warm-season Grassland", 
          "Cool-season grassland", "Forage opening")
# make columns for locations in 2011, 2012, 2013, 2014
all_rs_data$loc.2011 <- 0; all_rs_data$loc.2011[all_rs_data$year == 1] = 1
all_rs_data$loc.2012 <- 0; all_rs_data$loc.2012[all_rs_data$year == 2] = 1
all_rs_data$loc.2013 <- 0; all_rs_data$loc.2013[all_rs_data$year == 3] = 1
all_rs_data$loc.2014 <- 0; all_rs_data$loc.2014[all_rs_data$year == 4] = 1
##  Make separate table for intrinsic factors (sex, cohort, age)
# animal sex, id, and cohort
indiv_info <- as.data.frame(all_rs_data) %>% group_by(ind_id) %>% 
  filter(Use == 1) %>% summarise(
    cohort = unique(Cohort), sex = unique(female), count = sum(Use),
    count.2011 = sum(loc.2011), count.2012 = sum(loc.2012),
    count.2013 = sum(loc.2013), count.2014 = sum(loc.2014),
    Elk_ID = unique(Elk_ID)
  )
# Add string for sex
indiv_info$sex[indiv_info$sex == 1] = "Female"
indiv_info$sex[indiv_info$sex == -1] = "Male"
# Add capture info (merge)
indiv_info <- merge(indiv_info, capture_data, by = 'Elk_ID')
# Remove unnecessary columns
indiv_info[,c('Cohort', 'Sex')] <- c(NULL, NULL)
# Modify age class: if released as juvenile or adult: adult
# all else: subadult
indiv_info$Age.eff <- factor(levels = c("Sub-adult", "Adult"),
                             rep("Adult", 88))
indiv_info$Age.eff[indiv_info$Capture.Age == "Yearling"] = "Sub-adult"

## Prepare population results as an array
# subset output of "stanfit" class method for summary() to get stats
pop.results <- summary(fit1)$summary
# Collect mean, sd, and limits of 95% CI on the mean and sd's (mus & sigmas)
pop.means <- pop.results[1:17, c("mean", "sd", "2.5%", "97.5%")]
pop.sds <- pop.results[18:34, c("mean", "sd", "2.5%", "97.5%")]
# Calculate 90% CI (5 and 95 percentiles) for mus
pop.means <- cbind(pop.means, 
                   qnorm(0.05, pop.means[,'mean'], pop.means[,'sd']))
pop.means <- cbind(pop.means, 
                   qnorm(0.95, pop.means[,'mean'], pop.means[,'sd']))
# Calculate 90% CI for stdevs
pop.sds <- cbind(pop.sds, 
                   qnorm(0.05, pop.sds[,'mean'], pop.sds[,'sd']))
pop.sds <- cbind(pop.sds, 
                   qnorm(0.95, pop.sds[,'mean'], pop.sds[,'sd']))
# rename 
colnames(pop.means)[5:6] <- colnames(pop.sds)[5:6] <- c("5%", "95%")
# combine point esimates of mu and sigma
pop.results <- cbind(pop.means, pop.sds)
pop.results <- as.data.frame(pop.results)
# rename columns based on whether applies to mu or sigma
names(pop.results) <- c("mean.m","sd.m", "2.5%.m", "97.5%.m", "5%.m", "95%.m",
                        "mean.s","sd.s", "2.5%.s", "97.5%.s", "5%.s", "95%.s")
# Add variable names as column
pop.results$vars <- vars
# make limits for pop params: mu +- sigma
pop.results$ymax <- pop.results$mean.m+pop.results$mean.s
pop.results$ymin <- pop.results$mean.m-pop.results$mean.s

## Preapre individual results as a matrix
# setup array to hold mean coefficient results: 88 individuals,
#   17 parameters, 7 stats
mean.arr <- array(dim = c(88, 17, 7))
# fill array of mean coefficient results + posterior intervals
for(i in 1:88){
  temp_smmry <- summary(fitlst[[i]])
  # Summary of coefficients
  temp_betas <- temp_smmry$summary[1:17, ]
  # populate with coeff means, se's, 2.5, 5, 95, 97.5
  mean.arr[i,, 1] <- temp_betas[,'mean']
  mean.arr[i,, 2] <- temp_betas[,'se_mean']
  mean.arr[i,, 3] <- temp_betas[,'sd']
  mean.arr[i,, 4] <- qnorm(0.0275, mean.arr[i,, 1], temp_betas[, 'sd'])
  mean.arr[i,, 5] <- qnorm(0.05, mean.arr[i,, 1], temp_betas[, 'sd'])
  mean.arr[i,, 6] <- qnorm(0.95, mean.arr[i,, 1], temp_betas[, 'sd'])
  mean.arr[i,, 7] <- qnorm(0.975, mean.arr[i,, 1], temp_betas[, 'sd'])
}
# extract means to a data frame
mean.df <- as.data.frame(mean.arr[,,1])
# extract UPPER limits of 95% mean CI to data frame
meanCIh.df <- as.data.frame(mean.arr[,,7])
# extract LOWER limits of 95% mean CI to data frame
meanCIl.df <- as.data.frame(mean.arr[,,4])
# Add env variable names
names(mean.df) <- vars
names(meanCIh.df) <- vars
names(meanCIl.df) <- vars      
# Add individual id column
mean.df$idvar <- rownames(mean.df)
meanCIh.df$idvar <- rownames(meanCIh.df)
meanCIl.df$idvar <- rownames(meanCIl.df)
# Timevar (for making the data frame long format)
mean.df$timevar <- rep(1, 88)
meanCIh.df$timevar <- rep(1, 88)
meanCIl.df$timevar <- rep(1, 88)
# Melt data frames to long format
mean.df_long <- melt(mean.df)
meanCIh.df_long <- melt(meanCIh.df)
meanCIl.df_long <- melt(meanCIl.df)
# remove timevar?
mean.df_long <- mean.df_long[mean.df_long$variable != "timevar",]
meanCIh.df_long <- meanCIh.df_long[meanCIh.df_long$variable != "timevar",]
meanCIl.df_long <- meanCIl.df_long[meanCIl.df_long$variable != "timevar",]
# merge data frames 
mean.df_long <- merge(mean.df_long, meanCIl.df_long,
                     by = c("idvar", "variable"))
mean.df_long <- merge(mean.df_long, meanCIh.df_long,
                      by = c("idvar", "variable"))
# change names of variables
names(mean.df_long) <- c("idvar", "variable", "mean", "CI_L", "CI_H")
# Now, mean.df_long is a data frame ready for use in ggplot2

## Miscellaneous labels
# Summary statistic labels
summ.labs <- c("Mean", "SE", "sigma", "2.5%", "5%", "95%", "97.5%")
# Name the dimensions of the array
dimnames(mean.arr) <- list(paste("elk", 1:88, sep = ''),
                           vars, summ.labs)
# continuous predictor names
con.vars <- vars[1:9]
# habitat predictor naems
hab.vars <- vars[10:17]

## Bind indiv info columns (Id, cohort, sex) to indiv results
# vector of identity variable names
idvec <- c("Elk_ID", "sex", "cohort", "count", "Capture.Age", 
           "Capture.Location.in.KY", "Age.eff")
# bind to means
mean.mat <- cbind(mean.arr[,, 1], indiv_info[, idvec])
# bind to std errors (probably irrelevant)
se.mat <- cbind(mean.arr[,, 2], indiv_info[, idvec])
# bind to sigmas
sd.mat <- cbind(mean.arr[,, 3], indiv_info[, idvec])
# Give the proper variable names for each set of parameter ests
colnames(mean.mat)[1:17] <- vars
colnames(se.mat)[1:17] <- vars
colnames(sd.mat)[1:17] <- vars



################################################################################

### comparison of population and individual estimates

## Collect output for result inquiries in .Rout file
# Creates result directory for this step on the specified date if not created
ifelse(!dir.exists(file.path(outdir)), 
       dir.create(file.path(outdir)), FALSE)
# file to output results
outpath <- paste(outdir, "output.Rout", sep = "")
# initiate sink
sink(outpath)
# Collect output for how many parameters were included in top models
cat("What are the sample statistics on the number of parameters
    included in top models?"); cat("\n")
# create vector of number of parameters in top models
k.count <- unlist(lapply(lapply(top.mods, '[', 1:17), sum))
# provide a summary
summary(k.count); cat("\n", "\n")
# How many times was the habitat variable included in the top model?
cat("How many times were certain variables included in the 
    top model?"); cat("\n")
# fill binary vector indicating inclusion of habitat
hab.count <- sum(unlist(lapply(top.mods, '[', 10)))
# fill binary vector indicating inclusion of slope
slp.count <- sum(unlist(lapply(top.mods, '[', 3)))
# fill binary vector indicating inclusion of iji
iji.count <- sum(unlist(lapply(top.mods, '[', 8)))
# fill binary vector indicating inclusion of canopy cover
cancov.count <- sum(unlist(lapply(top.mods, '[', 6)))
# fill binary vector indicating inclusion of distance to paved road
dprd.count <- sum(unlist(lapply(top.mods, '[', 2)))
# fill binary vector indicating inclusion of Years since Rx Burn
rxburn.count <- sum(unlist(lapply(top.mods, '[', 9)))
# fill binary vector indicating inclusion of road density
rdens.count <- sum(unlist(lapply(top.mods, '[', 5)))
# fill binary vector indicating inclusion of distance to two-track
d2t.count <- sum(unlist(lapply(top.mods, '[', 1)))
# fill binary vector indicating inclusion of distance to edge
d2e.count <- sum(unlist(lapply(top.mods, '[', 7)))
# fill binary vector indicating inclusion of Aspect
aspect.count <- sum(unlist(lapply(top.mods, '[', 4)))
cat(
  paste("Habitat was included in", hab.count, "of top models"),"\n",
  paste("Slope was included in", slp.count, "of top models"),"\n",
  paste("IJI was included in", iji.count, "of top models"),"\n",
  paste("Canopy cover was included in", cancov.count, "of top models"),"\n",
  paste("Distance to paved rds was included in", dprd.count, 
        "of top models"),"\n",
  paste("Years since prescribed burn was included in", rxburn.count, 
        "of top models"),"\n",
  paste("Road density was included in", rdens.count, "of top models"),"\n",
  paste("Distance to 2track was included in", d2t.count, "of top models"),"\n",
  paste("distance to wooded edge was included in", d2e.count, 
        "of top models"),"\n",
  paste("Aspect was included in", aspect.count, "of top models"),"\n","\n"
  
)



# does zero fall within credible interval for each param?
cat("Which mean and sd hyperparameter credible intervals 
      overlapped zero?")
cat("\n")
# create vector to rename wit vars
mean.creds <- (pop.means[, 3] < 0 & pop.means[, 4] > 0)
# rename
names(mean.creds) <- vars
# make data frame
cred.ints <- as.data.frame(mean.creds)
# Add mean point estimates
cred.ints$mean.ests <- round(pop.means[, 1], digits = 2)
# Create vector to rename with vars
# This is the lower bound of the 95% CI on the mean minus
#   the upper bound of the 95% CI on the sd 
sd.creds <- ((pop.means[, 3] - pop.sds[, 4]) < 0 &
# This is the upper bound of the 95% CI on the mean plus
#   the upper bound of the 95% CI on the sd 
               (pop.means[, 4] + pop.sds[, 4]) > 0)
# add to df
cred.ints$sd.creds <- sd.creds
# Add sd estimates to df
cred.ints$sd.ests <- round(pop.sds[, 1], digits = 2)
# print it
cred.ints
cat("\n");cat("\n")


# Close sink for now
sink()



# create plots for each individual, 
library("arm")

ordering <- sort.list(var_mean)
# figure output path
fig6path <- paste(outdir, "ERS1_step7_figure6_",
                  datestr, ".jpeg", sep = "")
jpeg(fig6path, height=176, width=12)
par(mgp=c(3, 1, 0), tck=-.01,
    mfrow = c(30, 3))
for(i in 1:88){
  beta <- extract(fit.lst[[i]][[1]])$beta
  var_mean <- apply(beta, 2, mean)
  var_sd <- apply(beta, 2, sd)
  coefplot(var_mean, sds=var_sd, 
           varnames=vars, 
           main=paste("Elk", i), 
           xlim=c(-3,3), mar=c(0,7,5,2), cex.main=.9)
}


dev.off()


################################################################################



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


pop.p <- ggplot(pop.results, aes(x=vars, y=mean.m)) + 
  geom_hline(yintercept = 0, color = "navy") +
  geom_point(color = "black") + 
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
                    ymax=pop.results$'95%.m' + pop.results$'95%.s'), color = "gray50",
                width = 0.7, size = 1, fatten = 0.5) + 
  geom_errorbar(data = pop.results, 
                aes(ymin=mean.m-mean.s, ymax=mean.m+mean.s), color = "gray75",
                width = 0.7, size = 1) +
  geom_jitter(data = mean.df_long, aes(x=variable, y=mean), color = "gray10",
              width = 0.3, size = 1) 


################################################################################
# Scatterplot figure of continuous coefficients
# save original graphing parameters
.og <- par(no.readonly = TRUE)

all.mean.arr <- mean.arr
# remove individual 78 for now
mean.arr <- mean.arr[-78,,]

fig1path <- paste(resultdir, "ERS1_step5_figure1_",
                  datestr, ".jpeg", sep = "")
jpeg(fig1path, units = 'in', res = 300,
     width = 11, height = 8.5)

par(mfrow = c(3, 3),
    mar = c(0.5, 2.1, 1.5, 1),
    bty = "l",
    oma = c(2.3, 0, 2.3, 0))

for(i in c(7, 9, 8, 6, 5, 1, 2, 3, 4)){
  # re-order array according to estimates of mean
  mean_df <- data.frame(mean = mean.arr[, i, 1],
                        se = mean.arr[, i, 2],
                        lo.hpdo = mean.arr[, i, 3],
                        lo.hpdi = mean.arr[, i, 4],
                        hi.hpdi = mean.arr[, i, 5],
                        hi.hpdo = mean.arr[, i, 6])
  # Create indiv id variable for merging
  mean_df$ind_id <- ind_id <- as.integer(gsub('elk', '', rownames(mean_df)))
  # merge with indiv info
  mean_df <- merge(mean_df, indiv_info, 'ind_id')
  # re-order decreasing means
  mean_df <- mean_df[order(mean_df$mean, decreasing = TRUE), ]
  
  attr <- vars[i]
  # define lims
  ylo <- min(mean_df$lo.hpdo)
  yhi <- max(mean_df$hi.hpdo)
  plot(mean_df$mean,
       ylim = c(ylo, yhi),
       xlab="", ylab="", 
       xaxt = "n", yaxt = "n",
       type = "n")
  mtext(attr, side = 3, line = 0.25, outer = FALSE, cex = 0.75)
  # define tick marks
  ll <- ylo
  ul <- yhi
  adj <- (ul - ll)/15
  at1 <- round(seq(ll + adj, ul - adj, 
               length.out = 4), digits = 2)
  
  # Add population information
  # get axis limits
  ext <- par('usr')
  # add rectangles for 90% HPD of pop mean
  rect(ext[1], pop.means[i, 5], ext[2], pop.means[i, 6],
       col = "grey75", border = NA)
  # add shading for 90% hpd of sds
  rect(ext[1],
       qnorm(0.05, pop.means[i, 5], pop.sds[i, 6]),
       ext[2],
       qnorm(0.95, pop.means[i, 6], pop.sds[i, 6]), 
       col = "grey90", border = NA)
  
  # add y axis
  axis(2, at = c(ext[3], ext[4]), labels=c("",""), lwd.ticks=0)
  axis(side = 2, at = at1, labels = T, tck = -0.015, mgp = c(3, 0.4, 0),
       cex.lab = 0.75)
  # add mean estimate
  abline(h = pop.means[i,1], lty = 1, col = "black", lwd = 3)
  # Add individual 95% hpds
  for(j in 1:87){
    lines(c(j, j), c(mean_df$lo.hpdo[j], mean_df$hi.hpdo[j]), lwd = 0.25)
    lines(c(j, j), c(mean_df$lo.hpdi[j], mean_df$hi.hpdi[j]), lwd = 0.5)
  }
  # Add point estimates
  colvec <- rep('pink', 88)
  colvec[as.factor(mean_df$sex) == 'Male'] = 'blue'
  points(mean_df$mean,
         pch = 21,
         cex = 0.6, bg = colvec, col = colvec)
}
mtext("Selection strength", side = 2, line = 0.75, outer = TRUE)
mtext("Individual", side = 1, line = 0.75, outer = TRUE)
mtext("Individual variation in selection of resource attributes", 
      side = 3, line = 0.75, outer = TRUE)

dev.off()


################################################################################
# Scatterplot figure of discrete coefficients + rxfire
# save original graphing parameters
.og <- par(no.readonly = TRUE)

fig2path <- paste(resultdir, "ERS1_step5_figure2_",
                  datestr, ".jpeg", sep = "")
jpeg(fig2path, units = 'in', res = 300,
    width = 11, height = 8.5)

par(mfrow = c(3, 3),
    mar = c(0.5, 2.1, 1.5, 1),
    bty = "l",
    oma = c(2.3, 0, 2.3, 0))

for(i in c(18, 11, 12, 13, 14, 15, 16, 17, 10)){
  # re-order array according to estimates of mean
  mean_df <- data.frame(mean = mean.arr[, i, 1],
                        se = mean.arr[, i, 2],
                        lo.hpdo = mean.arr[, i, 3],
                        lo.hpdi = mean.arr[, i, 4],
                        hi.hpdi = mean.arr[, i, 5],
                        hi.hpdo = mean.arr[, i, 6])
  # Create indiv id variable for merging
  mean_df$ind_id <- ind_id <- as.integer(gsub('elk', '', rownames(mean_df)))
  # merge with indiv info
  mean_df <- merge(mean_df, indiv_info, 'ind_id')
  # re-order decreasing means
  mean_df <- mean_df[order(mean_df$mean, decreasing = TRUE), ]
  
  attr <- vars[i]
  # define lims
  ylo <- min(mean_df$lo.hpdo)
  yhi <- max(mean_df$hi.hpdo)
  plot(mean_df$mean,
       ylim = c(ylo, yhi),
       xlab="", ylab="", 
       xaxt = "n", yaxt = "n",
       type = "n")
  mtext(attr, side = 3, line = 0.25, outer = FALSE, cex = 0.75)
  # define tick marks
  ll <- ylo
  ul <- yhi
  adj <- (ul - ll)/15
  at1 <- round(seq(ll + adj, ul - adj, 
                   length.out = 4), digits = 2)
  
  # Add population information
  # get axis limits
  ext <- par('usr')
  # add rectangles for 90% HPD of pop mean
  rect(ext[1], pop.means[i, 5], ext[2], pop.means[i, 6],
       col = "grey75", border = NA)
  # add shading for 90% hpd of sds
  rect(ext[1],
       qnorm(0.05, pop.means[i, 5], pop.sds[i, 6]),
       ext[2],
       qnorm(0.95, pop.means[i, 6], pop.sds[i, 6]), 
       col = "grey90", border = NA)
  
  # add y axis
  axis(2, at = c(ext[3], ext[4]), labels=c("",""), lwd.ticks=0)
  axis(side = 2, at = at1, labels = T, tck = -0.015, mgp = c(3, 0.4, 0),
       cex.lab = 0.75)
  # add mean estimate
  abline(h = pop.means[i,1], lty = 1, col = "black", lwd = 3)
  # Add individual 95% hpds
  for(j in 1:87){
    lines(c(j, j), c(mean_df$lo.hpdo[j], mean_df$hi.hpdo[j]), lwd = 0.25)
    lines(c(j, j), c(mean_df$lo.hpdi[j], mean_df$hi.hpdi[j]), lwd = 0.5)
  }
  # Add point estimates
  colvec <- rep('pink', 88)
  colvec[as.factor(mean_df$sex) == 'Male'] = 'blue'
  points(mean_df$mean,
         pch = 21,
         cex = 0.6, bg = colvec, col = colvec)
}
mtext("Selection strength", side = 2, line = 0.75, outer = TRUE)
mtext("Individual", side = 1, line = 0.75, outer = TRUE)
mtext("Individual variation in selection of resource attributes", 
      side = 3, line = 0.75, outer = TRUE)

dev.off()

# END OF SECTION.
################################################################################


################################################################################
# SECTION 2: MULTIDIMENSIONAL PREDICTOR COEFFICIENTS FIGURE 
################################################################################

# Remove individual with extreme observation of rxfire effect
mean.mat <- mean.mat[-78,]
se.mat <- se.mat[-78,]
  
# split matrices into continuous (v1) /habitat coefficients (v2)
mean_v1 <- mean.mat[, c(con.vars, idvec)]
mean_v2 <- mean.mat[, c(hab.vars, idvec)]

se_v1 <- se.mat[, c(con.vars, idvec)]
se_v2 <- se.mat[, c(hab.vars, idvec)]

# change identity vars to factors
mean_v1$sex <- as.factor(mean_v1$sex)
mean_v2$sex <- as.factor(mean_v2$sex)
mean_v1$cohort <- as.factor(mean_v1$cohort)
mean_v2$cohort <- as.factor(mean_v2$cohort)



################################################################################
# Create an ordination plot based on the estimates, and think about how 
#   it might be possible to incorporate uncertainty around those esimtates in.
#   Could be it's own subfigure, don't think anyone's done dimensional reduction
#   on estimates of uncertainty before. 
#
# Not sure what the final ordination will look like as of yet. Could be NMS, 
#   could be metric scaling, or just plain PCA. 
#
# Definitely need one ordination based on continuous predictors, and one based 
#   on habitat. 

# conduct PCA on continuous predictor coefficient mean estimates
# Add legend function, to add legend after plots have been set
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
# setup
scl <- 3
# same pch for each cohort
cohpchvec <- 21
sexcohvec <- c(21, 25)

# coors for "by cohort" biplot
cohcolvec <- c("chartreuse", "cyan", "darkorchid1")
# colors for "by age" biplot
agecolvec <- c("antiquewhite", "gold")
# make a cohort vector of factors
cohort <- as.factor(mean_v1$cohort)
# make a sex vector of factors
sex <- as.factor(mean_v1$sex)
# Make an effective age vector of factors
age <- as.factor(mean_v1$Age.eff)
# save original plotting options
.og <- par(no.readonly = TRUE)
#set image directory
fig3path <- paste(outdir, "ERS1_step7_results_figure3.tiff", sep = "")

# initialize jpeg plotting device
tiff(fig3path, units = 'in', res = 300,
     width = 8.5, height = 5.5)
# continuous response coefficients
data.mat <- mean_v1[, 1:9]
# prep data for PCA - scale but don't center, make data frame
X <- as.data.frame(scale(data.mat, center = FALSE))
# perform PCA
mod.con <- rda(X)
# Proportions of variance explained by each axis
prop.exp.con <- round(mod.con$CA$eig/sum(mod.con$CA$eig), digits = 3)
# set plotting options
par(mfrow = c(1, 2),
    mar = c(0.5, 2.1, 3, 1),
    bty = "o",
    oma = c(5.3, 1.2, 4.3, 1.2))
# X axis label
x.lab.con <- paste("PC 1, % Variance Explained:",
                   round(prop.exp.con[1], digits = 3)*100,
                   sep = " ")
# Y axis label
y.lab.con <- paste("PC 2, % Variance Explained:",
                   round(prop.exp.con[2], digits = 3)*100,
                   sep = " ")
# plotting continuous results
plot(mod.con, type = "n", scaling = scl,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
with(X, points(mod.con, display = "sites", bg = agecolvec[age],
               scaling = scl, pch = sexcohvec[sex], 
               cex = 0.75))
# Add continuous resource coefficients
text(mod.con, display = "species", scaling = scl, cex = 0.75, col = "darkcyan")
# draw axes lines and add axis titles
axis(1, at = c(-1, 0, 1), tck = -0.01, mgp = c(3, 0.4, 0))
axis(2, at = c(-1, 0, 1), tck = -0.01, mgp = c(3, 0.4, 0))
mtext(x.lab.con, side = 1, line = 1.25, cex.lab = 0.1, font = 2)
mtext(y.lab.con, side = 2, line = 1.25, cex.lab = 0.1, font = 2)
mtext("Biplot of PCA scores on individual RSF selection coefficients", 
      side = 3, line = -0.25, outer = TRUE)
mtext("(Intrinsic factors)", 
      side = 3, line = 1.70)

# Now ordered by cohort
plot(mod.con, type = "n", scaling = scl,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
with(X, points(mod.con, display = "sites", bg = cohcolvec[cohort],
               scaling = scl, pch = 21,
               cex = 0.7))
text(mod.con, display = "species", scaling = scl, cex = 0.75, col = "darkcyan")
axis(1, at = c(-1, 0, 1), tck = -0.01, mgp = c(3, 0.4, 0))
axis(2, at = c(-1, 0, 1), tck = -0.01, mgp = c(3, 0.4, 0))
mtext(x.lab.con, side = 1, line = 1.25, cex.lab = 0.1, font = 2)
mtext(y.lab.con, side = 2, line = 1.25, cex.lab = 0.1, font = 2)
mtext("(Cohort)", 
      side = 3, line = 1.70)

# vectors for drawing legend
leg.vec <- c("Ad. F", "Ad. M", "Sub-ad. F", "Sub-ad. M")
leg.ptbg <- c("antiquewhite", "antiquewhite", "gold", "gold")
leg.pch <- c(21, 25, 21, 25)
# draw legend
with(X, add_legend('topleft', legend = leg.vec, bty = "n", horiz = T,
                   pch = leg.pch, inset = c(0.1, 0.21),
                   pt.bg = leg.ptbg, 
                   text.width = c(0.1, 0.1, 0.1, 0.12),
                   cex = .75, pt.cex = 1.25))


# vectors for drawing legend
leg.vec <- c("2011", "2012", "2013")
leg.ptbg <- c("chartreuse", "cyan", "darkorchid1")
# draw legend
with(X, add_legend('topleft', legend = leg.vec, bty = "n", horiz = T,
                   pch = 21, inset = c(0.64, 0.21),
                   pt.bg = leg.ptbg,
                   text.width = 0.08,
                   cex = .75, pt.cex = 1.25))
dev.off()


################################################################################
# SECTION 4: PROPORTION IN TOP 95% OF MODELS 
################################################################################


vars.mc <- c("Distance to 2-track", "Distance to paved road", "Slope", 
          "Aspect", "Road Density", 
          "Canopy cover", "Distance to edge", "IJI", 
          "Years since Rx burn", "Habitat")


# Code from zonination at https://github.com/zonination/skittles/blob/master/skittles.R
# Set up an official palette
# Skittles Palette obtained by: http://www.color-hex.com/color-palette/1146
#   Colors:     Sberry    Orange    Lemon     Apple     Grape
#skipalette<-c("#c0043f","#e64808","#f1be02","#048207","#441349",
#              )

skipalette <- c("#CC3399", "#FF6633", "#FFCC00", "#009900", 
                "#6633CC", "#FF3399", "#FF9933", "#99FF33",
                '#33CCCC', '#0066CC')

# top model checkerboard plot, using package qdap
# need nvar x 88 length matrix
top.mat <- matrix(unlist(top.mods), ncol = 88)
# reorder according to cohort, sex, count
indiv_info$ind_id <- reorder(indiv_info$ind_id, order(indiv_info$cohort,
                                                      indiv_info$sex,
                                                      indiv_info$Age.eff,
                                                      decreasing = F))
idvec <- as.numeric(levels(indiv_info$ind_id))
top.mat <- top.mat[,idvec]

# take only first 10 variables for top models
top.vars.mat <- top.mat[1:10,]
top.vars.mat <- matrix(as.logical(top.vars.mat), ncol = 88)
dimnames(top.vars.mat)[[1]] <- vars.mc
# order by figure 4 main plot
top.var.counts <- rowSums(top.vars.mat)
var.ord <- order(top.var.counts, decreasing = TRUE)
top.vars.mat <- top.vars.mat[var.ord,]


# calculate proportions
#var.props <- lapply(var.sel, function(x){x/(sum(x))})
# create matrix of proportions in top 5
var.mat <- matrix(unlist(var.sel), ncol = 88)
var.df <- as.data.frame(var.mat)
# name columns and rows
colnames(var.df) <- paste(1:88, sep = ' ')
var.df$vars <- vars.mc

var.df.long <- melt(var.df, id.vars = "vars")

ord_vars_df <- as.data.frame(top.var.counts)
ord_vars_df$vars <- rownames(ord_vars_df)
var.df.long <- merge(var.df.long, ord_vars_df, by = 'vars')
names(var.df.long)[2] <- 'ind_id'
var.df.long$ind_id <- as.integer(var.df.long$ind_id)
final_dat <- merge(var.df.long, indiv_info, by = 'ind_id')
final_dat$vars <- factor(final_dat$vars, 
                         levels = vars.mc[var.ord])
final_dat$vars <- reorder(final_dat$vars, final_dat$top.var.counts,
                          order = T, decreasing = F)
# Reorder individuals by cohort, then sex, then GPS count
final_dat <- final_dat[order(final_dat$cohort, final_dat$sex, 
                             final_dat$Age.eff),]


# Rename data headers, and grab the overall mean
names(skittles)<-vars
netmean<-mean(c(skittles$Strawberry,skittles$Orange,skittles$Lemon,skittles$Apple,skittles$Grape))

# Summarize data into frame "skitsum" for the final plot.
coeffsum<-data.frame("Covariate"=names(skittles)[2:6],
                    "mean"=sapply(skittles[,2:6],mean),
                    "sd"=sapply(skittles[,2:6],sd))
skitsum$flavor <- factor(skitsum$flavor,
                         c("Strawberry","Orange","Lemon","Apple","Grape"))

# Generate "skitbin". "skitbin" will allow us to plot a neat "Raw Results" diagram in the first plot
skitmelt<-melt(skittles,id="Pak")
skitbin<-data.frame("Pak"=NA,"variable"=NA,"value"=NA)
for(n in 1:nrow(skitmelt)){
  skitbin<-rbind(skitbin,skitmelt[rep(n,skitmelt$value[n]),])   }
skitbin<-skitbin[2:nrow(skitbin),]
skitbin$variable <- factor(skitbin$variable,
                           c("Strawberry","Orange","Lemon","Apple","Grape"))
rm(n);rm(skitmelt)



# Convert raw data into heatmap, mildly more pleasing than "skitbin"
# figure 7 path
final_dat$value[final_dat$value == 103] = 77
final_dat$prop <- round(final_dat$value/77, digits = 1)
final_dat$prop <- substr(final_dat$prop, 2, 4)

ggplot(final_dat,aes(x=ind_id,y=vars))+
  geom_tile(aes(fill=vars,alpha=value),color="white")+
  geom_text(aes(label=substr(prop, 2, 4)),size=3)+
  scale_fill_manual("vars",values=rep("#6633CC", 10))+
  guides(alpha="none",fill="none")+
  scale_x_discrete(breaks=1:88)+
  scale_y_discrete(limits=levels(final_dat$vars))+
  labs(
       subtitle="Top 95% models",
       x="Individual elk",y="Covariate")+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 18))
setwd(outdir)
ggsave("ERS1_step7_results_figure7.png",height=4,width=16,dpi=100,type="cairo-png")



top.var.df <- as.data.frame(top.vars.mat)
top.var.df$vars <- factor(rownames(top.var.df),
                          levels = rownames(top.var.df))
top.var.df[top.var.df == T] = 1
final_top.df <- melt(top.var.df, id = "vars")
final_top.df$value <- as.numeric(final_top.df$value)

ggplot(final_top.df, 
            aes(x=variable,y=vars))+
  geom_tile(aes(fill=vars,alpha=value),color="white")+
  #geom_text(aes(label=prop),size=3)+
  scale_fill_manual("vars",values=rep("#009900", 10))+
  guides(alpha="none",fill="none")+
  scale_x_discrete(breaks=1:88)+
  scale_y_discrete(limits=levels(final_dat$vars))+
  labs(title="Parameter inclusion",
       subtitle="Top models",
       x="",y="Covariate")+
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 18),
        plot.title = element_text(size = 20))

setwd(outdir)
ggsave("ERS1_step7_results_figure8.png",height=4,width=16,dpi=100,type="cairo-png")



# top model checkerboard plot, using package qdap
# need nvar x 88 length matrix
top.mat <- matrix(unlist(top.mods), ncol = 88)
# reorder according to cohort, sex, count
indiv_info$ind_id <- reorder(indiv_info$ind_id, order(indiv_info$cohort,
                                                      indiv_info$sex,
                                                      indiv_info$count,
                                                      decreasing = F))
idvec <- as.numeric(levels(indiv_info$ind_id))
top.mat <- top.mat[,idvec]

# take only first 11 variables for top models
top.vars.mat <- top.mat[1:11,]
top.vars.mat <- matrix(as.logical(top.vars.mat), ncol = 88)
dimnames(top.vars.mat)[[1]] <- vars
# order by figure 4 main plot
top.var.counts <- rowSums(top.vars.mat)
var.ord <- order(top.var.counts, decreasing = TRUE)
top.vars.mat <- top.vars.mat[var.ord,]


# calculate proportions
var.props <- lapply(var.sel, function(x){x/(sum(x))})
# create matrix of proportions in top 5
prop.mat <- matrix(unlist(var.props), ncol = 88)
prop.df <- as.data.frame(prop.mat)
# name columns and rows
colnames(prop.df) <- paste(1:88, sep = ' ')
prop.df$vars <- vars
  
prop.df.long <- melt(prop.df, id.vars = "vars")

ord_vars_df <- as.data.frame(top.var.counts)
ord_vars_df$vars <- rownames(ord_vars_df)
prop.df.long <- merge(prop.df.long, ord_vars_df, by = 'vars')
names(prop.df.long)[2] <- 'ind_id'
prop.df.long$ind_id <- as.integer(prop.df.long$ind_id)
final_dat <- merge(prop.df.long, indiv_info, by = 'ind_id')
final_dat$vars <- factor(final_dat$vars, 
                            levels = vars[var.ord])
final_dat$vars <- reorder(final_dat$vars, final_dat$top.var.counts,
                          order = T, decreasing = F)
# Reorder individuals by cohort, then sex, then GPS count
final_dat <- final_dat[order(final_dat$cohort, final_dat$sex, 
                             final_dat$count),]

# Some publication-quality code with ggplot
# https://rpubs.com/Koundy/71792

theme_Publication <- function(base_size=14, base_family="serif") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
  + theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.text = element_text(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.key.size= unit(0.6, "cm"),
          legend.spacing = unit(0, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(1,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
  ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",
                 manual_pal(values = c("#386cb0","#fdb462",
                                       "#7fc97f","#ef3b2c",
                                       "#662506","#a6cee3",
                                       "#fb9a99","#984ea3",
                                       "#ffff33","#ff6633",
                                       "#009900")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#386cb0","#fdb462",
                                       "#7fc97f","#ef3b2c",
                                       "#662506","#a6cee3",
                                       "#fb9a99","#984ea3",
                                       "#ffff33")), ...)
  
}


# define color pallete
cbPalette <- c("#CC3399", "#FF6633", "#FFCC00", "#009900", 
               "#6633CC", "#FF3399", "#FF9933", "#99FF33",
               '#33CCCC', '#0066CC', '#6699CC')
mainplot <- ggplot(final_dat, aes(x=ind_id, y=value)) + 
            geom_bar(aes(fill = vars), stat = 'identity') +
            labs(x = 'Individual elk', y = 'Cumulative proportion of inclusion') +
            #scale_fill_manual(values=cbPalette, name = "Model\n coefficients") + 
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank())+
            scale_fill_Publication(name = "Model \n coefficients") + 
            theme_Publication() 

# make checkerboard plot
topplot <- qheat(t(data.frame(top.vars.mat)), by.column=NULL, 
      low = 'white', high="mediumseagreen", text.color =NA,
      grid='azure') + guides(fill=FALSE) +
      theme_Publication() + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(5, 50, 5, 25), 'mm')) 
# setting up layout params
layout.mat <- matrix(c(3,1, 1, 4, 2, 2, 2, 2), 
                     nrow = 2, byrow = TRUE)
nf <- layout(layout.mat)
fig4path <- paste(resultdir, "ERS1_step7_figure4_",
                  datestr, ".jpeg", sep = "")
jpeg(fig4path, width = 9, height = 5.5, units = 'in',
     res = 300)
grid.arrange(topplot, mainplot, nrow = 2, heights = c(2, 5))

dev.off()

################################################################################
# Barplot depicting the proportion of individuals for which that
#   parameter was in the top 95% of models ranked by WAIC
# Since discrete variables were in all top models..



# END OF SECTION.
################################################################################



################################################################################
# SECTION 5: PROBABILITY OF USE MAPS
################################################################################
# OPEN LIBRARIES
require(rstan)
require(foreach)
require(raster)
require(rgdal)
require(rasterVis)
require(ggplot2)
require(qdap)
require(viridis)
### DEFINE THE RESOURCE SELECTION FUNCTION

# Resource seleciton function
# This function creates a resource selection function out of a
#   stanfit object
# Args:
#   fit: The stanfit object, where parameter estimates are
#          ordered the same as the order of rasters in layer arg
#   layer: a raster or raster stack to use for predicting
#   param: a character vector of names of the parameters to be used in
#            predicting from the stanfit object
#   par.sel: an integer vector specifying which predictors to use               
# Returns:
#   pred.layer: a raster layer with predicted relative probabilities of use
#
rsf.predict <- function(fit, layer, par.name = c("mu", "stdev"), 
                        par.sel = 1:17){
  # For comparing to rownames, to decide which estimates
  #   are of means and which are stdev
  coeff.vec <- paste(par.name[1], "[", 1:100, "]", sep = "")
  #err.vec <- paste(par.name[2], "[", 1:100, "]", sep = "")
  # Extract model point estimates by comparing to these
  #   vectors
  results <- summary(fit)$summary
  means <- results[rownames(results) %in% coeff.vec,
                   c("mean")]
  #sds <- results[rownames(results) %in% err.vec, 
  #               c("sd")]
  
  # create an array from the raster stack
  lay_array <- as.array(layer)
  # create empty matrix for predicted probabilities
  pred.dim.r <- dim(lay_array)[1]
  pred.dim.c <- dim(lay_array)[2]
  predict.mat <- matrix(nrow = pred.dim.r,
                        ncol = pred.dim.c)
  # function to apply over each observation in the array
  # X is the vector slice containing all the observed resource
  #   values for a particular cell
  # m is the vector of meanvalues to be used, 
  #   must be supplied in call to apply
  rsf <- function(x, m){
    exp(x %*% m)/(1 + exp(x %*% m))
  }
  
  predict.mat <- apply(lay_array[, , par.sel], # raster data subset
                       MARGIN = c(1, 2), rsf, 
                       m = means) # mean parameters subset
  # make a raster from this matrix
  pred.rast <- raster(predict.mat)
  # then resmple it to get proper parameters
  #pred.layer <- resample(pred.rast, layer)
  return(pred.rast)
}


# Load individual model results
indiv.stan.results.path <- paste(resultdir, "ERS1_step7_results_2017-04-09", 
                                 sep = "")
load(indiv.stan.results.path)
# assign new name 
top.fitlst <- fitlst

# load population model results
pop.stan.results.path <- paste(resultdir, "ERS1_step2_results_2017-02-16", 
                               sep = "")
load(pop.stan.results.path)

# load predictors used in individual models
top.mods <- foreach(i = 1:88) %do% {
  # filepath
  path <- paste(step5dir, 'job4_results_elk_', i, '.rdata', sep = '')
  # load it
  load(path)
  # extract model data frame from first part of each list
  mod.df.lst <- lapply(fitlst, '[[', 1)
  # Remove timer element
  mod.df.lst[length(mod.df.lst)] <- NULL
  # number of models in model set (minus 1 for time calculation)
  M <- length(mod.df.lst)
  # unlist the model df list and populate matrix
  mod.comp.mat <- matrix(unlist(mod.df.lst), 
                         nrow = M, byrow = TRUE)
  # change to numeric
  mod.comp.mat <- apply(mod.comp.mat, 2, as.numeric)
  # discard all rows containing dist_2_gravel parameter
  mod.comp.mat <- mod.comp.mat[mod.comp.mat[, 5] == 1,]
  mod.comp.mat <- mod.comp.mat[,-5]
  # order and keep
  keep <- mod.comp.mat[order(mod.comp.mat[, 21]),][1,]
  # turn all '1' to 'FALSE'
  keep[keep == 1] = 0
  # turn all NA to TRUE
  keep[is.na(keep)] = TRUE
  keep
}

# list of integer vectors for predictors
top.preds <- lapply(top.mods, "[", 1:17) %>% lapply("==", 1) %>% lapply(which)

## load raster data
# raster names 
ras.hab.names <- c("forest", "wood", "sav", "shrub", "glade", "wsg", "csg",
             "foodplot")
# absolute paths
ras.hab.paths <- paste(rasdatadir, ras.hab.names, sep = '')
# Habitat and resource attribute rasters are of slightly different extent.
#   Need to stack separately and resample.
# Stack in order of cefficients estimated in the RSF
# Habitat raster stack
hab.stack <- stack(ras.hab.paths)
# Resource attribute stack
# rasters already scaled
ras.attr.names <- c('2track', 'paved', 'slope', 'aspect',  
                    'rd_den', 'canopy', 'edge', 'iji', 'rxfire')
# absolute paths
ras.attr.paths <- paste(rasdatadir, ras.attr.names, sep = '')
attr.stack <- stack(ras.attr.paths)

# VARIABLES ORDERED CANONICALLY ABOVE
#   can use this ordering to subset raster data individually



# Resample attributes to habitats
attr.stack.rs <- raster::resample(attr.stack, hab.stack, method = "ngb")
# Combine to final resource raster stack
res.stack <- stack(attr.stack.rs, hab.stack)

# Read ERZ boundary layer
setwd(erzlayerdir)
ERZ <- readOGR(dsn = ".", layer = "ERZ_bndry")
setwd(workdir)

# Re-project and mask the resource stack
temp.res <- stack(res.stack)
crs(temp.res) <- "+proj=utm +zone=15 +datum=NAD83 +units=m 
+no_defs +ellps=GRS80 +towgs84=0,0,0"
res.sub <- crop(temp.res, extent(ERZ))
res.ERZ <- mask(res.sub, ERZ)

## population level predictions
# Make predictions using the population-level RSF stanfit
pop.predictions <- rsf.predict(fit1, res.ERZ)

## Individual level predictions 
# start raster stack (for latter plotting)
indiv.s <- rsf.predict(fitlst[[1]], res.ERZ, par.name = "beta", 
                       par.sel = top.preds[[1]])
for(i in 2:88){
  indiv.preds <- rsf.predict(fitlst[[i]], res.ERZ, par.name = "beta",
                             par.sel = top.preds[[i]])
  indiv.s <- stack(indiv.s, indiv.preds)
}

# save individual predictions
indiv.preds.path <- paste(resultdir, "ERS1_step7_predictions.rda", 
                          sep = '')
saveRDS(indiv.s, file = indiv.preds.path)
# save population predictions
pop.predictions.path <- paste(resultdir, "ERS1_step7_pop_preds.rda",
                              sep = '')
saveRDS(pop.predictions, file = pop.predictions.path)

pred_corr_df <- data.frame(id = 1:88, corr_w_pop = rep(0, 88))
for(i in 1:88){
  comp_i <- stack(pop.predictions, indiv.s[[i]])
  pred_corr_df$corr_w_pop[i] <- rasterCorrelation(pop.predictions, indiv.s[[i]],
                                                  'pearson', na.rm = T)$
                                        'pearson correlation coefficient'[1, 2]
}

# individual predictions
indiv_preds.path <- paste(resultdir, "ERS1_step7_predictions.rda", sep = '')
indiv_preds_stack <- readRDS(indiv_preds.path)
# population predictions
pop_preds.path <- paste(resultdir, "ERS1_step7_pop_preds.rda", sep = '')
pop_preds_ras <- readRDS(pop_preds.path)

pred_corr_df <- data.frame(id = 1:88, corr_w_pop = rep(0, 88))
for(i in 1:88){
  b <- stack(pop_preds_ras, indiv_preds_stack[[i]])
  x <- as.matrix(b)
  pred_corr_df$corr_w_pop[i] <- cor(x, method = "spearman", 
                                    use = "complete.obs")[1, 2]
}

# add character string for interval of correlation
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.9] = "grt_eq_9"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.8 & 
                        pred_corr_df$corr_w_pop < 0.9] = "grt_eq_8"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.7& 
                        pred_corr_df$corr_w_pop < 0.8] = "grt_eq_7"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.6& 
                        pred_corr_df$corr_w_pop < 0.7] = "grt_eq_6"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.5& 
                        pred_corr_df$corr_w_pop < 0.6] = "grt_eq_5"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.4& 
                        pred_corr_df$corr_w_pop < 0.5] = "grt_eq_4"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.3& 
                        pred_corr_df$corr_w_pop < 0.4] = "grt_eq_3"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop >= 0.2& 
                        pred_corr_df$corr_w_pop < 0.3] = "grt_eq_2"
pred_corr_df$corr_cat[pred_corr_df$corr_w_pop < 0.2] = "lsth_2"

# save correlation results
corr.result.path <- paste(resultdir, "ERS1_step7_corr_results.rda", sep = '')
saveRDS(pred_corr_df, file = corr.result.path)

# plot pop level

pop.plot <- levelplot(pop_preds_ras, 
                      margin=FALSE,                       
                      colorkey=list(
                        space='bottom',                   
                        labels=list(at=-5:5, font=4),
                        axis.line=list(col='white'),
                        width=0.75
                      ),    
                      par.settings=list(
                        trellis.par.set("background", list(col = "white")),
                        strip.border=list(col='white'),
                        panel.background=list(col='white'),
                        
                        axis.line=list(col='transparent')
                      ),
                      scales=list(draw=FALSE),            
                      col.regions=viridis, 
                      begin = 0.6,
                      at=seq(0, 1, len=101))


## Plot 20 individual predictions
# subsample 20 IDs
elk.samp <- sample(88, 20)
# sumbsample corresponding rasters
indiv.s <- indiv_preds_stack[[elk.samp]]
# assign to plotting 
indiv.plot <- levelplot(indiv.s, 
                        margin=FALSE,                       
                        colorkey=FALSE,    
                        par.settings=list(
                          trellis.par.set("background", list(col = "white")),
                          strip.border=list(col='white'),
                          panel.background=list(col='white'),
                          plot.background = list(col="white"),
                          axis.line=list(col='white')
                        ),
                        scales=list(draw=FALSE),            
                        col.regions=viridis,
                        begin = 0.7,
                        at=seq(0, 1, len=101),
                        names.attr=rep('', nlayers(indiv.s))) 

fig9path <- paste(outdir, "ERS1_step7_results_figure9.jpg", sep = "")
jpeg(fig9path, units = 'in', res = 300,
     width = 6, height = 3.5)
trellis.par.set("background", list(col = "white"))
pop.plot
dev.off()

fig10path <- paste(outdir, "ERS1_step7_results_figure10.jpg", sep = "")
jpeg(fig10path, units = 'in', res = 300,
     width = 18, height = 18)
trellis.par.set("background", list(col = "white"))
trellis.par.set("strip.background", list(col = "white"))
indiv.plot
dev.off()

# pull out correlations
outpath <- paste(outdir, "fig_corr.Rout", sep = '')
sink(outpath)
pred_corr_df[elk.samp,]
sink()
################################################################################
# 



# END OF SECTION.
################################################################################



