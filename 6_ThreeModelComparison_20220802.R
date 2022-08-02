################################################################################
# TITLE: Steller sea lion resource selection model fitting 
# DESCRIPTION: The purpose of this script is to create caterpillar plots for each
# indidivual Steller sea lion that compare the best fitting model for each model 
# type (weekly, seasonal, and seasonal with a weekly homerange). 
# CREATED BY: Kelly Kapsar
# DATE: 2022-08-02
# LAST UPDATED: 2022-08-02
################################################################################


nInd <- 11

# specify directories
resultdir <-  "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Results/ModelComparisons_2022-08-02/"
resultdirreg <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Results/SSL_IndlAllCombos_2022-07-12/"
resultdirseason <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Results/SSL_IndlAllCombos_Seasonal_2022-06-16/"
resultdirswhr <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Results/SSL_IndlAllCombos_Swhr_2022-07-15/"

# Load best fitting model results from each model type
topfitreg <- readRDS(paste0(resultdirreg, "TopModelFits_ChiSquare.rds"))
topfitseason <- readRDS(paste0(resultdirseason, "TopModelFits_ChiSquare.rds"))
topfitswhr <- readRDS(paste0(resultdirswhr, "TopModelFits_ChiSquare.rds"))

P_labellistreg <- readRDS(paste0(resultdirreg, "P_labellist.rds"))
P_labellistseason <- readRDS(paste0(resultdirseason, "P_labellist.rds"))
P_labellistswhr <- readRDS(paste0(resultdirswhr, "P_labellist.rds"))

################################################################################

M_labels = c("774PWS", "775PWS", "776PWS", "777PWS", 
             "781KOD", "782KOD", "783KOD", "784KOD", "785KOD", "786KOD", "788KOD")

# Get colors for individuals
getPalette = c("#0070ff", "#002673", "#b2df8a", "#33a02c",
               "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
               "#cab2d6", "#8967ae", "#d5d000")

# Covariate names (in order)
covar_names <- c("bathymetry", "dist_land", "dist_500m", "slope", "sst", "wind", "logship", 
                 "logfish", "prox_fish_km", "prox_ship_km")
covar_labels <- c("Bathymetry", "Distance to\nland", "Distance to\nshelf", 
                  "Slope ", "Sea surface\ntemp.", "Wind speed", "Non-fishing intensity", "Fishing\nintensity", 
                  "Distance to\nfishing", "Distance to\nnon-fishing")

###############################################################################


### figure 2: Population & individual global model point estimates 


topmodbetas <- data.frame()

for(i in 1:11){
  tempreg <- as.data.frame(rstan::summary(topfitreg[[i]], pars=c("beta"))$summary)
  tempreg$Parameter <- rownames(tempreg)
  tempreg$ind_id <- i
  tempreg$model_id <- 1 
  tempreg <- left_join(P_labellistreg[[i]], tempreg)
  topmodbetas <- rbind(topmodbetas, tempreg)
  
  tempseason <- as.data.frame(rstan::summary(topfitseason[[i]], pars=c("beta"))$summary)
  tempseason$Parameter <- rownames(tempseason)
  tempseason$ind_id <- i
  tempseason$model_id <- 2 
  tempseason <- left_join(P_labellistseason[[i]], tempseason)
  topmodbetas <- rbind(topmodbetas, tempseason)
  
  tempswhr <- as.data.frame(rstan::summary(topfitswhr[[i]], pars=c("beta"))$summary)
  tempswhr$Parameter <- rownames(tempswhr)
  tempswhr$ind_id <- i
  tempswhr$model_id <- 3 
  tempswhr <- left_join(P_labellistswhr[[i]], tempswhr)
  topmodbetas <- rbind(topmodbetas, tempswhr)
}

topmodbetas$sigpos <- ifelse(topmodbetas$`2.5%` > 0 & topmodbetas$`97.5%` > 0, 1, 0)
topmodbetas$signeg <- ifelse(topmodbetas$`2.5%` < 0 & topmodbetas$`97.5%` < 0, 1, 0)
topmodbetas$signonsig <- ifelse(topmodbetas$`2.5%` < 0 & topmodbetas$`97.5%` > 0, 1, 0)

# Change Names, labels, and individual ids to factors 
topmodbetas$Name <- as.factor(topmodbetas$Name)
topmodbetas$Label <- as.factor(topmodbetas$Label)
topmodbetas$ind_id <- as.factor(topmodbetas$ind_id)
topmodbetas$model_id <- as.factor(topmodbetas$model_id)

# Rename lower and upper CIs
topmodbetas <- topmodbetas %>% rename(CI_L = '2.5%', CI_H = '97.5%') 

# Reorder Label based on decreasing frequency
topmodbetas$Label <- reorder(topmodbetas$Label, topmodbetas$Name, FUN = length)
topmodbetas$Label <- factor(topmodbetas$Label, levels = rev(levels(topmodbetas$Label)))

# Number of parameters in top model for eachindividual 
topmodbetas %>% group_by(ind_id, model_id) %>% summarize(nparams=n())


## Create figure 2 with credible intervals on individual estimates

# Isolate out just caterpillar plots for presentation
for(i in 1:11){
  indbetas <- topmodbetas %>% filter(ind_id == i)
  # initialize ggplot
  pop.p <- ggplot(data = indbetas) +
    # add horizontal line for zero
    geom_hline(yintercept = 0, color = "navy", lwd = 1.5) + 
    geom_pointrange(aes(ymin = CI_L, ymax = CI_H, x = Label, y = mean, color=model_id),
                    size = 0.8, position = position_jitterdodge(dodge.width=0.9)) + 
    scale_color_manual(values=getPalette, labels=c("Weekly", "Seasonal", "Seasonal (weekly home range)")) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(color = "black", size = 16),
          axis.title = element_text(color = "black", size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black"), 
          legend.text = element_text(size=16), 
          legend.title = element_text(size=16)) + 
    geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), color= "lightgray") +
    ylab('Mean coefficient estimate') + 
    xlab('Resource covariate') + 
    labs(color = "Model") +
    scale_y_continuous(limits = c(-3, 3)) + 
    ggtitle(M_labels[[i]])
  
  # figure 2 filepath
  ggsave(paste0(resultdir, "AllModels_CaterpillarPlot_", M_labels[i], ".png"), width=12, height=8, units="in")
}
