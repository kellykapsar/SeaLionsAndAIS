


nInd <- 11

# specify directories
resultdirreg <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Results/"
resultdirreg <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Results/"


# Load best fitting model results from each model type
topfitreg <- readRDS(paste0(resultdir, "SSL_IndlAllCombos_2022-07-12/TopModelFits_ChiSquare.rds"))
topfitseason <- readRDS(paste0(resultdir, "SSL_IndlAllCombos_Seasonal_2022-06-16/TopModelFits_ChiSquare.rds"))
topfitswhr <- readRDS(paste0(resultdir, "SSL_IndlAllCombos_Swhr_2022-07-15/TopModelFits_ChiSquare.rds"))

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


P_labellist <- list() # Names of covariates in each top model (in order)
betanums <- list() # Order number for covariates in each top model (1-10)


# 
# Get list of all files in results directory
files <- list.files(resultdir)

megaloos <- list()
topmodnums <- list()

# Check Rhat values & neff
for(i in 1:nInd){
  # load stan model fit
  load(files[grepl(paste0("Ind", i, "_"), files)])
  
  # Extract summary info 
  loos <- lapply(fitlst, "[[", 2)
  megaloos[[i]] <- loos

  # Calculate top model numbers for each individual
  topmod <- as.data.frame(loo_compare(megaloos[[i]])) %>% tibble::rownames_to_column()
  topmod$modnum <- as.numeric(substr(topmod$rowname, 6,10))
  topmodnums[[i]] <- topmod$modnum[1]
}

# Identify covariates in each  top model 
for(i in 1:nInd){
  
  # Import variable combinations for this individual 
  mods <- read.csv(paste0(resultdir, "ModelVariableList_", i, ".csv")) %>% dplyr::select(-X)
  
  # Adjust column names 
  colnames(mods) <- covar_names
  
  # Identify covariates included in the top model 
  betatruefalse <- mods[topmodnums[[i]],]
  betanums[[i]] <- which(betatruefalse[1,] == TRUE)
  
  covar_l <- covar_labels[betanums[[i]],]
  covar_l$Parameter <- Parameters[1:length(betanums[[i]])]
  
  P_labellist[[i]] <- covar_l
}

################################################################################


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
  ggsave(paste0(resultdir, "CaterpillarPlot_", M_labels[i], ".png"), width=8, height=9, units="in")
}

### figure 2: Population & individual global model point estimates 

# Change Names, labels, and individual ids to factors 
topmodbetas$Name <- as.factor(topmodbetas$Name)
topmodbetas$Label <- as.factor(topmodbetas$Label)
topmodbetas$ind_id <- as.factor(topmodbetas$ind_id)

# Rename lower and upper CIs
topmodbetas <- topmodbetas %>% rename(CI_L = '2.5%', CI_H = '97.5%') 

# Reorder Label based on decreasing frequency
topmodbetas$Label <- reorder(topmodbetas$Label, topmodbetas$Name, FUN = length)
topmodbetas$Label <- factor(topmodbetas$Label, levels = rev(levels(topmodbetas$Label)))

# Number of parameters in top model for eachindividual 
topmodbetas %>% group_by(ind_id) %>% summarize(nparams=n())

# Identify number of individuals with each parameter in the top model 
topmodtotals <- topmodcounts %>% group_by(Label) %>% summarize(Count = sum(Count)) %>% arrange(-Count)


# saveRDS(topmodbetas, "./TopModBetas.rds")
# saveRDS(topmodcounts, "./TopModCounts.rds")
# saveRDS(topmodtotals, "./TopModTotals.rds")

## Create figure 2 with credible intervals on individual estimates


# initialize ggplot
pop.p <- ggplot(data = topmodbetas) +
  # add horizontal line for zero
  geom_hline(yintercept = 0, color = "navy", lwd = 1.5) + 
  geom_pointrange(aes(ymin = CI_L, ymax = CI_H, x = Label, y = mean, color=ind_id),
                  size = 0.8, position = position_jitterdodge(dodge.width=0.9)) + 
  scale_color_manual(values=getPalette, labels=M_labels) +
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
  annotate(geom="text", x = 1:10, y = rep(-3, 10), label=paste0("n = ",topmodtotals$Count), color="black", size=6) +
  
  ylab('Mean coefficient estimate') + 
  xlab('Resource covariate') + 
  labs(color = "Individual") +
  scale_y_continuous(limits = c(-3, 3)) 

# figure 2 filepath
fig1path <- "./JointCaterpillar.png"
ggsave(filename=fig1path, plot =  pop.p, width = 12, height = 8, units = "in")

