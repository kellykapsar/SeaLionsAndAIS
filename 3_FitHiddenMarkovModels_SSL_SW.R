# Title: 1c_SSL_HMM_draft
# Author: Sydney Waloven
# Date: 2025-05-22
# Description: Script for determining behavioral states of the SSL data using a Hidden Markov Model


# setup -------------------------------------------------------------------

rm(list = ls())
library(moveHMM)
library(ggmap)

library(tidyverse)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Read in covariates 
landmask <- raster("../Data_Processed/Landmask_GEBCO.tif")
dist_500m <- raster("../Data_Processed/Dist500m.tif") %>%
  raster::mask(landmask, maskvalue = 1)
bathy <- raster("../Data_Processed/Bathymetry.tif") %>%
  raster::mask(landmask, maskvalue = 1)
dist_land <- raster("../Data_Processed/DistLand.tif") %>%
  raster::mask(landmask, maskvalue = 1)
slope <- raster("../Data_Processed/slope.tif") %>%
  raster::mask(landmask, maskvalue = 1)

# Create a list of static covariates 
staticcovars <- raster::stack(dist_land, dist_500m, bathy, slope)


# data import and formatting --------------------------------------------

# # Here I'm bringing in cleaned wildebeest data that has been resampled to 3 hr steps in the amt package. This is resampled data we saved from the earlier script, since we are investigating the influence of time of day here, and we needed regular sampling throughout the 24 period (even though animals were sampled at hourly intervals during the day). Remember this file now includes only animals from one of the three study areas.

# Bring in cleaned SSL data that has been resampled to 30 min time steps using the amt package (from script 1b_SSLDataProcessing_V2_SW.R)

raw <- read_rds("../Data_Processed/ssl_ak_30min.rds")

# Here we extract the value of each raster file for each tracking location. The extract_covariates() function in amt calls the terra::extract function.
ssl_data_cov <- raw %>% 
  extract_covariates(staticcovars) %>% 
  # Remove points that are on land
  filter(!is.na(Bathymetry))

# In terms of data prep...there are a few keys processing/formatting steps. 

# First it is assumed that location error is relatively small. If you are working with a high error data set, you should likely process it before using this package. For example, you may try a state space approach (e.g. crawl, aniMotum). 

# Second, note that the package does not use date/time information, and so assumes the steps are regular. It also assumes the data is organized chronologically, and this makes sense since the functions do not use any date/time information. So you have to arrange the data, and potentially resample the data (as we have done here) to ensure the steps are regular.

# Beyond that, we should convert any factor variable that we plan to investigate into a numeric, ordered variable. For example season names should be switched to 1, 2, 3, and 4. 

# Finally the animal id column must be named "ID".

# Importantly, I could not get this package to work properly without dividing my coordinates by 1000 to put them into km, but I could not figure out why exactly. The vignettes do this conversion, but simply say it is for convenience of units. However it appears if the data are in UTM coordinates, as ours are, you need to divide by 1000 to get this package to work.

# Here I will assign the data to an object and perform all the necessary reformatting. I'm converting to a data frame here because this is actually a track object from amt since we used that package for resampling. I'm also creating a column for hour of the day, which we will use below for modeling.

ssl_data <- ssl_data_cov %>%
  as.data.frame() %>% 
  arrange(deploy_id, t_) %>%
  select(x_,
         y_,
         t_,
         type,
         ID = deploy_id,
         slope,
         Dist500m,
         DistLand,
         Bathymetry) %>% 
  
  # create column for hour of the day and divide UTMs by 1000 to get km units 
  mutate(hour = lubridate::hour(t_),
         x_ = x_/1000,
         y_ = y_/1000)


# Create move object ------------------------------------------------------

# The function to prepare data to be analyzed in this package is called "prepData()". Here we could have used the raw lat/long values here as well, and then we would indicate the type = "LL".

ssl_move <- ssl_data %>% 
  prepData(type = "UTM",
           coordNames = c("x_","y_"))

class(ssl_move)
summary(ssl_move) # error

# We now have movement data for each individual sea lion. There are two new columns "step" and "angle" which represent the step lengths and turning angles.

# Overall plot of dataset. Multiple animals are plotted in separate windows.
plot(ssl_move,
     ask = FALSE)
# if ask = TRUE each plot will come separately and you'll need to hit enter, compact = T plots all tracks at once

# Or plot them separately
plot(ssl_move[ssl_move$ID == "SSL2018774PWS", ],
     ask = FALSE)

# Look at step and turning angle
hist(ssl_move$step)
summary(ssl_move$step)

# This indicates that one animal moved 108 km in one step (?)

quantile(ssl_move$step, 
         probs = 0.90, 
         na.rm = T) 
# 90% of steps below 6772m.

hist(ssl_move$angle)

# Reminder these plots are using data from all animals.


# Fit Hidden Markov Models - Set starting values --------------------------

# We need to set starting values for all relevant parameters that will be estimated so that the optimization algorithm has a starting point. The initial parameters will be specified in two vectors, stepPar0 (step distribution) and anglePar0 (angle distribution). Start with basic default values that often work well for animal movement with 2 behavioral states.

# Decide what statistical distribution will be used to describe the step lengths and turning angles. Different distributions have different numbers of parameters.

# We'll start with a gamma for step lengths and von Mises distribution for turning angles. For a gamma distribution, we will need a mean and SD (and zero-mass if there are any step lengths = 0). For the von Mises, we will need a mean and concentration parameter. Zero-inflation must be included in the step length distribution if some steps are of length exactly zero. Add another parameter to the distribution: its mass on zero. 

# See if there are any animals with zero step lengths:
# Look at the 10 lowest values of step. If "with ties = T" the function will not count ties for step value and so you'll end up with more than 10 rows in the output.
slice_min(ssl_move,
          order_by = step,
          n = 30, # increase to 30 to see all 24 rows with zero step length 
          with_ties = F)
# Since we have just one value of zero here, we'll need to include this extra parameter. 

# Here, the initial values are chosen such that they correspond to the commonly observed pattern in 2-state HMMs for animal movement data, with state 1 involving relatively short steps and many turns (hence the choice of a small initial value for the mean of the gamma step length distribution and an initial value of pi radians (180 degrees) for the mean turning angle) and state 2 involving longer steps and fewer turns (hence the choice of a larger initial value for the mean of the gamma step length distribution and an initial value of 0 for the mean turning angle).

# Note that some authors, to ensure starting values did not overly impact results, will repeat model fitting with a range of starting values. This allows one to test the sensitivity of results to starting value decisions.

## Step distributions

# Step mean (two parameters: one for each state)
mu0 <- c(0.1, 1) 

# Step standard deviations
sigma0 <- c(0.1, 1)

# Zero-distribution term. 
zeromass0 <- c(0.1, 0.05) # Only one zero.
# zeromass0 <- c(0.15, 0.1) # Since there are more frequent zeros, we'll increase.

# Assigning the step distribution starting values
stepPar0 <- c(mu0,
              sigma0,
              zeromass0)

## Turning angle starting values

# Mean turning angle. In radians, pi, or 3.14 represents 180 degrees.
angleMean0 <- c(pi, 0)

# Angle concentration. This reflects the variance in the turning angle distribution.
kappa0 <- c(1, 1) 

# Setting selected parameters
anglePar0 <- c(angleMean0,
               kappa0)

## Below are parameters for a 3 state behavior model as well, one value for each behavior state.
mu0_3 <- c(0.1, 0.5, 3)
sigma0_3 <- c(0.05, 0.5, 1)
zeromass0_3 <- c(0.05, 0.0001, 0.0001)
stepPar0_3 <- c(mu0_3, sigma0_3, zeromass0_3)
angleMean0_3 <- c(pi, pi, 0)
kappa0_3 <- c(1, 1, 1)
anglePar0_3 <- c(angleMean0_3, kappa0_3)

# For numerical stability, we can standardize the covariate values before fitting the model. But we currently won't be combining them in any models.


# Fit Hidden Markov Models - Model Fitting --------------------------------

# When fitting models we can provide a formula with a covariate to inform transition from one state to the next. Specify how many behavioral states we want the model to identify with "nbStates = 2", which indicates we want to fit a two-state model to the data. Next, define the distributions that we want to use to characterize both the step lengths (stepDist = "gamma") and turning angles (angleDist = "vm").

# Distribution options for step length are: gamma (“gamma”), Weibull (“weibull”), exponential (“exp”) and log-normal (“lnorm”).

# For turning angle we can choose: von Mises (“vm”) and wrapped-Cauchy (“wrpcauchy”). It is also possible to specify angleDist = "none", if the angles are not modeled.

# Fit a Null model with two behavioral states but no covariates influencing transitions from one state to the next. Formula has an intercept only (~ 1) to follow standard regression formula specification.
ssl_m_null <- fitHMM(data = ssl_move, 
                      nbStates = 2, # number of behavioral states
                      stepPar0 = stepPar0, 
                      anglePar0 = anglePar0, 
                      stepDist = "gamma",
                      angleDist = "vm",
                      formula = ~ 1) # intercept only

ssl_m_null

# State 1 has a smaller mean step length than state 2 around 105m (may vary with each run), with mean for state 2 around 491. The zero-mass values are both very small since there were not that many zero step length values.

# The concentration represents a sort of variance around the mean turning angles. Both states have a mean turning angle close to 0. A relative turn angle of zero means continued movement in the same direction.

# There is only an intercept regression coefficient since we did not test a regression equation.

# The transition matrix shows the average probability of transitioning from state 1 to state 2 is 10.8% and transitioning from state 2 to state 1 is also 10.8%. This indicates that once an animal is in a state it tends to stay there and transitions between states are relatively rare.

# View plots. State 1 is orange and state 2 is blue. Path segments are colored according to predicted state in all the full trajectory plots of each animal. Longer distance straight movements are colored as state 2 (blue) and more resident/locally focused movements as state 1 (orange). It is likely here that state 1 is more of a local foraging state while state 2 is directed travel. 
plot(ssl_m_null,
     ask = F)
# hit the back button in the plotting window to see all plots


# State Assignments and Probabilities -------------------------------------

# To globally decode the state process, the Viterbi algorithm is implemented in the function viterbi(). This function outputs the most likely sequence of states to have generated the observations, under the fitted model. 

ssl_states_2 <- viterbi(ssl_m_null)

# Look at the most probable states for the first 25 observations of the first individual.
ssl_states_2[1:25]
# The first 6 are state 1 and the rest are state 2.

# ssl_states_2 gives us a predicted state for each location in the dataset. This allows us to look at the proportion of time spent in each state.
prop.table(table(ssl_states_2))

# Overall, these sea lions spent about half their time in a resident state and half their time in a more nomadic or migratory state.

stateProps <- ssl_data %>% 
  
  # Add predicted state information to initial data frame
  mutate(state = ssl_states_2) %>% 
  
  # Add new column that is the total locations for each animal
  mutate(locs = n(),
         .by = ID) %>% 
  
  # Summarize for each animal, and each state, the proportion of locations
  summarize(stateProp = n()/locs,
            .by = c(ID, state)) %>% 
  
  # Unique rows of information only
  distinct()

# You can compare the proportions of each states visually to the maps to see how it matches up.

# Plot this across individuals 
stateProps %>% 
  filter(state == 1) %>% 
  mutate(ID = fct_reorder(ID, # Order the ID column by the value of stateProp.
                          stateProp)) %>% 
  
  ggplot(aes(y = stateProp,
             x = ID,
             fill = ID)) +
  geom_col(col = "black",
           position = position_dodge()) +
  coord_flip() +
  labs(y = "Prop. time in State 1 (foraging/encamped)") +
  theme_bw()

# SSL2019785KOD spends the most time in this local resident state. Let's look at their trajectory again.
plot(ssl_move[ssl_move$ID == "SSL2019785KOD", ],
     ask = F)

# Differentiating between individual animals in terms of overall movement strategy can be better done with a net-squared displacement model or even mean-squared displacement. HMM are used more to categorize segments as particular behavioral states.

# * Could further specify this analysis by splitting different times of year or different seasons?

# Look at actual probabilities assigned to each location for each behavioral state
ssl_probs_2 <- stateProbs(ssl_m_null)

# Output is a matrix. One row for each location. Values across a row should add to 1.0.
head(ssl_probs_2)
nrow(ssl_probs_2)
# The state with the highest probability according to stateProbs() might not be the same as the state in the most probable sequence returned by the Viterbi algorithm. This is because the Viterbi algorithm performs "global decoding" and the state probabilities are "local decoding".

# Visualize the results of viterbi() and stateProbs(). This plot shows the plots of the most likely state sequence decoded by the Viterbi algorithm, as well as both columns of the matrix of state probabilities, for one individual "SSL2019785KOD".
plotStates(ssl_m_null,
           animals = "SSL2019785KOD",
           ask = F)

# This information can be associated with all kinds of variables (e.g., temperature, NDVI, breeding state, etc)

# Custom plot of same information
ssl_data %>% 
  mutate(state = as.factor(ssl_states_2)) %>% 
  filter(ID == "SSL2019785KOD") %>% 
  ggplot(aes(x = x_,
             y = y_,
             fill = state,
             col = state)) +
  geom_path(alpha = 0.5) +
  geom_point(shape = 21,
             alpha = 0.8,
             col = "black") +
  scale_fill_manual(values = c("orange",
                               "cornflowerblue")) +
  scale_color_manual(values = c("orange",
                                "cornflowerblue")) +
  theme_minimal() +
  labs(x = "x",
       y = "y",
       title = "SSL2019785KOD")


# Model Comparisons and Covariate Testing ---------------------------------

# Look at the influence of a specific covariate
ssl_m_bathy <- fitHMM(data = ssl_move,
                      nbStates = 2,
                      stepPar0 = stepPar0,
                      anglePar0 = anglePar0,
                      formula = ~ Bathymetry)
