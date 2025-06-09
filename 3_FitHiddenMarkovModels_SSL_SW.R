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


# data import and formatting --------------------------------------------

# # Here I'm bringing in cleaned wildebeest data that has been resampled to 3 hr steps in the amt package. This is resampled data we saved from the earlier script, since we are investigating the influence of time of day here, and we needed regular sampling throughout the 24 period (even though animals were sampled at hourly intervals during the day). Remember this file now includes only animals from one of the three study areas.

# Bring in cleaned SSL data that has been resampled to 30 min time steps using the amt package (from script 1b_SSLDataProcessing_V2_SW.R)

read_rds("../Data_Processed/ssl_ak_30min.rds")

# In terms of data prep...there are a few keys processing/formatting steps. 

# First it is assumed that location error is relatively small. If you are working with a high error data set, you should likely process it before using this package. For example, you may try a state space approach (e.g. crawl, aniMotum). 

# Second, note that the package does not use date/time information, and so assumes the steps are regular. It also assumes the data is organized chronologically, and this makes sense since the functions do not use any date/time information. So you have to arrange the data, and potentially resample the data (as we have done here) to ensure the steps are regular.

# Beyond that, we should convert any factor variable that we plan to investigate into a numeric, ordered variable. For example season names should be switched to 1, 2, 3, and 4. 

# Finally the animal id column must be named "ID".

# Importantly, I could not get this package to work properly without dividing my coordinates by 1000 to put them into km, but I could not figure out why exactly. The vignettes do this conversion, but simply say it is for convenience of units. However it appears if the data are in UTM coordinates, as ours are, you need to divide by 1000 to get this package to work.

# Here I will assign the data to an object and perform all the necessary reformatting. I'm converting to a data frame here because this is actually a track object from amt since we used that package for resampling. I'm also creating a column for hour of the day, which we will use below for modeling.

ssl_data <- read_rds("../Data_Processed/ssl_ak_30min.rds") %>%
  as.data.frame() %>% 
  arrange(deploy_id, t_) %>%
  select(x_,
         y_,
         t_,
         type,
         ID = deploy_id) %>% 
  
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
# zeromass0 <- c(0.1, 0.05) If there were only one zero. 
zeromass0 <- c(0.15, 0.1) # Since there are more frequent zeros, we'll increase.  

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