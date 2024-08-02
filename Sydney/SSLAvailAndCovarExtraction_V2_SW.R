# Title: SSLAvailAndCovarExtraction
# Author: Kelly Kapsar
# Date: 2024-07-22
# Description: R script version of 2_SSLAvailAndCovarExtraction.Rmd for readability


# setup -------------------------------------------------------------------

library(sf)
library(raster)
library(scales)
library(ggmap)
library(amt) 
library(tidyverse)

# Read in data -----

# Set seed
set.seed(1015)

# Projection information for WGS84/UTM Zone 5N (EPSG:32605)
prj <- 32605

# Create study area polygon
coords <- data.frame(lat = c(56, 62, 62, 56, 56), 
                     lon = c(-155, -155, -143, -143, -155),
                     id = "study")

study <- coords %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  group_by(id) %>% 
  summarize(geometry = st_combine(geometry)) %>%  
  st_cast("POLYGON") %>% 
  st_transform(prj)

# Read in clean, non-land ssl used locations
# ssl <- readRDS("../Data_Processed/Telemetry/watersealis.rds")
# ssl$used <- 1

# # Read in resampled track rds - each animal's data is nested in a df row. This contains clean, non-land ssl used locations transformed to an amt track.
trk <- read_rds("../Data_Processed/ssl_steps_resampled.rds") %>% 
  # Unnest data
  select(deploy_id, data, steps) %>% 
  unnest(cols = data) %>% 
  
  mutate(
    # Add column marking these points as used locations
    used = 1,
    
    # Create an ID column that identifies each unique individual/year/week combination
    weeklyhr_id = factor(paste0(.$deploy_id, 
                                substr(.$weekofyear, 1, 4), 
                                substr(.$weekofyear, 7, 8)))
  )

# Save data
trk %>% 
  # Remove steps column before saving
  select(-steps) %>% 
  st_write("../Data_Processed/Telemetry/UsedLocsTrk_Clean.shp")

# Create an ID column that identifies each unique individual/year/week combination
# ssl$weeklyhr_id <- factor(paste0(ssl$deploy_id, 
#                                  substr(ssl$weekofyear, 1, 4), 
#                                  substr(ssl$weekofyear, 7, 8)))
# st_write(ssl, "../Data_Processed/Telemetry/UsedLocs_Clean.shp")

# Read in covariates 
landmask <- raster("../Data_Processed/Landmask_GEBCO.tif")
dist_500m <- raster("../Data_Processed/Dist500m.tif") %>% 
  raster::mask(landmask, maskvalue = 1)
depth <- raster("../Data_Processed/Bathymetry.tif") %>% 
  raster::mask(landmask, maskvalue = 1)
dist_land <- raster("../Data_Processed/DistLand.tif") %>% 
  raster::mask(landmask, maskvalue = 1)
slope <- raster("../Data_Processed/slope.tif") %>% 
  raster::mask(landmask, maskvalue = 1)

ship <- readRDS("../Data_Processed/AIS_AllOther.rds") 
fish <- readRDS("../Data_Processed/AIS_Fishing.rds")
sst <- readRDS("../Data_Processed/sst_weekly.rds")
wind <- readRDS("../Data_Processed/wind_weekly.rds")

# Create a list of all covariate files and check resolution 
raslist <- list(depth, dist_land, dist_500m, slope, ship, fish, sst, wind)
rasres <- lapply(raslist, function(x) res(x)/1000)

# Extract dynamic covariate data at used locations -----
### Modified FROM AMT PACKAGE
# https://github.com/jmsigner/amt/blob/master/R/extract_covariates.R
# Also helpful: https://cran.r-project.org/web/packages/amt/vignettes/p3_rsf.html
# I modified it so that it works with weekly timescale data 
# DOES NOT work with holes in the data set any more 

# extract_covar_var_time_custom <- function(xy, t, covariates) {
#   
#   t_covar <- raster::getZ(covariates) # Get time slices from covariates
#   t_obs <- format(as.POSIXct(t), "%G-W%V") # Convert timestamp to week - same week of year as lubridate::isoweek(ssl_dates$date)
#   
#   wr <- sapply(t_obs, function(x) which(x == t_covar)) # Identify which slice to select
#   ev <- terra::extract(covariates, xy) # Extract covariate values for all time slices at all used locations 
#   cov_val <- ev[cbind(seq_along(wr), wr)] # Select only the relevant time slice for each location
#   
#   return(cov_val) # Return the extracted covariate values
# }

# Weekly KDE home ranges --------------------------------------------------

# # Isolate out variables of interest 
ssl_simple <- trk %>% 
  select(weeklyhr_id, 
         t_, # date
         x_, # lon
         y_) %>% # lat
  # Calculate time of day based on lat/lon and timestamp
  time_of_day() # adds tod_ column to end of df

# Create a list of static covariates 
staticcovars <- raster::stack(dist_land, dist_500m, depth, slope)

# Calculate the number of rows that share the same weeklyhr_id (used points)
pointcounts <- trk %>% 
  st_drop_geometry() %>% 
  group_by(weeklyhr_id) %>% 
  summarize(n = n()) %>%
  ungroup()

# Create empty dfs
ssl_hr <- data.frame()
avail_pts <- data.frame()


# Identifying available points --------------------------------------------

# Here we are generating 10 times the number of "use" points to get our "available" points as a starting point. We will need to assess what a robust minimum # of available points should be. 

# random_points() places the random/regular points within a 100% MCP, though not explicitly stated. No other options currently available. We will use hr_kde to manually set our home range.

# Function to generate random points for one animal
generate_random_points <- function(track_data) {
  
  # Calculate home range (KDE)
  hr_kde <- hr_kde(track_data, levels = c(0.95))
  
  # Generate random points within the home range
  random_pts <- random_points(hr_kde,
                              n = nrow(track_data) * 10, 
                              presence = track_data)
  
  
  # Add animal ID to the random points
  random_pts <- random_pts %>%
    mutate(id = unique(track_data$weeklyhr_id))
  
  return(random_pts)
}

# Split track data by animal ID
split_trk <- split(ssl_simple, ssl_simple$weeklyhr_id)

# Generate random points for each animal and combine the results. This keeps the initial data points in the same df, but keeps track of which points are available vs used by adding a "case_" column. TRUE = use. FALSE = random.
ssl_rsf_10 <- map_dfr(split_trk, generate_random_points)

table(ssl_rsf_10$case_)

# Plot use/availability
terra::plot(ssl_rsf_10)

# Assign weights to available points that are very high.
ssl_rsf_10 <- ssl_rsf_10 %>% 
  mutate(weight = ifelse(case_ == T,
                         1,
                         5000))

# The argument type = "regular" samples points in a regular grid-like pattern instead of completely randomly, beneficial for spatial analysis requiring uniform coverage. The argument presence = track_data generates points that reflect the sampling distribution or habitat patterns in the tracking data. Which should we use?
# When using type = "regular", the plots of simulated random vs use points has the random points completely covering the use points. When the presence argument is used, the random vs use points are both plotted visibly.

# Here we extract the value of each raster file for each tracking location and each random location. The extract_covariates() function in amt calls the terra::extract function.
ssl_rsf_10 <- ssl_rsf_10 %>%
  extract_covariates(staticcovars) %>% 
  # Remove points that are on land
  filter(!is.na(Bathymetry))

# Make vector of variable names
vars <- ssl_rsf_10 %>% 
  select(DistLand:slope) %>% 
  names()

# Histograms for each static variable
hist_plots <- map(vars, ~
                    
                    ssl_rsf_10 %>% 
                    select(case_, var = .x) %>% 
                    ggplot(aes(x = var,
                               after_stat(density),
                               col = case_)) + 
                    geom_freqpoly(size = .7,
                                  bins = 30) + 
                    labs(title = .x))

# Put all plots into one grid. Expand to view clearly.
cowplot::plot_grid(plotlist = hist_plots)


# Boxplots for each static variable
box_plots <- map(vars, ~
                   
                   ssl_rsf_10 %>% 
                   select(case_, var = .x) %>% 
                   ggplot(aes(y = var,
                              col = case_)) + 
                   geom_boxplot() + 
                   labs(title = .x))

cowplot::plot_grid(plotlist = box_plots)


# Assess Collinearity -----------------------------------------------------

# Collinear predictors influence the variance estimates of the model and makes it difficult to interpret model coefficients. Determining what covariates are important requires avoiding multi-collinearity. 

# Generally correlation values over 0.70 are problematic.

cor(as.data.frame(ssl_rsf_10[ , 6:9]))

# Check the Variance Inflation Factor. This assesses the multicollinearity across ALL the predictors, with higher values indicating high collinearity of a predictor with the rest of the set. Some sources indicate a VIF over 10 is something for concern, while other indicate 3.0 should be the threshold. 

usdm::vif(as.data.frame(ssl_rsf_10[ , 6:9])) # all < 3.0

# Sensitivity Analysis ----------------------------------------------------

# Determine how many available points are needed to get stable coefficient estimates.

# This will look at 1, 5, 20, 50, and 100 available points for each use location. 
n.frac <- c(1, 5, 20, 50, 100)

# Total available locations to be generated, based on the number of use points
n.pts <- nrow(trk) * n.frac

# Set the number of replicate runs of the model to estimate variability of the coefficient estimate across the runs
n.rep <- 100

# To run this we are creating a table that holds all the settings of each scenario. The result column holds the model results. This code is based off the Fieberg et al. 2021 publication code. For each scenario we are actually reextracting covariate values, since we have different sets of random locations each time, and running a full model.

# The model run itself is done in a piped chain, where we create the track object, create the random points, extract covariate values, scale our covariates, and run the model. We will scale our covariates in the model state later on. This takes approximately 5 hours to run.

wb.sim <- tibble(
  n.pts = rep(n.pts, n.rep),
  frac = rep(n.frac, n.rep),
  result = map(
    n.pts, ~
      trk %>% 
      random_points(n = .x) %>%
      mutate(w = ifelse(case_ == T, 
                        1, 
                        5000)) %>% 
      extract_covariates(staticcovars) %>%
      mutate(dist_land = scale(DistLand),
             dist_500m = scale(Dist500m),
             depth = scale(Bathymetry),
             slope = scale(slope)) %>%
      glm(case_ ~ dist_land + dist_500m + depth + slope,
          data = ., 
          weights = w,
          family = binomial(link = "logit")) %>%
      broom::tidy()))

# Good idea to save this so you have it and don't need to rerun in later.
# write_rds(wb.sim, file = "../Data_Processed/ssl_sa_sim.rds")

# Read in rds
wb.sim <- read_rds("../Data_Processed/ssl_sa_sim.rds")

# This is an object structure we already worked with in the amt package in week 1. We have a full list which contains the model result for each run, and that is saved for each of the 100 rows of our data frame. There are 100 rows (5 different amounts of available locations, 20 model runs each)

# We can look at the result of one model run like this
wb.sim$result[[1]]

# Now to plot the results and make some conclusions from all this. We now "unnest" the results to plot the coefficient estimates from individual model fits with different sets of available data (1, 5, 20). The unnesting results in a much longer data frame. Run that first part of the code to see the result. After that we are recoding the "term" column, essentially changing the values in the column. We are plotting here the estimate, which is the beta coefficient from the glm regression for each covariate.

wb.sim %>% unnest(cols = result) %>% 
  mutate(term = case_match(term, 
                           "(Intercept)" ~ "Intercept",
                           "anth_risk" ~ "Human disturbance", 
                           "woody_dist" ~ "Woody dist.",
                           "fence_dist" ~ "Fence dist.", 
                           "prirds_dist" ~ "Primary Rd dist.",
                           "river_dist" ~ "River dist.",
                           "secrds_dist" ~ "Sec Rd dist.",
                           "waterpts_dist" ~ "Water Point dist.")) %>% 
  
  ggplot(aes(x = factor(frac), 
             y = estimate)) +
  geom_boxplot() + 
  facet_wrap(~ term, 
             scale  ="free") +
  geom_jitter(alpha = 0.2) + 
  labs(x = "Avail. points per use location", 
       y = "Estimate") +
  theme_light()

# Pretty snazzy right? This is a key tool in making a final assessment of the minimum # of available points needed for each use location. We can see little change in the coefficient values once we have 20 points per use location, but we'll go with 50 to be safe.


# Kelly's code ------------------------------------------------------------

# # 
ssl_new <- trk %>% 
  select(deploy_id, 
         weeklyhr_id,
         t_,
         year, 
         weekofyear, 
         y_, # # northing
         x_, # # easting
         used) %>% # # 
  cbind(., raster::extract(staticcovars, .)) %>% 
  rename(dist_land = DistLand, 
         dist_500m = Dist500m)

# Extract week of year as date from original data 
ssl_dates <- ssl %>% 
  st_drop_geometry() %>% 
  group_by(weeklyhr_id) %>% 
  summarize(date = min(date))

avail_pts <- left_join(avail_pts, 
                       ssl_dates, 
                       by = c("weeklyhr_id" = "weeklyhr_id"))

# Check that week of year re-calculation worked okay 
# It works!
# ssl_dates$weekofyear <- lubridate::isoweek(ssl_dates$date)
# ssl_dates$weekofyear2 <- format(ssl_dates$date, "%G-W%V")
# ssl_datetest <- left_join(ssl, 
#                           ssl_dates, 
#                           by = "weeklyhr_id")
# length(which(ssl_datetest$weekofyear.x != ssl_datetest$weekofyear.y))

# Clean up column names in available points 
avail_pts <- avail_pts %>% 
  mutate(deploy_id = substr(avail_pts$weeklyhr_id, 1, 13), 
         year = lubridate::year(avail_pts$date),
         weekofyear = lubridate::isoweek(avail_pts$date), 
         used = 0) %>% 
  rename(northing = y_, 
         easting = x_, 
         dist_land = DistLand,
         dist_500m = Dist500m) %>%
  select(-case_)


# Merge used with available 
all_pts <- avail_pts %>% 
  st_as_sf(., coords = c("easting", "northing"), 
           crs = prj, 
           remove = FALSE) %>% 
  rbind(., ssl_new)

# Rename to match other naming schemes
ssl4 <- all_pts

# Extract dynamic weekly covariate values at used and available locations 
# Extract covariate data at used locations 
ssl4$sst <- extract_covar_var_time_custom(xy = st_coordinates(ssl4), 
                                          t = ssl4$date,
                                          covariates = sst)
ssl4$wind <- extract_covar_var_time_custom(xy = st_coordinates(ssl4), 
                                           t = ssl4$date, 
                                           covariates = wind)
ssl4$ship <- extract_covar_var_time_custom(xy = st_coordinates(ssl4), 
                                           t = ssl4$date, 
                                           covariates = ship)
ssl4$fish <- extract_covar_var_time_custom(xy = st_coordinates(ssl4),
                                           t = ssl4$date,
                                           covariates = fish)

# convert column names to lowercase 
colnames(ssl4) <- tolower(colnames(ssl4)) 

save(ssl4, file = "../Data_Processed/Telemetry/TEMP_20220706.rda")
# browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ&list=RDK1b8AhIsSYQ&start_radio=1")

