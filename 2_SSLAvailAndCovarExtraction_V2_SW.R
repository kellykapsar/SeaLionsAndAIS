# Title: SSLAvailAndCovarExtraction
# Author: Kelly Kapsar, Sydney Waloven
# Date: 2024-07-22
# Description: R script version of 2_SSLAvailAndCovarExtraction.Rmd for readability. Original script by Kelly Kapsar. Edits and updates made by Sydney Waloven.


# setup -------------------------------------------------------------------

library(sf)
library(raster)
library(scales)
library(ggmap)
library(amt) 
library(tidyverse)

# Load covariate data from source script
# source("Sydney/SourceScript_CovariateDataProcessing_V2_SW.R")

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
# trk %>% 
#   # Remove steps column before saving
#   select(-steps) %>% 
#   st_write("../Data_Processed/Telemetry/UsedLocsTrk_Clean.shp")


# Read in covariates 
# landmask <- raster("../Data_Processed/Landmask_GEBCO.tif")
# dist_500m <- raster("../Data_Processed/Dist500m.tif") %>%
#   raster::mask(landmask, maskvalue = 1)
dist_500m <- dist500m
# depth <- raster("../Data_Processed/Bathymetry.tif") %>%
#   raster::mask(landmask, maskvalue = 1)
depth <- bathy
# dist_land <- raster("../Data_Processed/DistLand.tif") %>%
#   raster::mask(landmask, maskvalue = 1)
dist_land <- distland2
# slope <- raster("../Data_Processed/slope.tif") %>%
#   raster::mask(landmask, maskvalue = 1)
# ship <- readRDS("../Data_Processed/AIS_AllOther.rds") 
ship <- shippingbrick
# fish <- readRDS("../Data_Processed/AIS_Fishing.rds")
fish <- fishingbrick
# sst <- raster("../Data_Processed/sst_weekly.tif")
# wind <- raster("../Data_Processed/wind_weekly.tif")

# Create a list of all covariate files and check resolution 
# raslist <- list(depth, dist_land, dist_500m, slope, ship, fish, sst, wind)
# rasres <- lapply(raslist, function(x) res(x)/1000)

# Remove unnecessary objects loaded from source script in the environment
rm(aisbound_sf, aisbound_sp, basemap, basemap.crop, bathy, bathy2, bathyDf,
   cargobrick, cargoras, deep, dist500m, dist500m.df, distland2, distland2.df,
   fillvalue, fishingbrick, fishingras, land, landmaskDf, otherbrick, otherras,
   shallow, shippingbrick, ships, slopeDf, sst_brick, sst_week, sst.df, sst.r,
   sst.slice, tankerbrick, tankerras, water, wind_brick, wind_week, wind.df,
   wind.r, wind.slice, cargolist, fishinglist, otherlist, shiplist, sst.array,
   start, lat, lon, studylatlon, t, t2)


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

# Remove individual rasters to save memory
rm(dist_land, dist_500m, depth, slope)

# Calculate the number of rows that share the same weeklyhr_id (used points)
pointcounts <- trk %>% 
  st_drop_geometry() %>% 
  group_by(weeklyhr_id) %>% 
  summarize(n = n()) %>%
  ungroup()


# Identifying preliminary available points --------------------------------

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

# Generate random points for each animal and combine the results. This keeps the initial data points in the same df, but keeps track of which points are available vs used by adding a "case_" column. TRUE = use. FALSE = random. This can give us a starting model to assess any collinearity between the variables.
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


# Assess collinearity -----------------------------------------------------

# Collinear predictors influence the variance estimates of the model and makes it difficult to interpret model coefficients. Determining what covariates are important requires avoiding multi-collinearity. 

# Generally correlation values over 0.70 are problematic.

cor(as.data.frame(ssl_rsf_10[ , 6:9]))

# Check the Variance Inflation Factor. This assesses the multicollinearity across ALL the predictors, with higher values indicating high collinearity of a predictor with the rest of the set. Some sources indicate a VIF over 10 is something for concern, while other indicate 3.0 should be the threshold. 

usdm::vif(as.data.frame(ssl_rsf_10[ , 6:9])) # all < 3.0


# Sensitivity analysis ----------------------------------------------------

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

# Now this is a full list which contains the model result for each run, and that is saved for each of the 500 rows of our data frame. There are 500 rows (5 different amounts of available locations, 100 model runs each)

# Result of a model run
wb.sim$result[[1]]

# Now to plot the results and make some conclusions from all this. We now "unnest" the results to plot the coefficient estimates from individual model fits with different sets of available data (1, 5, 20). The unnesting results in a much longer data frame. Run that first part of the code to see the result. After that we are recoding the "term" column, essentially changing the values in the column. We are plotting here the estimate, which is the beta coefficient from the glm regression for each covariate.

# Unnest results to plot the coefficient estimates from indv model fits with different sets of available data (1, 5, 20, 50, 100). 
p <- wb.sim %>% 
  unnest(cols = result) %>% 
  # Recode the "term" column, essentially changing the values in the column
  mutate(term = case_match(term, 
                           "(Intercept)" ~ "Intercept",
                           "dist_land" ~ "Distance to land", 
                           "dist_500m" ~ "Dist_500m",
                           "depth" ~ "Depth",
                           "slope" ~ "Slope")) %>% 
  # Plot the estimate (beta coefficient from the glm regression for each covariate)
  ggplot(aes(x = factor(frac), 
             y = estimate)) +
  geom_boxplot() + 
  facet_wrap(~ term, 
             scale = "free") +
  geom_jitter(alpha = 0.2) + 
  labs(x = "Avail. points per use location", 
       y = "Estimate") +
  theme_light()

# Save plot
ggsave("../Figures/ssl_covar_sensitivity_analysis.png", plot = p, width = 10, height = 8)

# This is a key tool in making a final assessment of the minimum number of available points needed for each use location. There is little change in the coefficient values once there are 20 points per use location, but choosing 50 would be safe.


# Random points and extracting covariates ---------------------------------

# Generate 50 available locations per use location.

# Modify function from earlier
# Function to generate random points for one animal
generate_random_points <- function(track_data) {
  
  # Calculate home range (KDE)
  hr_kde <- hr_kde(track_data, levels = c(0.95))
  
  # Generate random points within the home range
  random_pts <- random_points(hr_kde,
                              n = nrow(track_data) * 50, 
                              presence = track_data)
  
  # Add animal ID to the random points
  random_pts <- random_pts %>%
    mutate(id = unique(track_data$weeklyhr_id))
  
  return(random_pts)
}

# Split track data by animal ID
split_trk <- split(ssl_simple, ssl_simple$weeklyhr_id)


# Generate 50 random points for each animal and combine the results. This keeps the initial data points in the same df, but keeps track of which points are available vs used by adding a "case_" column. TRUE = use. FALSE = random.
ssl_rsf_50 <- map_dfr(split_trk, generate_random_points)

# Save rds
# write_rds(ssl_rsf_50, "../Data_Processed/ssl_rsf_50_random_points.rds")

# Read rds
ssl_rsf_50 <- read_rds("../Data_Processed/ssl_rsf_50_random_points.rds")

# Assign weights to available points that are very high.
ssl_rsf_50 <- ssl_rsf_50 %>% 
  mutate(weight = ifelse(case_ == T,
                         1,
                         5000))


## Extract static covariates -----
ssl_rsf_50 <- ssl_rsf_50 %>%
  extract_covariates(staticcovars) %>% 
  # Remove points that are on land
  filter(!is.na(Bathymetry))

# Make vector of static covariate names
vars <- ssl_rsf_50 %>% 
  select(DistLand:slope) %>% 
  names()

# Histograms for each static variable
hist_plots <- map(vars, ~
                    
                    ssl_rsf_50 %>% 
                    select(case_, var = .x) %>% 
                    ggplot(aes(x = var,
                               after_stat(density),
                               col = case_)) + 
                    geom_freqpoly(lwd = .7,
                                  bins = 30) + 
                    labs(title = .x))

# Put all plots into one grid. Expand to view clearly.
cowplot::plot_grid(plotlist = hist_plots) %>% ggsave("../Figures/ssl50_histplots_covar_sensitivity_analysis_20240805.png", plot = ., width = 10, height = 8)


# Boxplots for each static variable
box_plots <- map(vars, ~
                   
                   ssl_rsf_50 %>% 
                   select(case_, var = .x) %>% 
                   ggplot(aes(y = var,
                              col = case_)) + 
                   geom_boxplot() + 
                   labs(title = .x))

cowplot::plot_grid(plotlist = box_plots) %>% ggsave("../Figures/ssl50_boxplots_covar_sensitivity_analysis_20240805.png", plot = ., width = 10, height = 8)

# Save rds
# write_rds(ssl_rsf_50, "../Data_Processed/ssl_rsf_50_random_points_static.rds")

# Extract week of year as date from original data
ssl_dates <- ssl %>% 
  st_drop_geometry() %>% 
  group_by(weeklyhr_id) %>% 
  summarize(date = min(date)) %>% 
  # Add year and week of year columns
  mutate(year = lubridate::year(date),
         weekofyear = lubridate::isoweek(date))

# Join ssl_dates with ssl_rsf_50 to add the date, year, and weekofyear columns
ssl_rsf_50 <- ssl_rsf_50 %>%
  left_join(ssl_dates, 
            by = c("id" = "weeklyhr_id"))

# Save rds
# write_rds(ssl_rsf_50, "../Data_Processed/ssl_rsf_50_points_date_weeks.rds")

# Read in rds
ssl_rsf_50 <- read_rds("../Data_Processed/ssl_rsf_50_points_date_weeks.rds")


## Extract dynamic covariates -----

# Custom function
extract_covar_var_time_custom <- function(data, t, covariates) {
  
  # Pull values for coordinates
  xy <- st_coordinates(data)
  
  # Get time slices from covariates
  t_covar <- raster::getZ(covariates) # YYYY-WWW
  # Convert timestamp to week - same as weekofyear column
  t_obs <- format(data$date, "%G-W%V")
  
  # For each date in t_obs, identify which slice matches from the covariate times
  wr <- sapply(t_obs, function(x) which(x == t_covar))
  
  # Extract covariate values for all time slices at all used locations
  ev <- raster::extract(covariates, xy) 
  
  # Select only the relevant time slice for each location
  cov_val <- ev[cbind(seq_along(wr), wr)] 
  
  # Return the extracted covariate values
  return(cov_val) 
}

# Convert to sf df
ssl_rsf_50_sf <- ssl_rsf_50 %>% 
  st_as_sf(coords = c("x_", "y_"),
           crs = 32605)

## Test this function on part of the data -----

# Subset data to test
test <- ssl_rsf_50_sf[100000:100010, ]

# Test individual parts of function 
xy <- st_coordinates(test)
xy

# Get time slices from covariates
t_covar <- raster::getZ(wind) # YYYY-WWW
# Convert timestamp to week - same as weekofyear column
t_obs <- format(test$date, "%G-W%V")

# For each date in t_obs, identify which slice to select 
wr <- sapply(t_obs, function(x) which(x == t_covar))
wr

# Extract covariate values for all time slices at all used locations
ev <- raster::extract(wind, xy) 
ev

# Select only the relevant time slice for each location
cov_val <- ev[cbind(seq_along(wr), wr)] 
cov_val

# Run function on full test data
test$wind <- test %>% 
  extract_covar_var_time_custom(data = .,
                                t = date,
                                covariates = wind)

# Run function on total data for each dynamic covariate (wind, sst, ship, fish)

# Wind
ssl_rsf_50_sf$wind <- ssl_rsf_50_sf %>% 
  extract_covar_var_time_custom(data = .,
                                t = date,
                                covariates = wind)

# Sea surface temperature
ssl_rsf_50_sf$sst <- ssl_rsf_50_sf %>% 
  extract_covar_var_time_custom(data = .,
                                t = date,
                                covariates = sst)

# Other ships
ssl_rsf_50_sf$ship <- ssl_rsf_50_sf %>% 
  extract_covar_var_time_custom(data = .,
                                t = date,
                                covariates = ship)

# Fishing
ssl_rsf_50_sf$fish <- ssl_rsf_50_sf %>% 
  extract_covar_var_time_custom(data = .,
                                t = date,
                                covariates = fish)

# Convert column names to lowercase
colnames(ssl_rsf_50_sf) <- tolower(colnames(ssl_rsf_50_sf))


# Save outputs
write_rds(ssl_rsf_50_sf, "../Data_Processed/ssl_extracted_covariates_20240805.rds")
                      


