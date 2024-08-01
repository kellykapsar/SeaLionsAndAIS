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


# New function for extracting covariates at used locations ----
# Select once daily locations


# -----


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

# First set the KDE home range with 95% UD
hr_kde <- hr_kde(trk, levels = c(0.95))

ssl_rsf_10 <- random_points(hr_kde,
                            n = nrow(trk) * 10,
                            type = "regular")
# This keeps the initial data points in the same df, but keeps track of which points are available vs used by adding a "case_" column. TRUE = use. FALSE = random.
table(ssl_rsf_10$case_)

# Plot use/availability
terra::plot(ssl_rsf_10)

# Generate 5x the number of "use" points to compare against 10x
ssl_rsf_5 <- random_points(trk,
                           n = nrow(trk) * 5,
                           type = "regular")
terra::plot(ssl_rsf_5)
# There is no noticeable difference between the 10 and 5 random points per use location so let's stick with the default 10 random points.

# The argument type = "regular" samples points in a regular grid-like pattern instead of completely randomly, beneficial for spatial analysis requiring uniform coverage. The argument presence = "trk" generates points that reflect the sampling distribution or habitat patterns in the tracking data. Which should we use?




# # -----

# Create 5 random non-land points within each weekly MCP for each used point
for(i in 1:length(unique(trk$id))) {
  
  # Print iteration every 10th ID to keep track of progress
  if (i %% 10 == 0) { print(i) }
  
  # Subset df for current id
  t <- trk[which(trk$id == unique(trk$id)[i]), ]
  
  # Calculate number of random pts needed
  npts <- pointcounts$n[which(pointcounts$weeklyhr_id == unique(trk$id)[i])]*5
  
  # KDE home range with 95% UD
  hr_kde <- hr_kde(t, levels = c(0.95))
  
  # Generate random points within the KDE home range
  pts <- random_points(hr_kde, n = npts) %>%
    # Add hr_id to points
    mutate(weeklyhr_id = unique(trk$id)[i]) 
  
  # Convert to sf
  landpts <- st_as_sf(pts, coords = c("x_", "y_"), 
                      crs = prj,
                      remove = FALSE) %>% 
    # Extract static covariates
    cbind(., raster::extract(staticcovars, .)) %>%
    st_drop_geometry() %>% 
    # Remove points that are on land
    filter(!is.na(.data$Bathymetry))
  
  # Generate more random points if the required 5 points is not met yet
  while (length(landpts$case_) < npts) {
    
    newpts <- random_points(hr_kde, n = npts-length(landpts$case_)) %>%
      mutate(weeklyhr_id = unique(trk$id)[i])
    
    newlandpts <- st_as_sf(newpts, coords = c("x_", "y_"), 
                           crs = prj, 
                           remove = FALSE) %>% 
      cbind(., raster::extract(staticcovars, .)) %>%
      st_drop_geometry() %>% 
      filter(!is.na(.data$Bathymetry))
    
    landpts <- rbind(landpts, newlandpts)
  }
  
  # Convert KDE to isopleths
  hr_kde <- hr_isopleths(hr_kde) %>% 
    # Add hr_id
    mutate(weeklyhr_id = unique(trk$id)[i])
  
  # Add results to empty df
  ssl_hr <- rbind(ssl_hr, hr_kde)
  avail_pts <- rbind(avail_pts, landpts)
}

# Save output polygons 
st_write(ssl_hr, "../Data_Processed/Telemetry/Homerange_KDE_weekly_20220705.shp")


ssl_new <- ssl %>% 
  select(deploy_id, 
         weeklyhr_id,
         date,
         year, 
         weekofyear, 
         northing, 
         easting) %>% 
  cbind(., raster::extract(staticcovars, .)) %>% 
  rename(dist_land = DistLand, 
         dist_500m = Dist500m) %>%
  mutate(used = 1)

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

