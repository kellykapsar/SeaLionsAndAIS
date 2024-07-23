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
ssl <- readRDS("../Data_Processed/Telemetry/watersealis.rds")
ssl$used <- 1

# Create an ID column that identifies each unique individual/year/week combination
ssl$weeklyhr_id <- factor(paste0(ssl$deploy_id, 
                                 substr(ssl$weekofyear, 1, 4), 
                                 substr(ssl$weekofyear, 7, 8)))

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

# Extract dynamic covariate data at used locations 
### Modified FROM AMT PACKAGE
# https://github.com/jmsigner/amt/blob/master/R/extract_covariates.R
# Also helpful: https://cran.r-project.org/web/packages/amt/vignettes/p3_rsf.html
# I modified it so that it works with weekly timescale data 
# DOES NOT work with holes in the data set any more 

extract_covar_var_time_custom <- function(xy, t, covariates) {
  
  t_covar <- raster::getZ(covariates) # Get time slices from covariates
  t_obs <- format(as.POSIXct(t), "%G-W%V") # Convert timestamp to week - same week of year as lubridate::isoweek(ssl_dates$date)
  
  wr <- sapply(t_obs, function(x) which(x == t_covar)) # Identify which slice to select
  ev <- terra::extract(covariates, xy) # Extract covariate values for all time slices at all used locations 
  cov_val <- ev[cbind(seq_along(wr), wr)] # Select only the relevant time slice for each location
  
  return(cov_val) # Return the extracted covariate values
}

