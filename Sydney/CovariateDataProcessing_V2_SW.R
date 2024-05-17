# Title: Covariate Data Processing V2
# Author: Kelly Kapsar
# Date: 2024-05-17
# Description: R script version of 1b_CovariateDataProcessing.Rmd for readability


# setup -------------------------------------------------------------------

# Import libraries
library(sf)
library(raster)
library(marmap)
library(ncdf4)
library(anytime)
library(rgdal)
library(scales)
library(tidyverse)


# Study area boundaries ---------------------------------------------------

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
# st_write(study, "../Data_Raw/studyarea.shp")

basemap <- read_sf("../Data_Raw/BBmap.shp") %>%
  st_transform(prj) %>%  
  st_buffer(0)

# Crop basemap to buffered extent of study area 
study.buff <- st_buffer(study, 100000) # Buffer study area by 100 km
basemap.crop <- st_crop(basemap, study.buff)

# Map of study area 
ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", 
          color = "black", 
          lwd = 0.5) +
  geom_sf(data = study,
          fill = NA, 
          color = "red")

# Get latlong coordinates for study area for use in downloading other data sets
studylatlon <- study %>% 
  st_transform(4269) %>% 
  st_bbox()


# Bathymetry --------------------------------------------------------------

# Bathymetry data were collected from the [General Bathymetric Chart of the Oceans](gebco.net). You can use the marmap package to load the data and create a "bathy" object type, but I'm not sure exactly how to work with that kind of object, so I just imported the tif data and converted into a raster. 

# bathymarmap <- marmap::getNOAA.bathy(lon1 = -155, lon2 = -143, lat = 56, lat2 = 62, resolution = 3) 
# # resolution units are in minutes (3 minutes = 0.05 degrees)

# Import raster, reproject to 4336, and crop to study area
bathy <- raster::raster("../Data_Raw/GEBCO_16_Aug_2021/gebco_2021_n64.0_s54.0_w-157.0_e-141.0.tif") %>%
  projectRaster(crs = prj) %>%
  raster::crop(study)

bathyDf <- as.data.frame(bathy, xy = TRUE) %>%
  drop_na()

colnames(bathyDf) <- c("x", "y", "depth")


# Map of study area with bathymetric data 
ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", 
          color = "black",
          lwd = 0.5) +
  geom_raster(data = bathyDf, 
              aes(x = x, y = y, fill = depth), 
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")

# Save raster object
# writeRaster(bathy, "../Data_Processed/Bathymetry.tif")


# Landmask ----------------------------------------------------------------

# This landmask is at the same scale as the bathymetry data. Note: Zero values were included as land and not water.

landmask <- bathy
landmask[landmask > 0 | landmask == 0] <- 1
landmask[landmask < 0] <- NA

landmaskDf <- as.data.frame(landmask, xy = TRUE) %>%
  drop_na()

colnames(landmaskDf) <- c("x", "y", "land")

ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", 
          color = "black", 
          lwd = 0.5) +
  geom_tile(data = landmaskDf, 
            aes(x = x, y = y, fill = land), 
            alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")

# writeRaster(landmask, "../Data_Processed/Landmask_GEBCO.tif", overwrite = TRUE)


# Slope -------------------------------------------------------------------

slope <- terrain(bathy, opt = "slope", unit = "radians", neighbors = 8)

slopeDf <- as.data.frame(slope, xy = TRUE) %>% 
  drop_na()

colnames(slopeDf) <- c("x", "y", "slope")

# Map of study area with bathymetric data 
ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", 
          color = "black", 
          lwd = 0.5) +
  geom_raster(data = slopeDf, 
              aes(x = x, y = y, fill = slope), 
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")

# saveRDS(slope, "../Data_Processed/slope.rds")
# writeRaster(slope, "../Data_Processed/slope.tif")


# Distance to land --------------------------------------------------------

# Calculate distance to nearest NA cell for all cells in landmask raster
distland2 <- bathy
values(distland2)[values(distland2) <= 0] = NA # set all non-positive values in distland2 to "NA" 

# ID cells on land and water
land <- distland2 > 0
water <- distland2 <= 0

# Land = 1, isobath = 2
distland2[water] <- NA
distland2[land] <- 2

# Calculate distance to land 
distland2 <- gridDistance(distland2, origin = 2)/1000

# Plot results
distland2.df <- as.data.frame(distland2, xy = TRUE) %>% 
  drop_na()

colnames(distland2.df) <- c("x", "y", "distland")

ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray",
          color = "black",
          lwd = 0.5) +
  geom_raster(data = distland2.df, 
              aes(x = x, y = y, fill = distland), 
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")

# Save output
# writeRaster(distland2, "../Data_Processed/DistLand.tif")


# Distance to shelf break -------------------------------------------------

# Calculate distance to nearest cell between 400 and 500 m depth.

# Based on this tutorial: 
# https://www.r-bloggers.com/2020/02/three-ways-to-calculate-distances-in-r/

bathy2 <- bathy

# ID cells on land, shallower than isobath, and deeper than isobath
land <- bathy2 > 0
shallow <- bathy2 > -500
deep <- bathy2 < -510 

# Land = 1, isobath = 2
bathy2[!shallow & !deep] <- 2
bathy2[raster::coordinates(bathy2)[,2] > 6650000] <- 0 # This removes the deeper parts of PWS 
bathy2[land] <- 1

# Calculate distance to nearest cell with a value of 2 (going around any cell with a value of 1)
dist500m <- gridDistance(bathy2, 
                         origin = 2,
                         omit = 1)/1000

# Crop to study area 
dist500m <- dist500m # %>% raster::crop(study) %>% raster::mask(study)

# Save results
# writeRaster(dist500m, "../Data_Processed/Dist500m.tif")


# Plot results
dist500m.df <- as.data.frame(dist500m, xy = TRUE) %>% 
  drop_na()

colnames(dist500m.df) <- c("x", "y", "dist500")

ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", color = "black", lwd = 0.5) +
  geom_raster(data = dist500m.df, 
              aes(x = x, y = y, fill = dist500), 
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")


# AIS Data - Intensity ----------------------------------------------------

start <- proc.time()

# Boundary for AIS data acquisition
aisbound_sf <- st_read("../Data_Raw/AIS_Bounds_FromAP/ais_reshape.shp") %>% 
  st_transform(prj) %>% 
  st_crop(study)

aisbound_sp <- aisbound_sf %>% as("Spatial")

# Import raster, reproject to 4336, and crop to study area
fishinglist <- list.files("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", 
                          pattern = "Fishing")
otherlist <- list.files("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", 
                        pattern = "Other")
cargolist <- list.files("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/",
                        pattern = "Cargo")
tankerlist <- list.files("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", 
                         pattern = "Tanker")

# Create raster bricks for each ship type -----

fishingras <- lapply(fishinglist, 
                     # Create a raster by combining file path with current list item
                     function(x) {projectRaster(raster::raster(
                       paste0("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", x)), 
                       crs = prj)} )

# Combine the list of raster objects into one raster brick
fishingbrick <- raster::brick(fishingras) %>%
  # Mask the raster brick using 'aisbound_sp'
  raster::mask(aisbound_sp)

otherras <- lapply(otherlist, 
                   function(x) {projectRaster(raster::raster(
                     paste0("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", x)),
                     crs = prj)} )
otherbrick <- raster::brick(otherras) %>% 
  raster::mask(aisbound_sp)

cargoras <- lapply(cargolist, 
                   function(x) {projectRaster(raster::raster(
                     paste0("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", x)), 
                     crs = prj)} )
cargobrick <- raster::brick(cargoras) %>%
  raster::mask(aisbound_sp)

tankerras <- lapply(tankerlist, 
                    function(x) {projectRaster(raster::raster(
                      paste0("../Data_Raw/AIS_SSLWeeklySubset/Raster/5000m/", x)), 
                      crs = prj)})
tankerbrick <- raster::brick(tankerras) %>% 
  raster::mask(aisbound_sp)

shippingbrick <- otherbrick + cargobrick + tankerbrick

# Convert units to km 
shippingbrick <- shippingbrick/1000
fishingbrick <- fishingbrick/1000

# Create list of weeks 
yearmon <- seq(as.Date("2018/11/1"), as.Date("2020/7/31"), "week")
yearmon <- format(yearmon, "%G-W%V")

names(shippingbrick) <- yearmon
names(fishingbrick) <- yearmon

shippingbrick <- raster::setZ(shippingbrick, yearmon)
fishingbrick <- raster::setZ(fishingbrick, yearmon)


# Save raster bricks as netcdf files 

# Save output file
writeRaster(shippingbrick, "../Data_Processed/AIS_AllOther.nc",
            overwrite = TRUE, 
            format = "CDF",
            varname = "intensity", 
            varunit = "km",
            longname = "Shipping intensity (km travelled per cell) -- raster brick to netCDF",
            xname = "lon", 
            yname = "lat",
            zname = "time",
            zunit = "numeric")

# saveRDS(shippingbrick, "../Data_Processed/AIS_AllOther.rds")

writeRaster(fishingbrick, "../Data_Processed/AIS_Fishing.nc",
            overwrite = TRUE, 
            format = "CDF",
            varname = "intensity", 
            varunit = "km",
            longname = "Shipping intensity (km travelled per cell) -- raster brick to netCDF",
            zname = "time", 
            xname = "lon",
            yname = "lat")

# saveRDS(fishingbrick, "../Data_Processed/AIS_Fishing.rds")

# Plot results
fish.df <- as.data.frame(fishingbrick[[1]], xy = TRUE) %>% 
  drop_na()
colnames(fish.df) <- c("x", "y", "intensity")

# Fun animation of the raster brick - may cause R to crash
# animate(fishingbrick, pause = 0.5, n = 1)

# Map of study area with bathymetric data 
ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", 
          color = "black", 
          lwd = 0.5) +
  geom_raster(data = fish.df, 
              aes(x = x, y = y, fill = intensity),
              alpha = 0.9) +
  # geom_sf(data = aisbound_sf, fill = NA, color = "blue" ) +
  geom_sf(data = study, fill = NA, color = "red") +
  scale_fill_gradient2(trans = "log10", low = "black", mid = "gray", high = "red")



# AIS Data - Proximity ----------------------------------------------------

# Import ship data, and generate year value
# These data were generated by the 0_AIS_Vectorization_SSLSubset_20210824.R script
shiplist <- list.files("../Data_Raw/AIS_SSLWeeklySubset/Vector", 
                      pattern = ".shp")

# Create sf objects
ships <- lapply(shiplist, 
                function(x) {st_read(
                  paste0("../Data_Raw/AIS_SSLWeeklySubset/Vector/", x)) %>% 
                    st_transform(prj)})

# Combine all ships into one sf object
ships <- do.call(rbind, ships)

# Extract year from AIS_ID
ships <- ships %>% 
  mutate(year = substr(AIS_ID, 11, 14))

st_write(ships, "../Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/AllVessels_Reprojected.shp")
saveRDS(ships, "../Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/AllVessels_Reprojected.rds")

