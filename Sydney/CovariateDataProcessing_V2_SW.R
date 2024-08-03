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
# writeRaster(bathy, "../Data_Processed/Bathymetry.tif", overwrite = TRUE)


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
# writeRaster(slope, "../Data_Processed/slope.tif", overwrite = TRUE)


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
# writeRaster(distland2, "../Data_Processed/DistLand.tif", overwrite = TRUE)


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
writeRaster(dist500m, "../Data_Processed/Dist500m.tif", overwrite = TRUE)


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

saveRDS(shippingbrick, "../Data_Processed/AIS_AllOther.rds")
writeRaster(shippingbrick, "../Data_Processed/AIS_AllOther.tif")

writeRaster(fishingbrick, "../Data_Processed/AIS_Fishing.nc",
            overwrite = TRUE, 
            format = "CDF",
            varname = "intensity", 
            varunit = "km",
            longname = "Shipping intensity (km travelled per cell) -- raster brick to netCDF",
            zname = "time", 
            xname = "lon",
            yname = "lat")

saveRDS(fishingbrick, "../Data_Processed/AIS_Fishing.rds")
writeRaster(fishingbrick, "../Data_Processed/AIS_Fishing.tif")

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
  dplyr::mutate(year = substr(AIS_ID, 11, 14))

# st_write(ships, "../Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/AllVessels_Reprojected.shp", append = FALSE)
# saveRDS(ships, "../Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/AllVessels_Reprojected.rds")


# Wind Speed --------------------------------------------------------------

# Wind speed data from the Copernicus Marine Service (Met-Op A and B satellites).

# Data set name: [GLOBAL OCEAN WIND L4 NEAR REAL TIME 6 HOURLY OBSERVATIONS](https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004). 

# Open netcdf file 
wind <- nc_open("../Data_Raw/CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_1626911972037.nc")
# Save metadata to a text file
{
  sink('../Data_Raw/CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_16269119720374.txt')
  print(wind)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(wind, "lon")
lat <- ncvar_get(wind, "lat", verbose = F)
t <- ncvar_get(wind, "time")

head(lon)

# Read in data from the wind variable and verify the dimensions of the array
wind.array <- ncvar_get(wind, "wind_speed") # 3dim array
dim(wind.array)

# Identify fill value and replace with NA
fillvalue <- ncatt_get(wind, "wind_speed", "_FillValue")
fillvalue

wind.array[wind.array == fillvalue$value] <- NA

# Close netcdf file
nc_close(wind)

# Isolate and plot a random time step to check
wind.slice <- wind.array[,,1]

dim(wind.slice) #2dim

wind.r <- raster(t(wind.slice), xmn = min(lon), xmx = max(lon), 
                 ymn = min(lat), ymx = max(lat),
                 # Found projection on the website
                 crs = CRS("+proj=longlat +datum=WGS84 +no_defs")) %>%
  flip(direction = "y")%>%
  raster::projectRaster(crs = prj)

wind.df <- as.data.frame(wind.r, xy = TRUE) %>%
  drop_na()
colnames(wind.df) <- c("x", "y", "windspeed")

plot(wind.r)


# Map of study area with wind data
ggplot() +
  geom_sf(data = basemap.crop, fill = "gray", color = "black", lwd = 0.5) +
  geom_raster(data = wind.df, 
              aes(x = x, y = y, fill = windspeed), 
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")


test <- aperm(wind.array, c(2, 1, 3)) # resize array
# Make a raster brick of all values 
wind_brick <- brick(test, xmn = min(lon), xmx = max(lon), 
                    ymn = min(lat), ymx = max(lat), 
                    crs = CRS("+proj=longlat +datum=WGS84 +no_defs")) %>% 
  flip(direction = "y") %>%
  raster::projectRaster(crs = prj) 

# Convert date from seconds since 01/01/1970 to yyyy-mm-dd format
t2 <- as.POSIXct("1900-01-01 00:00") + as.difftime(t, units = "hours")
t2 <- format(t2, "%G-W%V")

# Name raster layers after the date that they portray
names(wind_brick) <- t2
wind_brick <- raster::setZ(wind_brick, t2)

# Save output file
# writeRaster(wind_brick, "../Data_Processed/wind_AOOS_cropped_4336.tif", overwrite = T)

# Fun animation of the raster brick
# animate(wind_brick, pause=0.5, n=1)

# Calculate mean monthly wind rasters 
wind_week <- zApply(wind_brick, t2, fun = mean)

# animate(wind_week, pause=0.5, n=1)

# plot monthly average wind data for November, 2018
windweek.df <- as.data.frame(wind_week[[1]], xy = TRUE) %>%
  drop_na()
colnames(windweek.df) <- c("x", "y", "windspeed")

ggplot() +
  geom_sf(data = basemap.crop, fill = "gray", color = "black", lwd = 0.5) +
  geom_raster(data = windweek.df, 
              aes(x = x, y = y, fill = windspeed),
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")

# animate(wind_week, pause=0.5, n=1)

# writeRaster(wind_week, "../Data_Processed/wind_weekly.tif", options = "INTERLEAVE=BAND", overwrite = T)
# saveRDS(wind_week, "../Data_Processed/wind_weekly.rds")
# wind <- wind_week
# save(wind, file = "../Data_Processed/wind_weekly.rda")


# Save output file
# writeRaster(wind_week, "../Data_Processed/wind_weekly.nc",
# overwrite = TRUE, format = "CDF",
# varname = "wind_speed", varunit = "m/s",
# longname = "Wind Speed -- raster brick to netCDF",
# xname = "lon", yname = "lat", zname = "time",
# zunit = "numeric")

# Test to make sure that the netcdf file saved correctly. #

test <- nc_open("../Data_Processed/wind_weekly.nc")

# Save metadata to a text file
{
  sink('../Data_Processed/WIND_Monthly.txt')
  print(test)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(test, "lon")
lat <- ncvar_get(test, "lat", verbose = F)
t <- ncvar_get(test, "time")

head(lon)

# Read in data from the wind variable and verify the dimensions of the array
test.array <- ncvar_get(test, "wind_speed") # 3dim array
dim(test.array)

# Identify fill value and replace with NA
fillvalue <- ncatt_get(test, "wind_speed", "_FillValue")
fillvalue

test.array[test.array == fillvalue$value] <- NA

nc_close(test)

test_brick <- brick(test.array, xmn = min(lon), xmx = max(lon), 
                    ymn = min(lat), ymx = max(lat), 
                    crs = prj) 

plot(test_brick[[1]])


test.df <- as.data.frame(test_brick[[1]], xy = TRUE)
colnames(test.df) <- c("x", "y", "windspeed")

plot(test.df) # previously test.r but test.r is not called anywhere else?


# Map of study area with wind data
ggplot() +
  geom_sf(data = basemap.crop, fill = "gray", color = "black", lwd = 0.5) +
  geom_raster(data = test.df, 
              aes(x = x, y = y, fill = windspeed), 
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")


# Sea Surface Temperature -------------------------------------------------

# [OSTIA: Operational Sea Surface Temperature and Sea Ice Analysis](https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001) from the Copernicus Marine Service. 

# Variables = analysed_sst, analysis_error, mask, sea_ice_fraction
 
# North = 62
# South = 56
# West = -155
# East = -143
# Start = 01 Nov 2018
# End = 31 July 2020


# Open netcdf file 
sst <- nc_open("../Data_Raw/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1630513353420.nc")
# Save metadata to a text file
{
  sink('../Data_Raw/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1630513353420.txt')
  print(sst)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(sst, "lon")
lat <- ncvar_get(sst, "lat", verbose = F)
t <- ncvar_get(sst, "time")

head(lon)

# Read in data from the SST variable and verify the dimensions of the array
sst.array <- ncvar_get(sst, "analysed_sst") # 3dim array
dim(sst.array)

# Identify fill value and replace with NA
fillvalue <- ncatt_get(sst, "analysed_sst", "_FillValue")
fillvalue

sst.array[sst.array == fillvalue$value] <- NA

# Close netcdf file
nc_close(sst)

# Convert from kelvin to celsius
sst.array <- sst.array - 273.15

# Isolate and plot a random time step to check
sst.slice <- sst.array[,,1]

dim(sst.slice) #2dim

sst.r <- raster(t(sst.slice), xmn = min(lon), xmx = max(lon),
                ymn = min(lat), ymx = max(lat),
                # Found projection on the website
                crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")) %>%
  flip(direction = "y") %>%
  raster::projectRaster(crs = prj) %>%
  raster::crop(study)

sst.df <- as.data.frame(sst.r, xy = TRUE) %>% 
  drop_na()
colnames(sst.df) <- c("x", "y", "sst")

plot(sst.r)

# Map of study area with sst data
ggplot() +
  geom_sf(data = basemap.crop, fill = "gray", color = "black", lwd = 0.5) +
  geom_raster(data = sst.df, 
              aes(x = x, y = y, fill = sst),
              alpha = 0.9) +
  scale_fill_gradient2() +
  geom_sf(data = study, fill = NA, color = "red")

test <- aperm(sst.array, c(2, 1, 3))
# Make a raster brick of all values 
sst_brick <- brick(test, xmn = min(lon), xmx = max(lon), 
                   ymn = min(lat), ymx = max(lat), 
                   crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  flip(direction = "y") %>%
  raster::projectRaster(crs = prj) %>% 
  raster::crop(study)

# Convert date from seconds since 01/01/1970 to yyyy-mm-dd format
t2 <- as.POSIXct("1981-01-01 00:00") + as.difftime(t, units = "secs")
t2 <- format(t2, "%G-W%V")

# Name raster layers after the date that they portray
names(sst_brick) <- t2
sst_brick <- raster::setZ(sst_brick, t2)


# Calculate the change between any two layers
# tempchange <- sst_brick[[2]]-sst_brick[[1]]
# tempchange

# Get the date from the names of the layers and extract the month
sst_week <- zApply(sst_brick, by = t2, fun = mean)

# animate(sst_week, pause=0.5, n=1)

# writeRaster(sst_week, "../Data_Processed/sst_weekly.tif", options = "INTERLEAVE=BAND", overwrite = T)
# saveRDS(sst_week, "../Data_Processed/sst_weekly.rds")
sst <- sst_week
# save(sst, file = "../Data_Processed/sst_weekly.rda")

# writeRaster(sst_week, "../Data_Processed/sst_weekly.nc",
#       overwrite = TRUE, format = "CDF",
#       varname = "analysed_sst", varunit = "m/s",
#       longname = "Sea Surface Temperature -- raster brick to netCDF",
#       xname = "lon", yname = "lat", zname = "time",
#       zunit = "numeric")
