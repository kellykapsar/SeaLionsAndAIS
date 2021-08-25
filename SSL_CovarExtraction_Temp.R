
library(tidyr)
library(dplyr)
library(sf)
library(raster)


# Read in covariates 
landmask <- raster("../Data_Processed/Landmask_GEBCO.tif")
dist500m <- raster("../Data_Processed/Dist500m.tif") %>% raster::mask(landmask, maskvalue=1)
depth <- raster("../Data_Processed/Bathymetry.tif") %>% raster::mask(landmask, maskvalue=1)
distland <- raster("../Data_Processed/DistLand.tif") %>% raster::mask(landmask, maskvalue=1)
covarstack <- raster::stack(distland, dist500m, depth)


used <- readRDS("../Data_Processed/watersealis.rds")
usedbuff <- readRDS("../Data_Processed/watersealibuffs.rds")

# ID available points based on previously calculated radii
# Accounts for points on land 
custom_extract_avail <- function(usedbuff, covarstack, npts = 5){
  df <-  as.data.frame(sampleRandom(x=covarstack, size = npts, na.rm = TRUE, ext = as(usedbuff, "Spatial"), xy = TRUE))
  df$Used <- 0
  df$Date <- usedbuff$Date
  df$DeployID <- usedbuff$DeployID
  return(df)
}

# Runing function as a loop takes ~15 hrs.... 
start <- proc.time()
q <- lapply(1:length(usedbuff$DeployID), 
            function(x){custom_extract_avail(usedbuff[x,], covarstack=covarstack)})
avail <- do.call(rbind, q)
proc.time()-start

writeRDS(avail, "../Data_Processed/AvailableLocsWithCovars.rds")