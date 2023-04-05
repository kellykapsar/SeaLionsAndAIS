

# Load libraries
library(lubridate)
library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(fasterize)
library(foreach)
library(doParallel)

# Load functions
source("./S0_Functions.R")
# source("./S0_Functions.R") # Home computer path

# weekly homerange polygons
hr <- st_read("../Data_Processed/Telemetry/Homerange_KDE_weekly_20220705.shp")


# Sea lion data from 2_SSLAvailAndCovarExtraction.Rmd
# ssl5 <- readRDS(file="../Data/ssl5.rds")
load(file="../Data_Processed/Telemetry/TEMP_20220705.rda")
ssl4$season <-  getSeasonYear(ssl4$date)
ssl4$seasonhr_id <- factor(paste0(ssl4$deploy_id,"-", ssl4$season))

# Vessel tracklines
# ships <- readRDS( "../Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/AllVessels_Reprojected.rds") # Home computer path
ships <- readRDS( "../Data_Processed/AIS/AllVessels_Reprojected.rds") # HPCC path

# Land raster
land <- raster::raster("../Data_Processed/Bathymetry.tif")
land[values(land) > 0] <- 1
land[values(land) < 0] <- NA



options(mc.cores = parallel::detectCores())

# Start cluster
# cl <- 28
# clus <- makeCluster(cl)
# registerDoParallel(clus)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

##################################################################
############# CALCULATE PROXIMITY TO VESSELS ############# 
##################################################################
# Standardize dates 
ships$season <- as.Date(substr(ships$AIS_ID, 11, 18),format = "%Y%m%d") %>% getSeasonYear()

# ID unique SSL/week combos 
weeklyhr_ids <- unique(ssl4$weeklyhr_id)

# Separate fishing vessel from other vesels 
nofish <- ships[ships$AIS_Typ != "Fishing",]
fish <- ships[ships$AIS_Typ == "Fishing",]

# Afunction to create a cost distance raster tot he nearest ship using a set of points
# with their associated homerange polygon, vector lines for vessels, and a raster
# of land pixels 
extractcostdist <- function(pts, lines, land, hrs){
  # Id correct weekly homerange polygon 
  polyhr <- hrs[which(hrs$wklyhr_ == pts$weeklyhr_id[1]),]
  polyhr <- st_buffer(polyhr, 100000)
  # Id lines within the correct season
  lines <- lines[lines$season == pts$season[1],]
  # Isolate lines that spatially intersect homerange
  linesin <- lines[st_intersects(polyhr, lines, sparse=FALSE),]
  if(length(linesin$AIS_ID) == 0){return(NA)}
  # Create template raster at 250 m resolution
  template <- raster(polyhr, res=250, crs=st_crs(polyhr)) 
  values(template) <- 0
  # Convert template to sf obj
  templatesf <- stars::st_as_stars(template) %>% st_as_sf()
  # Id cells with intersecting lines in sf obj
  cellswithtraff <- st_intersects(linesin, templatesf, sparse=F)
  templatesf$traff <- colSums(cellswithtraff) > 0
  # Convert back to raster
  test <- fasterize(templatesf, template, field="traff")
  # Resample land raster to match ship raster 
  landcrop <- raster::resample(land, test, method="ngb")
  # Mask land out of ship raster 
  newtest <- raster::mask(test, landcrop, maskvalue = 1)
  # Calculate least cost path around land and convert to km 
  distshipras <- gridDistance(newtest, origin = 1, omit=NA)/1000
  # Convert back to sf object and buffer 
  distshipsf <- stars::st_as_stars(distshipras) %>% st_as_sf()
  distshipbuff <- st_buffer(distshipsf, dist=125)
  return(distshipbuff)
}

################ FISHING ################ 
# Extract pixel values from cost distance raster
start <- proc.time()

ssl5 <- foreach(i = 201:length(weeklyhr_ids), .packages = c("raster", "sf", "dplyr", "tidyr", "stars")) %dopar%{
  print(i)
  pts <- ssl4[which(ssl4$weeklyhr_id == weeklyhr_ids[i]),]
  
  # FISH
  costdistfish <- extractcostdist(pts, fish, land, hr)
  pts$prox_fish_km <- NA
  
  if(class(costdistfish) == "logical"){
    print(paste0(weeklyhr_ids[i], " removed."))
  }
  valsfish <- st_intersects(pts, costdistfish)
  
  tempfish  <- lapply(1:length(pts$easting), function(x){mean(costdistfish$layer[valsfish[[x]]])})
  pts$prox_fish_km <- unlist(tempfish)
  
  # SHIPS 
  costdist <- extractcostdist(pts, nofish, land, hr)
  pts$prox_ship_km <- NA
  
  if(class(costdist) == "logical"){
    next
  }
  vals <- st_intersects(pts, costdist)
  
  temp  <- lapply(1:length(pts$easting), function(x){mean(costdist$layer[vals[[x]]])})
  pts$prox_ship_km <- unlist(temp)
  
  saveRDS(pts, paste0("../ssl5_seasonwhr/",weeklyhr_ids[i],"_swhr.rds"))
}
