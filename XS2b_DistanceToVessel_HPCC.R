library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(fasterize)
library(foreach)
library(doParallel)


# time begin
print("Time begin"); print(start); cat("\n"); 


# Sea lion data from 2_SSLAvailAndCovarExtraction.Rmd
# load(file="../Data/TEMP_seasonal.rda")
load(file="../Data_Processed/Telemetry/TEMP_seasonal.rda")

# Vessel cost distance rasters
# shiplist <- list.files("../Data/AIS/", pattern="SeasonalShipDist")
# fishlist <- list.files("../Data/AIS/", pattern="SeasonalFishDist")

shiplist <- list.files("../Data_Processed/AIS/", pattern="SeasonalShipDist")
fishlist <- list.files("../Data_Processed/AIS/", pattern="SeasonalFishDist")


options(mc.cores = parallel::detectCores())

# Start cluster
# cl <- 28
# clus <- makeCluster(cl)
# registerDoParallel(clus)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

# Extract value for each used/available point as the average of the four nearest cells


# ssl4$prox_ship_km_new <- NA
# ssl4$prox_fish_km_new <- NA

newssl <- foreach(i = 1:length(unique(ssl4$season)), .packages = c("raster", "sf", "dplyr", "tidyr", "stars")) %dopar%{
  
  season <- unique(ssl4$season)[i]
  
  print(paste0("Processing ", season))
  
  sslssn <- ssl4[which(ssl4$season == season),]
  
  print(paste0(length(sslssn$ship), " SSL locs in ", season))
  
  fishwk <- raster(paste0("../Data_Processed/AIS/", fishlist[grepl(x = fishlist, pattern=season)]))
  
  print(paste0(season, ": Loaded fishing traffic"))
  
  shipwk <- raster(paste0("../Data_Processed/AIS/", shiplist[grepl(x = shiplist, pattern=season)]))
  
  print(paste0(season, ": Loaded shipping traffic"))
  
  distshp <- function(pt, shipshp){
    vals <- shipshp$layer[st_intersects(pt, shipshp, sparse=F)]
    ifelse(length(vals) == 0, return(NA), 
           ifelse(length(vals) > 1, return(mean(vals)), return(vals)))
  }
  
  distshipsf <- stars::st_as_stars(shipwk) %>% st_as_sf()
  distshipbuff <- st_buffer(distshipsf, dist=125)
  
  print(paste0(season, ": ship buffer complete"))
  
  prox_ship_km_new  <- lapply(sslssn$geometry, function(x){distshp(x, distshipbuff)})
  sslssn$prox_ship_km_new[which(ssl4$season == season)] <- prox_ship_km_new
  
  print(paste0(season, ": ship proximity complete"))
  
  distfishsf <- stars::st_as_stars(fishwk) %>% st_as_sf()
  distfishbuff <- st_buffer(distfishsf, dist=125)
  
  print(paste0(season, ": fish buffer complete"))
  
  prox_fish_km_new  <- lapply(sslssn$geometry, function(x){distshp(x, distfishbuff)})
  sslssn$prox_fish_km_new[which(ssl4$season == season)] <- prox_fish_km_new
  
  print(paste0(season, ": fish proximity complete"))
  
  return(sslssn)
}

saveRDS(newssl, "../Data/newssl.rds")

newssl$X <- st_coordinates(ssl4[,1])
newssl$Y <- st_coordinates(ssl4[,2])
newssldf <- st_drop_geometry(ssl4) 

write.csv(newssl, "../Data/newssl.csv")
