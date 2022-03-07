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
load(file="../Data/TEMP.rda")

# Vessel cost distance rasters
shiplist <- list.files("../Data/AIS/", pattern="ShippingDistance")
fishlist <- list.files("../Data/AIS/", pattern="FishingDistance")

# file.rename(paste0("../Data_Processed/AIS/",shiplist), gsub(paste0("../Data_Processed/AIS/",shiplist), pattern="_Fishing", replacement = ""))
# file.rename(paste0("../Data_Processed/AIS/",fishlist), gsub(paste0("../Data_Processed/AIS/",fishlist), pattern="_Fishing", replacement = ""))

# Boundary for AIS data acquisition
# aisbound_sf <- st_read("../Data_Raw/AIS_Bounds_FromAP/ais_reshape.shp") %>% st_transform(32605) 
# aisbound_sp <- as(aisbound_sf, "Spatial")


options(mc.cores = parallel::detectCores())

# Start cluster
# cl <- 28
# clus <- makeCluster(cl)
# registerDoParallel(clus)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

# Extract value for each used/available point as the average of the four nearest cells


# ssl4$prox_ship_km_new <- NA
# ssl4$prox_fish_km_new <- NA

newssl <- foreach(i = 1:length(unique(ssl4$date)), .packages = c("raster", "sf", "dplyr", "tidyr", "stars")) %dopar%{
  
  week <- unique(ssl4$date)[i]
  
  print(paste0("Processing ", week))
  
  sslwk <- ssl4[which(ssl4$date == week),]
  
  print(paste0(length(sslwk$date_outer), " SSL locs in ", week))
  
  fishwk <- raster(paste0("../Data/AIS/", fishlist[grepl(x = fishlist, pattern=week)]))
  
  print(paste0(week, ": Loaded fishing traffic"))
  
  shipwk <- raster(paste0("../Data/AIS/", shiplist[grepl(x = shiplist, pattern=week)]))
  
  print(paste0(week, ": Loaded shipping traffic"))
  
  distshp <- function(pt, shipshp){
    vals <- shipshp$layer[st_intersects(pt, shipshp, sparse=F)]
    ifelse(length(vals) == 0, return(NA), 
           ifelse(length(vals) > 1, return(mean(vals)), return(vals)))
  }
  
  distshipsf <- stars::st_as_stars(shipwk) %>% st_as_sf()
  distshipbuff <- st_buffer(distshipsf, dist=125)
  
  print(paste0(week, ": ship buffer complete"))
  
  prox_ship_km_new  <- lapply(sslwk$geometry, function(x){distshp(x, distshipbuff)})
  sslwk$prox_ship_km_new[which(ssl4$date == week)] <- prox_ship_km_new
  
  print(paste0(week, ": ship proximity complete"))
  
  distfishsf <- stars::st_as_stars(fishwk) %>% st_as_sf()
  distfishbuff <- st_buffer(distfishsf, dist=125)
  
  print(paste0(week, ": fish buffer complete"))
  
  prox_fish_km_new  <- lapply(sslwk$geometry, function(x){distshp(x, distfishbuff)})
  sslwk$prox_fish_km_new[which(ssl4$date == week)] <- prox_fish_km_new
  
  print(paste0(week, ": fish proximity complete"))
  
  return(sslwk)
}

saveRDS(newssl, "../Data/newssl.rds")

newssl$X <- st_coordinates(ssl4[,1])
newssl$Y <- st_coordinates(ssl4[,2])
newssldf <- st_drop_geometry(ssl4) 

write.csv(newssl, "../Data/newssl.csv")
