
# Load Libraries
library(spatstat)
library(rgdal)
library(raster)
library(maptools)
library(sf)
library(sp)
library(dplyr)

# Define rasterization function 
SSL.Rasta <- function(vectorName, dsn, savedsn, cellsize=1000, outproj=3338){
  # Convert study area to spatial polygon 
  coords <- data.frame(lat=c(56, 62, 62, 56, 56), lon=c(-155, -155, -143, -143, -155), id="study")
  studyAA <- coords %>% 
    st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
    group_by(id) %>% 
    summarize(geometry = st_combine(geometry)) %>%  
    st_cast("POLYGON") %>% 
    st_transform(outproj) %>% as("Spatial")
  
  #create bounding extent for all area
  extentAOI <- as.owin(list(xrange=c(studyAA@bbox[1,1],studyAA@bbox[1,2]),yrange=c(studyAA@bbox[2,1],studyAA@bbox[2,2])))
  
  #read line shapefile
  # AISjoined <- st_read(paste0(dsn, vectorName, ".shp"))
  AISjoined <- vectorName
  
  # Loop through each ship type, rasterize, and save shp file as well 
  allTypes <- unique(AISjoined$AIS_Typ)
  weeks <- unique(AISjoined$week)
  
  # Set up data frame for weekly metadata collection
  weeklystats <- data.frame(week=c(), ship_type=c(), naisids=c(), nmmsis=c(), dist_km=c())
  
  # Iterate through weeks and ship types to create raster files 
  for(j in 1:length(unique(AISjoined$week))){
    for (k in 1:length(allTypes)){
      
      print(paste0("Processing: ",weeks[j],"-",allTypes[k]))
      
      # Subset data by ship type and week of the month
      AISfilteredType <- AISjoined %>%
        filter(AIS_Typ==allTypes[k] & week==weeks[j])
      
      # Check to make sure that there were ships there. If not, then fill data frame with 0s 
      if(length(AISfilteredType$AIS_ID) == 0){
        wkdata <- data.frame(week=j, 
                             ship_type=allTypes[k], 
                             naisids=0, 
                             nmmsis=0, 
                             dist_km=0)
        weeklystats <- rbind(weeklystats, wkdata)
        next
      }
      
      # Fill data frame with attributes
      wkdata <- data.frame(week=j, 
                           ship_type=allTypes[k], 
                           naisids=length(unique(AISfilteredType$AIS_ID)), 
                           nmmsis=length(unique(AISfilteredType$MMSI_x)), 
                           dist_km=sum(st_length(AISfilteredType))/1000)
      
      # convert to sp object
      moSHP <- as(AISfilteredType, "Spatial")
      
      #convert to spatial lines format - this step takes the longest
      moPSP <- as.psp(moSHP)
      
      # clip vessel tracks to custom AOI; faster and safer (with geometry) than gIntersect, but uses bounding box so it's a big rectangle regardless of input geometry.
      # this is technically optional, but I thought it might help expedite things if AOI is small and cellsize is tiny
      clippedMoPSP <- clip.psp(moPSP,extentAOI)
      
      #create mask with pixel size
      allMask <- as.mask(extentAOI,eps=cellsize) # Cell size of 1 km
      
      #run pixellate with mask
      moPXL <- pixellate.psp(clippedMoPSP,W=allMask)
      
      #render as raster
      moRAST <- raster(moPXL)
      # we forgot to set spatial reference for raster before; I didn't realize it gets automatically wiped out
      crs(moRAST) <- outproj
      
      # if line 56 (rasterize) gives you trouble, try uncommenting line 55 and switching "moRast" to "blankR" in line 56
      #  blankR <- raster(extent(studyAA),resolution=cellsize)
      AOIr <- rasterize(studyAA,moRAST,1)
      crs(AOIr) <- outproj
      # the above steps create a raster with values of 1 everywhere there was an AOI polygon filled in and NAs everywhere else
      
      # so when you multiply it by the vessel density raster, you just get results where you want
      resRast <- moRAST*AOIr
      
      #write output
      writeRaster(resRast,paste0(savedsn, "Raster/Raster", weeks[j],"-",allTypes[k],"_",cellsize,"m",".tif"))

      # Append metadata for study area
      weeklystats <- rbind(weeklystats, wkdata)
  }
}

write.csv(weeklystats, paste0(savedsn,"Metadata/Metadata.csv"))
}



wd <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/"
# Setup to see syntax - yours will vary based on input/output file structure
files <- list.files(paste0(wd, "Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/"), pattern='.shp')

# the "-4" part gets rid of ".shp"
inS <- substr(files,1,nchar(files)-4)
dsn <- paste0(wd, "Data_Raw/AIS_SSLWeeklySubset/Vector/")
savedsn <- paste0(wd, "Data_Raw/AIS_SSLWeeklySubset/")

# temp <- lapply(files, function(x)st_read(paste0(dsn, x)))
# t <- do.call(rbind, temp) 
t <- readRDS("../Data_Raw/AIS_SSLWeeklySubset/Vector/EPSG32605/AllVessels_Reprojected.rds")

cellsize <- 250
outproj <- 32605

t <- st_transform(t, outproj)

start <- proc.time()
# lapply(inS, function(x){SSL.Rasta(vectorName=x, dsn=dsn, savedsn=savedsn, cellsize=5000, outproj=3338)}) 
SSL.Rasta(vectorName=t, dsn=dsn, savedsn=savedsn, cellsize=250, outproj=32605)
tottime <- proc.time()-start
browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")








