library(spatstat)
library(rgdal)
library(raster)
library(maptools)
library(sp)
library(sf)
library(dplyr)

# outProj ( = projection you want to end up with) must be in proj4string format
# # (also it sounds like proj4string format is getting cancelled by the sf team, BTW - we may soon have to learn yet another coordinate system format)


AIS.Rasta <- function(vectorName, dsn, savedsn, studydsn, cellsize=25000, outProj=AA){
  # start timer 
  starttime <- proc.time()
  
  study <- readOGR(studydsn) 
  studyProj <- study@proj4string

  #read line shapefile
  moSHP <- readOGR(dsn, vectorName)

  # define AA coordinate system (this is the lazy way, otherwise we could copy/paste from spatialreference.org)
  AA <- moSHP@proj4string
  
  # convert AOI to AA
  studyAA <- spTransform(study,AA)
  
  #convert to spatial lines format - this step takes the longest
  moPSP <- as.psp(moSHP)
  
  #create bounding extent for all area
  extentAOI <- as.owin(list(xrange=c(studyAA@bbox[1,1],studyAA@bbox[1,2]),yrange=c(studyAA@bbox[2,1],studyAA@bbox[2,2])))
  
  # clip vessel tracks to custom AOI; faster and safer (with geometry) than gIntersect, but uses bounding box so it's a big rectangle regardless of input geometry.
  # this is technically optional, but I thought it might help expedite things if AOI is small and cellsize is tiny
  clippedMoPSP <- clip.psp(moPSP,extentAOI)
  
  #create mask with pixel size
  allMask <- as.mask(extentAOI,eps=cellsize)
  
  #run pixellate with mask
  moPXL <- pixellate.psp(clippedMoPSP,W=allMask)
  
  #render as raster
  moRAST <- raster(moPXL)
  # we forgot to set spatial reference for raster before; I didn't realize it gets automatically wiped out
  crs(moRAST) <- AA
  
  # if line 56 (rasterize) gives you trouble, try uncommenting line 55 and switching "moRast" to "blankR" in line 56
  #  blankR <- raster(extent(studyAA),resolution=cellsize)
  AOIr <- rasterize(studyAA,moRAST,1)
  crs(AOIr) <- AA
  # the above steps create a raster with values of 1 everywhere there was an AOI polygon filled in and NAs everywhere else
  
  # so when you multiply it by the vessel density raster, you just get results where you want
  resRast <- moRAST*AOIr
  
  # in case you want to switch coordinate systems; but I'd recommend doing this post hoc, one result raster at a time. 
  # I think we want to deliver all data, both vector and raster, in Alaska Albers for simplicity's sake
  moRAST <- projectRaster(moRAST,res=cellsize,crs=outProj)
  #write output
  writeRaster(resRast,paste0(savedsn,"/AISRasterSubset",substr(vectorName,7,nchar(vectorName)),"_",cellsize,"m",".tif"))
  
  runtime <- proc.time() - starttime 
  print(runtime)
}



dsnstudy <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Data_Raw"
studyarea <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Data_Raw/studyarea.shp"

# Setup to see syntax - yours will vary based on input/output file structure
files2018 <- list.files("D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/2018/Vector/", pattern='.shp')[41:48]
files2019 <- list.files("D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/2019/Vector/", pattern='.shp')
files2020 <- list.files("D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/2020/Vector/", pattern='.shp')[1:28]

dsn2018 <- "D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/2018/Vector"
dsn2019 <- "D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/2019/Vector"
dsn2020 <- "D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/2020/Vector"

# the "-4" part gets rid of ".shp"
inS2018 <- substr(files2018,1,nchar(files2018)-4)
inS2019 <- substr(files2019,1,nchar(files2019)-4)
inS2020 <- substr(files2020,1,nchar(files2020)-4)

savedsn <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/SeaLionsAndAIS/Data_Processed/AIS"

utm5 <- 32605 # WGS84 UTM Zone 5
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

start <- proc.time()
lapply(inS2018, function(x){AIS.Rasta(vectorName=x, dsn=dsn2018, savedsn=savedsn, cellsize=1000, studydsn=studyarea, outProj=utm5)})
lapply(inS2019, function(x){AIS.Rasta(vectorName=x, dsn=dsn2019, savedsn=savedsn, cellsize=1000, studydsn=studyarea, outProj=utm5)})
lapply(inS2020, function(x){AIS.Rasta(vectorName=x, dsn=dsn2020, savedsn=savedsn, cellsize=1000, studydsn=studyarea, outProj=utm5)})
tottime <- proc.time()-start

# TEST RUN ON ONE MONTH. 
# AIS.Rasta(vectorName=inS2018[[1]], dsn=dsn2018, savedsn=savedsn, studydsn=studyarea, cellsize=1000, outProj=utm5)


######################################################################
# sample implementation:
setwd("/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Data/Unzipped/Sandbox")
files <- list.files(pattern='.shp')
inS <- substr(files,1,nchar(files)-4)
dsn1 <- "/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Data/Unzipped/Sandbox"
savedsn1 <- "/Users/bensullender/Documents/KickstepApproaches_Projects/FWS_AIS/Data/Unzipped/Sandbox"
AK5km <- c("/Users/bensullender/Documents/GIS/Basedata/AK_5km_buffer.shp")


AIS.Rasta(inS,dsn1,savedsn1,AK5km,cellsize=1000)
######################################################################

