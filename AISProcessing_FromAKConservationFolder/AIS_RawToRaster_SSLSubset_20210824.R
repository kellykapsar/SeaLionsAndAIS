# Revised Vessel Traffic modeling
# exactEarth CSV Data
# Created by Ben Sullender & Kelly Kapsar, 2021

# Start timer
start <- proc.time()

# Load libraries 
library(maptools)
library(rgdal)
library(dplyr)
library(tidyr)
library(tibble)
library(sf)
library(foreach)
library(doParallel)
library(spatstat)
library(raster)

####################################################################
##################### AIS PROCESSING FUNCTION ######################
####################################################################

# INPUTS: A list of lists containing all daily csv file names for one year of AIS data.  
## Inner list = file paths/names for daily AIS csvs
## Outer list = a list of months 

# OUTPUTS: 
## Monthly shapefiles for each ship type containing vectorized daily ship transit segments
## Text file with output information (number of unique ships, rows excluded, etc.)

SSL.AIS <- function(csvList, minlon, maxlon, minlat, maxlat, cellsize=1000, outproj=3338){
  # start timer 
  starttime <- proc.time()
  
  # clear variables
  AIScsv <- NA
  AIScsvDF <- NA
  AISlookup <- NA
  temp <- NA
  
  start <- proc.time()
  # read in csv (specifying only columns that we want based on position in dataframe)
  temp <-  lapply(csvList, read.csv, header=TRUE, na.strings=c("","NA"),
                  colClasses = c(rep("character", 2), "NULL", "character", "NULL", "NULL", 
                                 "character", rep("NULL", 6), "character", "NULL", rep("character",8),
                                 "NULL", "character", "NULL", "character", "NULL", rep("character", 2), rep("NULL", 109)
                  ))
  
  
  AIScsv <- do.call(rbind , temp)
  
  importtime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Create AIS_ID field
  AIScsv <- AIScsv %>% add_column(AIS_ID = paste0(AIScsv$MMSI,"-",substr(AIScsv$Time,1,8)))
  
  orig_aisids <- length(unique(AIScsv$AIS_ID))
  orig_mmsis <- length(unique(AIScsv$MMSI))
  
  
  # Convert character columns to numeric as needed
  numcols <- c(1:2, 6:12, 14:17)
  AIScsv[,numcols] <- lapply(AIScsv[,numcols], as.numeric)
  
  # Remove points outside of (a 2 degree buffered) study area 
  AIScsvDF <- AIScsv %>% filter(Latitude > minlat-2) %>% 
    filter(Latitude < maxlat+2) %>% 
    filter(Longitude > minlon-2) %>% 
    filter(Longitude < maxlon+2)
  
  inbounds_aisids <- length(unique(AIScsvDF$AIS_ID))
  inbounds_mmsis <- length(unique(AIScsvDF$MMSI))
  
  # this will come in handy later. chars 28 to 34 = "yyyy-mm"
  MoName <- substr(csvList[[1]][1],45, 51)
  yr <- substr(MoName, 1, 4) 
  mnth <- substr(MoName, 6, 7)
  print(paste0("Processing",yr, mnth))
  
  # create df
  # we only care about lookup col, time, lat + long, and non-static messages
  AIScsvDF <- AIScsvDF %>%
    dplyr::select(MMSI,Latitude,Longitude,Time,Message_ID,AIS_ID,SOG) %>%
    filter(Message_ID!=c(5,24)) %>%
    filter(!is.na(Latitude)) %>%
    filter(!is.na(Longitude)) %>%
    filter(nchar(trunc(abs(MMSI))) > 8)
  
  # Remove remaining duplicate rows of data 
  AIScsvDF <- AIScsvDF %>% distinct(.keep_all=TRUE)
  
  filt_aisids <- length(unique(AIScsvDF$AIS_ID))
  filt_mmsis <- length(unique(AIScsvDF$MMSI))
  
  dftime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Filter out points > 100 km/hr 
  AISspeed <- AIScsvDF %>%
    mutate(Time = as.POSIXct(Time, format="%Y%m%d_%H%M%OS")) %>% 
    st_as_sf(coords=c("Longitude","Latitude"),crs=4326) %>%
    # project into Alaska Albers (or other CRS that doesn't create huge gap in mid-Bering with -180W and 180E)
    st_transform(crs=outproj) %>%
    mutate(speed = NA) %>% arrange(Time)
  
  # Function from: https://www.reddit.com/r/rstats/comments/8czqni/i_have_spatiotemporal_movement_data_how_would_i/
  # Yes. I found it on reddit. Please don't judge. I'm pretty sure it works. (See SpeedFilterTesting_20210433.R for verification
  # and comparison with previous method)  
  euclidean_speed <- function(lat2, lat1, long2, long1, time2, time1) {
    latdiff <- lat2 - lat1
    longdiff <- long2 - long1
    distance <- sqrt(latdiff^2 + longdiff^2)/1000
    timediff <- as.numeric(difftime(time2,time1,units=c("hours")))
    return(distance / timediff)
  }
  
  # Calculate instantaneous speed between consecutive points 
  AISspeed[, c("long", "lat")] <- st_coordinates(AISspeed)
  AISspeed <- AISspeed %>% 
    group_by(AIS_ID) %>%
    arrange(AIS_ID, Time) %>% 
    mutate(speed = euclidean_speed(lat, lag(lat), long, lag(long), Time, lag(Time)))
  
  toofast_pts <- length(which(AISspeed$speed >= 100))
  AISspeed <- AISspeed %>% filter(speed < 100 | is.na(speed))
  
  slowenough_aisids <- length(unique(AISspeed$AIS_ID))
  slowenough_mmsis <- length(unique(AISspeed$MMSI))
  
  # Remove AIS_IDs with only one point 
  SingleAISid <- AISspeed %>% st_drop_geometry() %>% group_by(AIS_ID) %>% summarize(n=n()) %>% filter(n <= 3)
  AISspeed <- AISspeed[!(AISspeed$AIS_ID %in% SingleAISid$AIS_ID),] 
  
  enoughpts_aisids <- length(unique(AISspeed$AIS_ID))
  enoughpts_mmsis <- length(unique(AISspeed$MMSI))
  
  speedtime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # Evaluate week of year for each point
  AISspeed$week <- lubridate::week(AISspeed$Time)
  
  # Save points for potential movement animation later one
  st_write(AISspeed, paste0("../Data_Processed_SSLWeekly/Points/Points_", MoName, ".shp"))
  
  # create sf lines by sorted / grouped points
  AISsf <- AISspeed %>%
    arrange(Time) %>%
    # create 1 line per AIS ID
    group_by(AIS_ID, week) %>%
    # keep MMSI for lookup / just in case; do_union is necessary for some reason, otherwise it throws an error
    summarize(MMSI=first(MMSI),do_union=FALSE, npoints=n()) %>%
    st_cast("LINESTRING") %>% 
    st_make_valid()
  
  # Figure out which rows aren't LINESTRINGS and remove from data 
  notlines <- AISsf[which(st_geometry_type(AISsf) != "LINESTRING"),]
  AISsf <- AISsf[-which(st_geometry_type(AISsf) != "LINESTRING"),]
  
  line_aisids <- length(unique(AISsf$AIS_ID))
  line_mmsis <- length(unique(AISsf$MMSI))
  
  linetime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
  # create lookup table
  # we only care about 7 columns in total: lookup col (MMSI + date), name + IMO (in case we have duplicates / want to do an IMO-based lookup in the future),
  #     ship type, and size (in 3 cols: width, length, Draught)
  #     I'm including Destination and Country because that would be dope! We could do stuff with innocent passage if that's well populated.
  AISlookup <- AIScsv %>%
    # add_column(DimLength = AIScsv$Dimension_to_Bow+AIScsv$Dimension_to_stern, DimWidth = AIScsv$Dimension_to_port+AIScsv$Dimension_to_starboard) %>%
    dplyr::select(MMSI, Message_ID, Country, Vessel_Name, IMO, Ship_Type, Draught, Destination, Navigational_status,
                  SOG, AIS_ID, Dimension_to_Bow, Dimension_to_stern, Dimension_to_port, Dimension_to_starboard) %>%
    filter(Message_ID==c(5,24)) %>%
    filter(nchar(trunc(abs(MMSI))) > 8) %>% 
    distinct(AIS_ID, .keep_all=TRUE)
  
  
  # # Calculate number of ships with 0 values in dimensions and convert to NA
  # nolength <- length(which(is.na(AISlookup$Dimension_to_Bow | AISlookup$Dimension_to_stern)))
  # zerolength <- which(AISlookup$Dimension_to_Bow == 0 | AISlookup$Dimension_to_stern == 0)
  # pctzerolength <- round((nolength + length(zerolength))/length(AISlookup$Dimension_to_Bow)*100,2)
  # 
  # nowidth <- length(which(is.na(AISlookup$Dimension_to_port | AISlookup$Dimension_to_starboard)))
  # zerowidth <- which(AISlookup$Dimension_to_port == 0 | AISlookup$Dimension_to_starboard == 0)
  # pctzerowidth <- round((nowidth + length(zerowidth))/length(AISlookup$Dimension_to_port)*100,2)
  # 
  # # Remove zero value rows for consideration of respective measurement
  # # (i.e. if either bow or stern is zero, then both get NA 
  # # and if either port or starboard is zero then both get NA)
  # AISlookup$Dimension_to_Bow[zerolength] <- NA
  # AISlookup$Dimension_to_stern[zerolength] <- NA
  # AISlookup$Dimension_to_port[zerowidth] <- NA
  # AISlookup$Dimension_to_starboard[zerowidth] <- NA
  # 
  # # Create length and width values from dimensions 
  # AISlookup <- AISlookup %>% 
  #   add_column(DimLength = AISlookup$Dimension_to_Bow+AISlookup$Dimension_to_stern, 
  #              DimWidth = AISlookup$Dimension_to_port+AISlookup$Dimension_to_starboard) %>% 
  #   dplyr::select(MMSI, Message_ID, Country, Vessel_Name, IMO, Ship_Type, Draught, Destination, Navigational_status,
  #                 SOG, AIS_ID, DimLength, DimWidth)
  # 
  InSfNotLookup_aisids <- length(AISsf$AIS_ID[!(AISsf$AIS_ID %in% AISlookup$AIS_ID)])
  InLookupNotSf_aisids <- length(AISlookup$AIS_ID[!(AISlookup$AIS_ID %in% AISsf$AIS_ID)])
  
  # step 2: join lookup table to the lines 
  AISjoined <- AISsf %>%
    left_join(AISlookup,by="AIS_ID")
  
  # Count how many NA ship types were removed
  pctmissingshiptype <- round((length(AISjoined$AIS_ID[is.na(AISjoined$Ship_Type)])/length(AISjoined$AIS_ID))*100,2)
  
  # Remove NA ship type - mostly signals from static platforms (e.g., terrestrial AIS receivers)
  AISjoined <- AISjoined[!is.na(AISjoined$Ship_Type),]
  
  # Convert other ship type to strings 
  AISjoined$AIS_Type <- ifelse(substr(AISjoined$Ship_Type,1,1)==7, "Cargo",
                               ifelse(substr(AISjoined$Ship_Type,1,1)==8, "Tanker", 
                                      ifelse(substr(AISjoined$Ship_Type,1,2)==30, "Fishing", "Other")))
  
  
  # De spoof the data (originally from "SpoofRemover.R")
  AISjoined$length_km <- as.numeric(st_length(AISjoined)/1000)
  # Remove impossibly long lines (and count how many were removed)
  spoofremoved <- sum(AISjoined$length_km > 10000)
  AISjoined <- AISjoined %>% filter(length_km < 10000)
  
  # Number of ais_ids and mmsis by ship type
  ntank_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Tanker")]))
  ntank_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Tanker")]))
  
  nfish_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Fishing")]))
  nfish_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Fishing")]))
  
  ncargo_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Cargo")]))
  ncargo_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Cargo")]))
  
  nother_mmsis <- length(unique(AISjoined$MMSI.x[which(AISjoined$AIS_Type == "Other")]))
  nother_aisids <- length(unique(AISjoined$AIS_ID[which(AISjoined$AIS_Type == "Other")]))
  
  ntotal_mmsis <- length(unique(AISjoined$MMSI.x))
  ntotal_aisids <- length(unique(AISjoined$AIS_ID))
  
  jointime <- (proc.time() - start)[[3]]/60
  start <- proc.time()
  
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
  
  # Loop through each ship type, rasterize, and save shp file as well 
  allTypes <- unique(AISjoined$AIS_Type)
  weeks <- unique(AISjoined$week)
  
  # Set up data frame for weekly metadata collection
  weeklystats <- data.frame(week=c(), ship_type=c(), naisids=c(), nmmsis=c(), dist_km=c())
  
  # Iterate through weeks and ship types to create raster files 
  for(j in 1:length(unique(AISjoined$week))){
    for (k in 1:length(allTypes)){
      
      print(paste0("Processing: ", MoName,"-wk",weeks[j],"-",allTypes[k]))
      
      # Subset data by ship type and week of the month
      AISfilteredType <- AISjoined %>%
        filter(AIS_Type==allTypes[k] & week==weeks[j])
      
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
                           nmmsis=length(unique(AISfilteredType$MMSI.x)), 
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
      writeRaster(resRast,paste0("../Data_Processed_SSLWeekly/AISRasterSubset", MoName,"-wk",weeks[j],"-",allTypes[k],"_",cellsize,"m",".tif"))
      
      # Append metadata for study area
      weeklystats <- rbind(weeklystats, wkdata)
      
      }
  }
  
  write.csv(weeklystats, paste0("../Data_Processed_SSLWeekly/Metadata/Metadata_Weekly_",MoName,".csv"))
  # Save processing info to text file 
  rastertime <- (proc.time() - start)[[3]]/60
  
  runtime <- proc.time() - starttime 
  runtime_min <- runtime[[3]]/60 
  summarystats <- data.frame(cbind(yr, mnth, runtime_min, orig_aisids, orig_mmsis, 
                                   inbounds_aisids,inbounds_mmsis, 
                                   filt_aisids,filt_mmsis,
                                   toofast_pts, enoughpts_aisids, enoughpts_mmsis,
                                   slowenough_aisids, slowenough_mmsis,
                                   line_aisids, line_mmsis, 
                                   InSfNotLookup_aisids, InLookupNotSf_aisids,
                                   pctmissingshiptype, 
                                   ntank_aisids, ntank_mmsis,
                                   nfish_aisids, nfish_mmsis,
                                   ncargo_aisids, ncargo_mmsis,
                                   nother_aisids, nother_mmsis,
                                   ntotal_mmsis, ntotal_aisids))
  write.csv(summarystats, paste0("../Data_Processed_SSLWeekly/Metadata/Metadata_",MoName,".csv"))
  
  
  runtimes <- data.frame(cbind(yr, mnth, runtime_min, importtime, dftime, speedtime, linetime, jointime, rastertime)) 
  write.csv(runtimes, paste0("../Data_Processed_SSLWeekly/Metadata/Runtimes_",MoName,".csv"))
  print(runtimes)
  return(runtimes)
}

################################################
# Set parameters
minlat <- 56
maxlat <- 62
minlon <- -155
maxlon <- -143
cellsize <- 1000
outproj <- 3338




#########################################################
####################### TEST CODE ####################### 
#########################################################
# Pull up list of AIS files
# files <- paste0("../Data_Raw/2020/", list.files("../Data_Raw/2020", pattern='.csv'))
# 
# jul <- files[grepl("-07-", files)]
# csvList <- jul[c(1,10,21)]
# 

# SSL.AIS(csvList, minlon, maxlon, minlat, maxlat, cellsize, outproj)

####################################################################
####################### PARALLELIZATION CODE ####################### 
####################################################################

# Pull up list of AIS files
files <- paste0("../Data_Raw/2018/", list.files("../Data_Raw/2018", pattern='.csv'))

# Separate file names into monthly lists
# jan <- files[grepl("-01-", files)]
# feb <- files[grepl("-02-", files)]
# mar <- files[grepl("-03-", files)]
# apr <- files[grepl("-04-", files)]
# may <- files[grepl("-05-", files)]
# jun <- files[grepl("-06-", files)]
# jul <- files[grepl("-07-", files)]
# aug <- files[grepl("-08-", files)]
# sep <- files[grepl("-09-", files)]
# oct <- files[grepl("-10-", files)]
nov <- files[grepl("-11-", files)]
dec <- files[grepl("-12-", files)]

# Create a list of lists of all csv file names grouped by month
# csvsByMonth <- list(jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec)
csvsByMonth <- list(nov, dec)



## MSU HPCC: https://wiki.hpcc.msu.edu/display/ITH/R+workshop+tutorial#Rworkshoptutorial-Submittingparalleljobstotheclusterusing{doParallel}:singlenode,multiplecores
# Request a single node (this uses the "multicore" functionality)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

# create a blank list to store the results (I truncated the code before the ship-type coding, and just returned the sf of all that day's tracks so I didn't 
#       have to debug the raster part. If we're writing all results within the function - as written here and as I think we should do - the format of the blank list won't really matter.)
res=list()

# foreach and %dopar% work together to implement the parallelization
# note that you have to tell each core what packages you need (another reason to minimize library use), so it can pull those over
# I'm using tidyverse since it combines dplyr and tidyr into one library (I think)
res=foreach(i=1:2,.packages=c("maptools", "rgdal", "dplyr", "tidyr", "tibble", "stars", "raster", "foreach", "doParallel", "spatstat", "raster"),
            .errorhandling='pass',.verbose=T,.multicombine=TRUE) %dopar% SSL.AIS(csvList=csvsByMonth[[i]],
                                                                                 minlon, maxlon, minlat, maxlat, cellsize, outproj)
# lapply(csvsByMonth, FWS.AIS)

# Elapsed time and running information
tottime <- proc.time() - start
tottime_min <- tottime[[3]]/60

cat("Time elapsed:", tottime_min, "\n")
cat("Currently registered backend:", getDoParName(), "\n")
cat("Number of workers used:", getDoParWorkers(), "\n")
