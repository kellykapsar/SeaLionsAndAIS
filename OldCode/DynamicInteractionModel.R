# Calculating coefficient of association between fishing vessels and sea lions

# Based on this vignette: 
# https://cran.r-project.org/web/packages/wildlifeDI/vignettes/wildlifeDI-vignette.html

# Load libraries
library(sf)
library(wildlifeDI)
library(adehabitatLT)
library(tidyr)
library(dplyr)
library(ggplot2)
library(xts)

# Location of study area shape file on my computer
study <- read_sf("../Data_Raw/DraftSeaLionAOI_Envelope_20kmBuff_6334.shp")


# Read in SSL data 
## Sea Lion location geodatabase
# Import sea lion location geodatabase (downloaded from Google Drive folder)
seali <- st_layers("../Data_Raw/SSL Adult Female Analysis 2018-20.gdb")
seali$name # list of layers 

# Determine all layers of interest within gdb 
lyrs <- grep("_loc", seali$name)

# Create initial sf object with data from one sea lion
sealis <- st_read("../Data_Raw/SSL Adult Female Analysis 2018-20.gdb",
                  layer=seali$name[lyrs[1]])

# Remove that layer from the layers of interest list
lyrs <- lyrs[2:length(lyrs)]

# Append all other layers of interest onto the main location data set 
for(i in 1:length(lyrs)){
  temp <- st_read("../Data_Raw/SSL Adult Female Analysis 2018-20.gdb",
                  layer=seali$name[lyrs[i]])
  colnames(temp) <- colnames(sealis)
  sealis <- rbind(sealis, temp)
}

# Adjust duplicate time stamps by 10 seconds 
# Recommended by https://jmlondon.github.io/crawl-workshop/crawl-practical.html#duplicate-times
make_unique <- function(x) {
  xts::make.time.unique(x$Date,eps = 10)
}

sealis <-  sealis  %>% 
  st_transform(crs=CRS("+proj=utm +zone=5 +datum=WGS84 +units=m +no_defs")) %>% 
  mutate(Date = as.POSIXct(Date), 
         Longitude = st_coordinates(.)[,1], 
         Latitude = st_coordinates(.)[,2]) %>% 
  dplyr::arrange(DeployID, Date) %>% 
  dplyr::group_by(DeployID) %>% tidyr::nest() %>% 
  dplyr::mutate(unique_time = purrr::map(data, make_unique)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-Date) %>% rename(Date = unique_time)

# Create trajectory

SSLtraj <- as.ltraj(data.frame(y = sealis$Latitude, x = sealis$Longitude), sealis$Date,
                    id = sealis$DeployID, typeII = TRUE)
SSL774 <- SSLtraj[1]
head(SSL774[[1]])
SSL775 <- SSLtraj[2]
# Temporal overlap 
checkTO(SSL774, SSL775)
# Identify simulatenous fixes  
sim774775 <- GetSimultaneous(SSL774, SSL775, tc=60*60) # Units for tc = seconds
SSL774.sim <- sim774775[1]
SSL775.sim <- sim774775[2]
# Proximity analysis
Prox(SSL774, SSL775, tc=60*60, dc=1000)
# 0.003 -- 0.3% of fixes were within one hour and one km of each other

# Coefficient of association 
Ca(SSL774, SSL775, tc=60*60, dc=1000)

# Doncaster's non-parametric test of interaction
Don(SSL774, SSL775, tc=60*60, dc=1000)

# Read in ship data
# Pull up list of AIS files
files <- paste0("../Data_Raw/Ships/", list.files("../Data_Raw/Ships", pattern='.csv'))

# Create a list of lists of all csv file names grouped by month
temp <-  lapply(files, read.csv, header=TRUE, na.strings=c("","NA"),
                colClasses = c(rep("character", 2), "NULL", "character", "NULL", "NULL", 
                               "character", rep("NULL", 6), "character", "NULL", rep("character",8),
                               "NULL", "character", "NULL", "character", "NULL", rep("character", 2), rep("NULL", 109)
                ))

AIScsv <- do.call(rbind , temp)


test <- AIScsv %>% mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>% 
  filter(Latitude > 56, Latitude < 60, AIScsv$Longitude > -155, AIScsv$Longitude < -148)

