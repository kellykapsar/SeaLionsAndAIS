# Title: SSL Data Processing V2
# Author: Kelly Kapsar, Sydney Waloven
# Date: 2024-05-16
# Description: R script version of 1a_SSLDataProcessing.Rmd for readability. Original script by Kelly Kapsar. Edits made by Sydney Waloven to regularize sea lion locations to a constant time step.


# Setup -------------------------------------------------------------------

library(sf)
library(raster)
library(ggplot2)
library(scales)
library(ggmap)
library(leaflet)
library(RColorBrewer)
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

basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>%
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


# Plot labels and color palette -------------------------------------------

sealilabels <- data.frame(names = c("SSL2018774PWS", 
                                    "SSL2018775PWS", 
                                    "SSL2018776PWS",
                                    "SSL2018777PWS", 
                                    "SSL2019781KOD", 
                                    "SSL2019782KOD", 
                                    "SSL2019783KOD", 
                                    "SSL2019784KOD",
                                    "SSL2019785KOD", 
                                    "SSL2019786KOD",
                                    "SSL2019788KOD"), 
                          labels = c("774PWS",
                                     "775PWS", 
                                     "776PWS", 
                                     "777PWS",
                                     "781KOD", 
                                     "782KOD", 
                                     "783KOD", 
                                     "784KOD",
                                     "785KOD", 
                                     "786KOD",
                                     "788KOD"), 
                          colors = c("#0070ff", "#002673", "#b2df8a", "#33a02c",
                                     "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                                     "#cab2d6", "#8967ae", "#d5d000"))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))


# Sea lion location geodatabase -------------------------------------------

# Import sea lion location geodatabase (downloaded from Google Drive folder)
seali <- list.files("../Data_Raw/raw data files - complete - all tags-20210803T143355Z-001")
# seali$name # list of layers 

# Determine all layers of interest within gdb
# lyrs <- grep("_loc", seali$name)
lyrs <- grep("S-Locations", seali)
lyrs2 <- grep("D-Locations", seali)

lyrs <- append(lyrs, lyrs2)

# Create initial sf object with data from one sea lion
lyrs <- lapply(lyrs, function(x) {
  st_read(paste0("../Data_Raw/raw data files - complete - all tags-20210803T143355Z-001/", seali[x]))})

# Make each table into an sf object
# lyrs <- lapply(lyrs, function(x) {st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326)})

# Separate first table layer
sealis <- lyrs[[1]]

# Remove that layer from the layers of interest list
lyrs <- lyrs[2:length(lyrs)]

# Append all other layers of interest onto the main location data set 
for(i in 1:length(lyrs)) {
  temp <- lyrs[[i]]
  sealis <- rbind(sealis, temp)
}

# Change all column names to lowercase
colnames(sealis) <- tolower(colnames(sealis))

# Save clean seali data 
sealis <- rename(sealis, 
                 deploy_id = deployid,
                 error_radius = error.radius, 
                 error_major = error.semi.major.axis, 
                 error_minor = error.semi.minor.axis,
                 error_ellipse = error.ellipse.orientation)

# Fix time field 
# (Have to do it separately for gps and argos)
sealis$date_old <- sealis$date
gps <- sealis[sealis$type == "FastGPS", ]
argos <- sealis[sealis$type == "Argos", ]

argos$date <- as.POSIXct(argos$date, 
                         format = c("%Y/%m/%d %H:%M:%S"),
                         tz = "GMT")
# All but 353 gps points are in the format of days since 12/30/1899
# Need to convert those to dates and then also fix the other 300ish points 
# Internet said it should be days since 1/1/1900, but that didn't work. No idea why. 
# But this matches up with the Microsoft Access database
gps$date <- as.POSIXct("1899-12-30 00:00:00", tz = "GMT") + 
  as.difftime(as.numeric(gps$date), units = "days")

gps$date[is.na(gps$date)] <- as.POSIXct(gps$date_old[is.na(gps$date)], 
                                        format = c("%Y/%m/%d %H:%M:%S"),
                                        tz = "GMT")

sealis_old <- sealis 

# Rejoin Argos and GPS data 
sealis <- rbind(gps, argos)

# Round date to minute scale 
sealis$date <- round(sealis$date, "mins")

# Create various date reference categories for future modeling 
sealis$fortnight <- ceiling(lubridate::week(sealis$date) / 2)
sealis$weekofyear <-  format(sealis$date, "%G-W%V")
sealis$month <- lubridate::month(sealis$date)
sealis$year <- lubridate::year(sealis$date)
sealis$dayofyear <- lubridate::date(sealis$date)

# Remove low quality points 
sealis <- sealis[-which(sealis$quality %in% c("A", "B", "0", "Z")), ]

# 40 duplicated rows (excluding geometry column)
test <- duplicated(data.frame(sealis))
# Remove duplicated rows 
sealis <- sealis[which(duplicated(data.frame(sealis)) == FALSE), ]

# Remove points from same time and SSL (keep smaller error radius)
sealis$ptID <- 1:length(sealis$deploy_id)

test <- sealis %>% 
  filter(type == "Argos") %>% 
  group_by(deploy_id, date) %>% 
  slice(which.min(error_radius))

test2 <- filter(sealis, type != "Argos")

sealis <- rbind(test, test2)

# Put back in temporal order by sea lion
sealis <- sealis[order(sealis$deploy_id, sealis$date), ]

# Convert to spatial object
sealis <- st_as_sf(sealis,
                   coords = c("longitude", "latitude"), 
                   crs = 4326)

# Keep latitude and longitude columns
sealis$lat <- st_coordinates(sealis)[,"Y"]
sealis$lon <- st_coordinates(sealis)[,"X"]

# Transform geometry to correct projection 
st_geometry(sealis) <- st_transform(st_geometry(sealis), prj)

# Projected coordinate columns
sealis$northing <- st_coordinates(sealis)[,"Y"]
sealis$easting <- st_coordinates(sealis)[,"X"]

# Remove blank columns 
sealis <- subset(sealis, 
                 select = -c(offset, 
                             offset.orientation, 
                             gpe.msd,
                             gpe.u,
                             count))

# Implement speed filter of 3 m/s on all locations using McConnell algorithm
## time between sampling intervals has not been standardized yet so we should not yet decide about removing points based on speed ## ~SW
sealis <- sealis %>% 
  arrange(deploy_id, date) %>% 
  group_by(deploy_id) %>% 
  mutate(spdfilt = 
           argosfilter::vmask(lat = lat, 
                              lon = lon,
                              dtime = date,
                              vmax = 3)) %>% 
  ungroup()
# vmask will flag points as "valid" or "removed" based on if the speed threshold is exceeded

# Remove speed filtered points that exceed the 3 m/s threshold
sealis_dirty <- sealis 
sealis <- sealis %>% 
  dplyr::filter(spdfilt != "removed") %>% # removes invalid points
  dplyr::arrange(date)

# Plot change between speed filtered and original, separated by ssl (Steller sea lions data)
sealis_dirty %>% st_drop_geometry() %>% 
  dplyr::arrange(date) %>% 
  ggplot() + 
  geom_path(aes(x = lon, y = lat)) + 
  geom_path(data = sealis, aes(x = lon, y = lat), col = "red") + 
  facet_wrap(~ deploy_id, scales = "free")

# Plot clean data by ssl - plot change between speed filtered and original, separated by ssl
sealis %>% st_drop_geometry() %>% 
  dplyr::arrange(date) %>% 
  ggplot() + 
  geom_path(aes(x = lon, y = lat)) + 
  facet_wrap(~ deploy_id, scales = "free")

# Map of study area 
ggplot() +
  geom_sf(data = basemap.crop, 
          fill = "gray", 
          color = "black",
          lwd = 0.5, 
          alpha = 0.8) +
  geom_sf(data = study, fill = NA, color = "brown") +
  # geom_path is much faster than geom_sf
  geom_path(data = sealis, 
            aes(x = easting, y = northing, color = deploy_id),
            key_glyph = "rect") +
  scale_color_manual(values = sealilabels$color, 
                     name = "Individual", 
                     labels = sealilabels$labels) +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 15)) 

ggsave("../Figures/PathMap.png", 
       width = 10, height = 8, units = "in")

# Tag duration timeline 
timeline <- sealis %>% 
  st_drop_geometry() %>% 
  group_by(deploy_id) %>% 
  summarize(starttag = as.Date(min(date)), 
            endtag = as.Date(max(date)))

ggplot(timeline, 
       aes(x = starttag, y = deploy_id, color = deploy_id)) +
  geom_linerange(aes(xmin = starttag, xmax = endtag), 
                 lwd = 2, 
                 show.legend = FALSE) + 
  scale_color_manual(values = sealilabels$color) +
  scale_x_date(breaks = date_breaks(width = "1 month"), 
               date_labels = "%b %Y") +
  scale_y_discrete(labels = sealilabels$labels) +
  ylab("Individual") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15))

ggsave("../Figures/TelemetryTimeline.png", 
       width = 10, height = 6, units = "in")


# Available locations - radius method -------------------------------------

# Calculating a (average hrly movement rate in km/hr) for each SSL and b (sd of movement rate)
euclidean_speed <- function(lat2, lat1, long2, long1, time2, time1) {
  latdiff <- lat2 - lat1
  longdiff <- long2 - long1
  distance <- sqrt(latdiff^2 + longdiff^2)/1000
  timediff <- as.numeric(difftime(time2, time1, units = c("hours")))
  return(distance / timediff)
}

# Recalculate Euclidean speed
sealis <- sealis %>% 
  group_by(deploy_id) %>%
  arrange(deploy_id, date) %>% 
  mutate(timediff = as.numeric(difftime(date, lag(date), units = c("mins"))),
         speed_kmhr = euclidean_speed(northing, lag(northing), 
                                      easting, lag(easting), 
                                      date, lag(date)))


# Need to clean out inf and NA speed as well as those above a certain threshold
# How to determine threshold?
sealispeed <- sealis %>% 
  st_drop_geometry() %>% 
  group_by(deploy_id) %>% 
  summarize(speed_avg = mean(speed_kmhr, na.rm = T), 
            speed_sd = sd(speed_kmhr, na.rm = T), 
            timediff = mean(timediff, na.rm = T)/60) # convert to hourly

# radius = c(a + 2b)
sealispeed$radius <- sealispeed$timediff*(sealispeed$speed_avg + 2*sealispeed$speed_sd) 

sealis <- left_join(sealis, sealispeed, by = "deploy_id")


# Average number of non-land points per sea lion per time period ----------
# Run script 1b study area, bathymetry, and landmask sections

# Remove locations at sea lion haul-out sites using landmask

# Read in landmask 
landmask <- raster("../Data_Processed/Landmask_GEBCO.tif")

# Extract landmask values at each sea lion location
sealis$land <- raster::extract(landmask, sealis)
# Calculate the proportion of points that are on land
sum(sealis$land, na.rm = T) / length(sealis$land)*100 # 20% of locations on land

# Subset sea lion locations to only those in the water (where landmask = NA)
watersealis <- sealis[is.na(sealis$land), ]

# Now calculate the mean number of seal locations in the water monthly, biweekly, weekly, and daily
ptcts_month <- watersealis %>% 
  st_drop_geometry() %>%
  group_by(deploy_id, year, month) %>% 
  summarize(n = n())
ptcts_month <- ptcts_month %>%
  group_by(deploy_id) %>% 
  summarize(meanmonthlypts = mean(n)) # calculate the mean number of points per month for each id

ptcts_biweek <- watersealis %>% 
  st_drop_geometry() %>% 
  group_by(deploy_id, year, fortnight) %>% 
  summarize(n = n())
ptcts_biweek <- ptcts_biweek %>% 
  group_by(deploy_id) %>% 
  summarize(meanbiweekpts = mean(n)) # calculate the mean number of points biweekly for each id

ptcts_week <- watersealis %>% 
  st_drop_geometry() %>% 
  group_by(deploy_id, year, weekofyear) %>%
  summarize(n=n())
ptcts_week <- ptcts_week %>% 
  group_by(deploy_id) %>% 
  summarize(meanweekpts = mean(n)) # calculate the mean number of points weekly for each id

ptcts_day <- watersealis %>%
  st_drop_geometry() %>%
  group_by(deploy_id, dayofyear) %>% 
  summarize(n = n())
ptcts_day <- ptcts_day %>% 
  group_by(deploy_id) %>% 
  summarize(meandaypts = mean(n)) # calculate the mean number of points daily for each id

ptcts <- left_join(ptcts_month, ptcts_biweek, by = "deploy_id")
ptcts <- left_join(ptcts, ptcts_week, by = "deploy_id")
ptcts <- left_join(ptcts, ptcts_day, by = "deploy_id")

# write.csv(ptcts, "../Data_Processed/SSL_PtCts.csv")

# Save clean data ---------------------------------------------------------

# Add in unique id for each used point
watersealis$choice_id <- 1:length(watersealis$deploy_id)

# Save clean data output 
# saveRDS(watersealis, "../Data_Processed/Telemetry/watersealis.rds")

# amt: creating tracks for all sea lions ----------------------------------


# Convert used locations to a track object using amt package
library(amt)

# Read in watersealis data
watersealis <- read_rds("../Data_Processed/Telemetry/watersealis.rds")

# Referencing Animal Space Use and Behavior SMSC with Dr. Joe Kowalski -----
seali_data <- watersealis %>% 
  # Convert sf object into a tibble
  as_tibble() %>% 
  # Remove the geometry column
  select(-geometry)
# This returns our sf object to its original data frame, the projected coordinates are in the "northing" (Y) and "easting" (X) columns.

# Here I'm adding nsd (net squared displacement) in the same code chain that I make the track, and I'm also specifying that we have an id column, so each animal can be kept separate.
seali_tracks <- seali_data %>%
  # Specify x and y coordinate columns and date column
  make_track(.x = easting,
             .y = northing,
             .t = date,
             id = deploy_id, # differentiates each animal
             crs = 32605,
             all_cols = TRUE) %>%  # ensures all columns are brought forward
  add_nsd() %>% # net squared displacement
  # Add speed information to track file
  mutate(speed = speed(.))
  
  
# Nest data except for deploy_id
seali_tracks2 <- seali_tracks %>% 
  nest(-deploy_id)

# Get sampling rate for each animal
sr_all <- seali_tracks2 %>% 
  # Create a new nested column to save the sampling rate summary info
  mutate(sr = map(data, summarize_sampling_rate)) %>% 
  # Select the id and sr col, remove all other info in the data column
  select(deploy_id, sr) %>% 
  # Unnest the data in the sr column to make it easily viewable
  unnest(cols = c(sr))

# median sampling rate for each sea lion
# 29 29 28 25 20 30 24 30 18 27 37

# Convert to df to visualize sr for each animal
sr_df <- data.frame(deploy_id = sr_all$deploy_id,
                    median = sr_all$median,
                    mean = sr_all$mean,
                    n = sr_all$n,
                    unit = sr_all$unit)
sr_df

# Visualize the number of locations per day for all animals
seali_tracks %>% 
  group_by(date = as_date(t_)) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = date,
             y = count)) +
  geom_col() +
  ylim(c(0, 20))


# amt: resample tracks ----------------------------------------------------


# Resample track information to regularized sampling interval

ssl_steps <- seali_tracks2 %>% 
  # Apply a function to x (each item in nested 'data' column)
  mutate(steps = map(data, function(x)
    # Tags were programmed to transmit every 15 min. Median sr was about 30 min so we'll resample to 30.
    x %>% track_resample(rate = minutes(30),
                         tolerance = minutes(3)) %>% 
      # Only keep bursts (strings of locations) with 3 consecutive locations
      filter_min_n_burst(min_n = 3) %>% 
      # Keep attribute associated with the end of the step
      steps_by_burst(keep_cols = "end")))

str(ssl_steps, width = 80, strict_width = "no", nchar.max = 80, give.attr = FALSE)

# Take a look at what was lost by choosing this sampling rate. This is showing the number of step lengths that were ultimately calculated for each day of the study period. If the sea lions had periods of sampling at different intervals, then this would show what time periods were lost from the resampling decision.

ssl_steps %>% 
  select(deploy_id, steps) %>% 
  unnest(cols = steps) %>% 
  group_by(date = as_date(t2_),
           deploy_id) %>% 
  summarize(count = n()) %>% # number of step lengths per day of study period
  ggplot(aes(x = date,
             y = count)) +
  geom_col() +
  facet_wrap(~ deploy_id) +
  labs(y = "steps per day")

# Select deploy_id and steps, unnest the data frame and create a plot of step-length distributions

ssl_steps %>% 
  # Select deploy_id and steps
  select(deploy_id, steps) %>% 
  # Unnest new df
  unnest(cols = steps) %>% 
  # Generate plot of step-length distributions
  ggplot(aes(x = sl_,
             fill = deploy_id)) +
  geom_density(alpha = 0.4) +
  labs(x = "Step length [m]",
       y = "Density") +
  xlim(c(0, 5000)) +
  theme_light()

# Save data
saveRDS(seali_tracks, "../Data_Processed/seali_tracks.rds")
saveRDS(ssl_steps, "../Data_Processed/ssl_steps_resampled.rds")

# Plot individual SSL locs ------------------------------------------------

# # Using ggplot without a background

# Use separate axes for each individual (add scales = "free" to facet_wrap
# + fig.height = 12, fig.width = 12
ggplot(sealis,
       aes(x = lon, y = lat)) +
  geom_point() +
  facet_wrap(~ deploy_id, scales = "free")

# Test individual seals
ssl781 <- sealis %>%
  filter(deploy_id == "SSL2019781KOD")

# Leaflet map - works but there are no points on the map
leaflet(ssl781) %>% 
  addTiles() %>% 
  addCircles(~lon, ~lat)

# ggmap - doesn't work
map <- get_map(location = c(lon = mean(ssl781$lon),
                            lat = mean(ssl781$lat)),
               zoom = 7,
               maptype = "hybrid",
               source = "google")

ggmap(map) +
  geom_point(data = ssl781,
             aes(x = lon, y = lat), size = 2.5)
