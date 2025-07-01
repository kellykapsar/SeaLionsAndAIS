

library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(fasterize)
library(foreach)
library(doParallel)

# weekly homerange polygons
raw <- read_rds("../Data_Processed/ssl_ak_30min.rds") %>% 
  mutate(day = lubridate::date(t_))

# Vessel tracklines
daily_ships <- list.files("/mnt/research/CSIS/AIS/Data_Processed_V4/Points", full.names = T)

t <- read.csv(daily_ships[1]) %>%
  mutate(Time = as.POSIXct(Time, tz = "GMT")) %>% 
  arrange(Time)

raw <- raw %>% nest_by(day)

r <- raw$data[[2]] %>%
  st_as_sf(coords = c("x_", "y_"), crs = 32605) %>% 
  st_buffer(10000)

# Step 2: Create a bounding box from the points
bbox <- r %>% 
  st_bbox() %>%
  st_as_sfc()


# Step 4: Transform to EPSG:3338
bbox_buffered_3338 <- st_transform(bbox, crs = 3338)

# You can now visualize or use this bounding box as needed
print(bbox_buffered_3338)
st_coordinates(bbox_buffered_3338)


u <- as.numeric(difftime(t$Time, lag(t$Time), units = "mins")) %>% max()



