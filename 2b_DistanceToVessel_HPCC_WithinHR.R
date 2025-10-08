

library(tidyr)
library(readr)
library(dplyr)
library(sf)
library(raster)
library(stars)
library(fasterize)
library(foreach)
library(doParallel)
library(lubridate)
library(terra)
library(gdistance)

# Load Functions 
#' Load a CSV File by Date from a List of Filenames
#'
#' This function searches a list of filenames for a file containing a given date
#' (formatted as "YYYY_MM_DD") and loads it as a data frame.
#'
#' @param date_input A `Date` object or character string (e.g., "2017-09-26")
#'   representing the date to match within the filenames.
#' @param file_list A character vector of file paths or filenames that follow
#'   a naming convention including the date in "YYYY_MM_DD" format.
#'
#' @return A data frame read from the matched CSV file.
#'
#' @examples
#' \dontrun{
#' file_list <- list.files(path = "data/", pattern = "Clean_Points_.*\\.csv$", full.names = TRUE)
#' df <- load_csv_by_date("2017-09-26", file_list)
#' head(df)
#' }
#'
#' @export
load_csv_by_date <- function(date_input, file_list) {
  # Ensure date_input is in Date format
  date_input <- as.Date(date_input)
  
  # Format date to match file pattern: "YYYY_MM_DD"
  date_str <- format(date_input, "%Y_%m_%d")
  
  # Find the file that contains the date string
  matched_file <- grep(date_str, file_list, value = TRUE)
  
  if (length(matched_file) == 0) {
    stop(paste("No file found for date:", date_str))
  } else if (length(matched_file) > 1) {
    warning(paste("Multiple files found for date:", date_str, "- using the first match"))
  }
  
  # Load the matched file
  df <- read.csv(matched_file[1])
  
  return(df)
}

#' Filter points within a bounding box
#'
#' This function removes rows from a dataframe that have coordinates
#' outside a specified bounding box.
#'
#' @param df A data frame containing `x` and `y` coordinate columns in EPSG:3338.
#' @param bbox An `sf` or `sfc` bounding box object defining the spatial filter area.
#'
#' @return A filtered data frame containing only rows within the bounding box.
#'
#' @examples
#' \dontrun{
#' filtered_df <- filter_points_in_bbox(df, bbox_buffered_3338)
#' }
#'
#' @export
filter_points_in_bbox <- function(df, bbox) {
  # Ensure bounding box is in the form of coordinates
  bbox_coords <- sf::st_bbox(bbox)
  
  # Filter by coordinate limits
  df_filtered <- df %>%
    dplyr::filter(
      x >= bbox_coords["xmin"],
      x <= bbox_coords["xmax"],
      y >= bbox_coords["ymin"],
      y <= bbox_coords["ymax"]
    )
  
  return(df_filtered)
}




################################################################################

# Regularized sea lion locations 
raw <- read_rds("../Data_Processed/ssl_ak_30min.rds") %>% 
  mutate(day = lubridate::date(t_), 
         t_ = with_tz(t_, tzone = "UTC")) %>% 
  nest_by(day)

r <- raw$data[[2]] %>%
  st_as_sf(coords = c("x_", "y_"), crs = 32605, remove=F) 

r_buff <- r %>% 
  st_buffer(100000) # 100 km buffer around all sea lion points within a day %>% 

# Step 2: Create a bounding box from the points
bbox <- r_buff %>% 
  st_bbox() %>%
  st_as_sfc()


# Step 4: Transform to EPSG:3338
bbox_buffered_3338 <- st_transform(bbox, crs = 3338)

# Vessel tracklines
daily_ships <- list.files("/mnt/research/CSIS/AIS/Data_Processed_V4/Points", full.names = T)

df <- load_csv_by_date(raw$day[2], daily_ships) 
df_in <- df %>% filter_points_in_bbox(., bbox_buffered_3338) %>% 
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>% 
  st_as_sf(coords = c("x", "y"), crs = 3338) %>% 
  st_transform(32605) %>% 
  mutate(x = st_coordinates(.)[,1], 
         y = st_coordinates(.)[,2])

plot(df_in$geometry)
plot(r$geometry, add=T, col="red")

################################################################################
# Testing out linearly interpolated ship positions 

sl <- r %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% 
  rename(x = x_, y = y_, t_sl = t_) %>% 
  dplyr::select(deploy_id, t_sl, x, y)

ais <- df_in %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% 
  rename(t_ais = Time, mmsi = scramblemmsi) %>% 
  dplyr::select(mmsi, t_ais, x,y)


# Get unique ship IDs
ships <- unique(ais$mmsi)

# Cross join: one row per deploy_id x t_sl x mmsi
sl_exp <- sl %>%
  tidyr::expand_grid(mmsi = ships)

# Join all AIS points to the expanded SL fixes
joined <- sl_exp %>%
  left_join(ais, by = "mmsi") %>%
  filter(!is.na(t_ais))



# Get previous AIS fix (t_ais <= t_sl)
prev <- joined %>%
  group_by(deploy_id, t_sl, mmsi) %>%
  filter(t_ais <= t_sl) %>%
  filter(t_ais == max(t_ais)) %>%
  rename(t0 = t_ais, x0 = x.y, y0 = y.y) %>%
  ungroup()

# Next AIS fix (t_ais >= t_sl)
next_ <- joined %>%
  group_by(deploy_id, t_sl, mmsi) %>%
  filter(t_ais >= t_sl) %>%
  filter(t_ais == min(t_ais)) %>%
  rename(t1 = t_ais, x1 = x.y, y1 = y.y) %>%
  ungroup()


interp_data <- prev %>%
  inner_join(next_, by = c("deploy_id", "t_sl", "mmsi", "x.x", "y.x")) %>%
  filter(difftime(t1, t0, units = "secs") <= 60 * 60) %>%
  mutate(
    w = as.numeric(difftime(t_sl, t0, units = "secs")) /
      as.numeric(difftime(t1, t0, units = "secs"))) %>% 
  mutate(w = ifelse(is.nan(w), 0, w)) %>% 
  mutate(x_ship = x0 + w * (x1 - x0),
    y_ship = y0 + w * (y1 - y0),
    euc_dist_m = sqrt((x.x - x_ship)^2 + (y.x - y_ship)^2)
  )

# Land raster
land <- raster("../Data_Processed/Bathymetry.tif")
land[values(land) > 0] <- NA
land[values(land) <= 0] <- 1




interp_data$water_dist_m <- NA

for (i in 1:length(sl$deploy_id)) {
  print(i)
  
  this_sl <- sl[i,]
  r_cost <- land
  
  # Create single-point origin raster
  cell_i <- cellFromXY(r_cost, xy= as.double(sl[i, c("x", "y")]))
  
  r_cost[cell_i] <- 2
  
  # Compute grid distance
  d <- gridDistance(r_cost, origin = 2, omit=NA)
  
  # Subset interp_data to all ships for this sea lion fix
  sl_idx <- which(
    interp_data$deploy_id == this_sl$deploy_id &
      interp_data$t_sl == this_sl$t_sl
  )
  
  # Extract to the ship point
  ship_xy <- interp_data[sl_idx, c("x_ship", "y_ship")]
  interp_data$water_dist_m[sl_idx] <- extract(d, matrix(unlist(ship_xy), ncol = 2))
  
  rm(r_cost)
}

write.csv(interp_data, "../temp_one_day_distance_test.csv")

st_transform(st_as_sf(bbox_buffered_3338), 32605) -> b

land_cropped <- terra::crop(land, b)
bbox_buffered_3338

test <- interp_data[!is.na(interp_data$water_dist_m),]
test <- test[!is.na(test$euc_dist_m),]
cor(test$euc_dist_m, test$water_dist_m)
plot(test$euc_dist_m, test$water_dist_m)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2) 

# Number of ships within 100 km of SL that day
# Number of ships unable to calculate the distance (e.g., points intersected with land)
# Closest ship distance that day 