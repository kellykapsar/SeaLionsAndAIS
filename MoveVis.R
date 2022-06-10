#
#
#
# Movevis tutorial 
# Copied from: http://movevis.org/
# 
# 
# 


library(moveVis)
library(move)
library(sf)
library(dplyr)
library(RColorBrewer)
library(RStoolbox)

# Projection information for WGS84/UTM Zone 5N (EPSG:32605)
prj_latlon <- "+proj=longlat +datum=WGS84 +no_defs"
prj_utm <- "+proj=utm +zone=5 +datum=WGS84 +units=m +no_defs"

ssl <- readRDS("../Data_Processed/Telemetry/watersealis.rds") %>% 
          as.data.frame() %>% 
          dplyr::select(deploy_id, date, lat, lon, speed_kmhr) %>% 
          mutate(deploy_id=as.factor(deploy_id))
ssl <- ssl[-which(is.nan(ssl$speed_kmhr)),]


move_kod <- df2move(ssl[which(grepl("KOD", ssl$deploy_id) == TRUE),], 
                    proj= prj_latlon, x="lon", y="lat", time = "date", track_id="deploy_id")
move_pws <- df2move(ssl[which(grepl("PWS", ssl$deploy_id) == TRUE),], 
                    proj= prj_latlon, x="lon", y="lat", time = "date", track_id="deploy_id") 

sealilabels<- data.frame(names = c("SSL2018774PWS", "SSL2018775PWS", "SSL2018776PWS", "SSL2018777PWS", 
                                   "SSL2019781KOD", "SSL2019782KOD", "SSL2019783KOD", "SSL2019784KOD",
                                   "SSL2019785KOD", "SSL2019786KOD", "SSL2019788KOD"), 
                         labels=c("774PWS", "775PWS", "776PWS", "777PWS", 
                                  "781KOD", "782KOD", "783KOD", "784KOD", 
                                  "785KOD", "786KOD", "788KOD"), 
                         colors= c("#0070ff", "#002673", "#b2df8a", "#33a02c",
                                   "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                                   "#cab2d6", "#8967ae", "#d5d000"))

# data("move_data", package = "moveVis") # move class object
# if your tracks are present as data.frames, see df2move() for conversion

# lags <- timeLag(move_kod, unit = "mins")
# lapply(lags, mean)

# lags <- timeLag(move_pws, unit = "mins")
# lapply(lags, mean)

# # Create custom basemap based on fishing data 
# fish <- readRDS("../Data_Processed/AIS_Fishing.rds") %>% raster::projectRaster(crs=prj_latlon)
# fish[fish < 0 ] <- 0
# 
# # Calculated weeks manually because it's really difficult to convert "2018-W44" to a date format
# fish_times <- seq(from = ISOdate(2018, 11, 01), to = ISOdate(2020, 07, 31), by="week")
# 
##################
### KODIAK GIF ###
##################

# # align move_data to a uniform time scale
# m_kod <- align_move(move_kod, res = 12, unit = "hours")
# 
# # Get colors for individuals
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# # create spatial frames with a OpenStreetMap watercolour map
# frames_kod <- frames_spatial(m_kod, path_colours = sealilabels$colors[5:11],
# #                            r_list = fish_kod[[1]], r_times=fish_kod[[2]], r_type="gradient")
# #                            map_service = "mapbox", map_type = "satellite",
# #                            map_token = "pk.eyJ1IjoibXlkb2dzd2Fsa21lIiwiYSI6ImNsMWR5MXA1ejAyY2czaXFrbHloc3F3ZmsifQ.Ja8O59-J527AKSPkEQAtVg") %>%
#                          map_service = "osm", map_type = "watercolor", alpha = 0.5) %>%
#   add_labels(x = "", y = "", title = "") %>% # Steller Sea Lion Movements, Nov. 2018 - Jul. 2019
#   add_scalebar(colour="black", distance=80, position = "bottomright") %>%
#   add_timestamps(m_kod, type = "label") %>%
#   add_progress()
# 
# # Add raster layer for marine vessel trafficon top
# # ts_kod <- sort(unique(timestamps(m_kod)))
# # r_pos_kod <- sapply(ts_kod, function(x) which.min(abs(difftime(x, fish_times))))
# # fish_kod <- lapply(r_pos_kod, function(i) fish[[i]])
# #
# #
# # frames_kod_fish <- add_gg(frames_kod, gg = expr(RStoolbox::ggR(data, alpha = 0.5, ggLayer = T)),
# #                           stretch="log", data = fish_kod)
# 
# # animate frames
# # animate_frames(frames_kod_fish, fps = 10, out_file = paste0("../Figures/moveVis_KOD_fish.gif"))
# animate_frames(frames_kod, fps = 10, out_file = paste0("../Figures/moveVis_KOD_10fps_osm.gif"))
# browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ&ab_channel=RHINO")
###############
### PWS GIF ###
###############

# align move_data to a uniform time scale
m_pws <- align_move(move_pws, res = 12, unit = "hours")

# Get colors for individuals
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# create spatial frames with a OpenStreetMap watercolour map
frames_pws <- frames_spatial(m_pws, path_colours = sealilabels$colors[1:4],
  #                            map_service = "mapbox", map_type = "satellite",
  #                            map_token = "pk.eyJ1IjoibXlkb2dzd2Fsa21lIiwiYSI6ImNsMWR5MXA1ejAyY2czaXFrbHloc3F3ZmsifQ.Ja8O59-J527AKSPkEQAtVg") %>%
                               map_service = "osm", map_type = "watercolor", alpha = 0.5) %>%
  add_labels(x = "", y = "", title = "") %>% # Steller Sea Lion Movements, Nov. 2018 - Jul. 2019
  add_scalebar(colour="black", distance=80) %>%
  add_timestamps(m_pws, type = "label") %>%
  add_progress()

# Add raster layer for marine vessel trafficon top
# Code from: https://gist.github.com/16EAGLE/4bfb0ca589204c53041244aa705b456b
# ts_pws <- sort(unique(timestamps(m_pws)))
# r_pos_pws <- sapply(ts_pws, function(x) which.min(abs(difftime(x, fish_times))))
# fish_pws <- lapply(r_pos_pws, function(i) fish[[i]])
# 
# frames_pws_fish <- add_gg(frames_pws, gg = expr(RStoolbox::ggR(data, alpha = 0.5, ggLayer = T)), 
#                      stretch="log", data = fish_pws)

# animate frames
# animate_frames(frames_pws_fish, fps = 10, out_file = paste0("../Figures/moveVis_PWS_fish.gif"))
animate_frames(frames_pws, fps = 10, out_file = paste0("../Figures/moveVis_PWS_10fps_osm.gif"))
browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ&ab_channel=RHINO")