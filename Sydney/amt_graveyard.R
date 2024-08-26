# amt supplemental material from Signer et al. 2019 -----
# #  Optional section for easier visualization # # 

# Create tracks using amt code from Signer et al. 2019. This provides a more readable and easier visualization of the sampling summaries.
seali_data_amt <- seali_data %>% 
  nest(-deploy_id)

seali_tracks_sr <- seali_data_amt %>% 
  mutate(trk = map(data, function(d) {
    make_track(d,
               .x = easting,
               .y = northing,
               .t = date,
               crs = 32605,
               all_cols = TRUE)
  }))

# Visualize sampling rate summaries for each animal
seali_tracks_sr %>% 
  mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  select(deploy_id, sr) %>% 
  unnest(cols = c(sr)) 

# end amt supplemental section -----