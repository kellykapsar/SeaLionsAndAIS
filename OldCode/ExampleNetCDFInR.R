# Example netcdf code from https://rpubs.com/boyerag/297592

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

nc_data <- nc_open('../Data_Dummy/gimms3g_ndvi_1982-2012.nc4')
# Save the print(nc) dump to a text file
{
  sink('../Data_Dummy/gimms3g_ndvi_1982-2012_metadata.txt')
  print(nc_data)
  sink()
}

# The following code reads the latitudes, longitudes, 
# and time of each NDVI observation and saves them in memory.

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

head(lon) # look at the first few entries in the longitude vector

# Read in the data from the NDVI variable and verify the dimensions of the array. 
ndvi.array <- ncvar_get(nc_data, "NDVI") # store the data in a 3-dimensional array
dim(ndvi.array) 

# Other pertinent information about the NDVI variable: 
# Lets’s see what fill value was used for missing data.
fillvalue <- ncatt_get(nc_data, "NDVI", "_FillValue")
fillvalue

nc_close(nc_data) 

# Replace fill value with NA
ndvi.array[ndvi.array == fillvalue$value] <- NA

# Time is the third dimension of the “ndvi.array”. 
# The first time slice represents the growing season of 1982.
ndvi.slice <- ndvi.array[, , 1] 

# Double check dimensions (should be 2d)
dim(ndvi.slice)

# Ok, everything checks out, so we can go ahead and save this data in a raster. 
# Note that we provide the coordinate reference system “CRS” in the standard well-known 
# text format. For this data set, it is the common WGS84 system.
r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# We will need to transpose and flip to orient the data correctly. 
# The best way to figure this out is through trial and error, 
# but remember that most netCDF files record spatial data from the bottom left corner.
r <- flip(r, direction='y')

# Plot the raster
plot(r)

# Save to GeoTIFF
# writeRaster(r, "GIMMS3g_1982.tif", "GTiff", overwrite=TRUE)

###############################
# Extract data at a study site
###############################

r_brick <- brick(ndvi.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# note that you may have to play around with the transpose (the t() function) and flip() 
# before the data are oriented correctly. In this example, the netcdf file recorded 
# latitude on the X and longitude on the Y, so both a transpose and a flip in the y 
# direction were required.
r_brick <- flip(t(r_brick), direction='y')

# Extract timeseries of data at the Toolik Lake study location from the raster brick 
# using the ‘extract()’ function.
toolik_lon <- -149.5975
toolik_lat <- 68.6275
toolik_series <- extract(r_brick, SpatialPoints(cbind(toolik_lon,toolik_lat)), method='simple')

# This timeseries is in a simple vector indexed only by the raster layer ID, so let’s 
# put it in an easier-to-use dataframe form and then plot the timeseries.
toolik_df <- data.frame(year= seq(from=1982, to=2012, by=1), NDVI=t(toolik_series))
ggplot(data=toolik_df, aes(x=year, y=NDVI, group=1)) +
  geom_line() + # make this a line plot
  ggtitle("Growing season NDVI at Toolik Lake Station") +     # Set title
  theme_bw() # use the black and white theme

###########################################
# Difference in NDVI between 2 time periods
###########################################

# The ‘ndvi.slice’ array has the data from 1982. 
# Let’s get the data from 2012, the 31st time slice.
ndvi.slice.2012 <- ndvi.array[, , 31] 

# Now take the difference between them 
ndvi.diff <- ndvi.slice.2012 - ndvi.slice

# Save the difference map as a raster
r_diff <- raster(t(ndvi.diff), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Re-orient the data for geotiff.
r_diff <- flip(r_diff, direction='y')

# Plot
plot(r_diff)
