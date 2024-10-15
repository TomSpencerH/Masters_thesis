library(tidyverse)
library(metR) 
library(zoo)
library(rgdal)
library(sf)
library(sp)
library(lwgeom)
library(FNN)
library(viridis)
library(rnaturalearthdata)
library(rnaturalearth)
library(rgeos)
library(maps)

bind <- readRDS("c:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/Processed_data/Wind/Francois_wind_14h00.Rds")

buffer_func <- function(buffer){
  
  # Download the coastline data for the entire world
  world_coastline <- ne_download(scale = "large", type = "coastline", category = "physical", returnclass = "sf")
  
  # Define the bounding box for Africa
  africa_bbox <- st_bbox(c(xmin = -20, ymin = -35, xmax = 52, ymax = 38), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box for Africa
  africa_coastline <- st_crop(world_coastline, africa_bbox)
  
  # The BCLME region:
  bclme_bbox <- st_bbox(c(xmin = 7.91755, ymin = -36.61979, xmax = 19.788742, ymax = -5.811113), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box
  bclme_coastline <- st_crop(africa_coastline, bclme_bbox)
  
  # Create a buffer of 50 nautical miles (1 nautical mile = 1852 meters)
  buffer_50nm <- st_make_valid(st_union(st_buffer(bclme_coastline, dist = buffer * 1852)))
  
  # Ensure coastline extends beyond buffer: Extend the coastline line if necessary
  # (I added/subtracted 3 degrees to the bbox to get a longer coastline)
  # Create a bbox slightly larger than the BCLME extent to get a longer coastline
  # that will insersect the buffer and thus can be used for splitting it lengthwise
  extended_bbox <- st_bbox(c(xmin = 7.91755 - 3,
                             ymin = -36.61979 - 3,
                             xmax = 19.788742 + 3,
                             ymax = -5.811113 + 3),
                           crs = st_crs(world_coastline))
  
  # Crop the coastline data to the outer bbox; this is the coastline that
  # will be used for splitting the buffer
  extended_coastline <- st_crop(africa_coastline, extended_bbox)
  
  # Convert the outer coastline to a single LINESTRING object
  extended_coastline_line <- st_union(st_cast(extended_coastline, "LINESTRING"))
  extended_coastline_line <- st_make_valid(extended_coastline_line)
  #extended_coastline_line_extended <- st_segmentize(extended_coastline_line, dfMaxLength = 1000000) # Adjust the max length as needed
  
  # Split the buffer using the extended coastline line
  split_buffers <- st_split(buffer_50nm, extended_coastline_line)
  
  # Extract the resulting polygons
  split_buffers_sf <- st_collection_extract(split_buffers)
  
  # Separate into two objects
  offshore_buffer <- split_buffers_sf[1, drop = FALSE]
  inland_buffer <- split_buffers_sf[2, drop = FALSE]
  
  # Ensure buffers are in the original CRS
  offshore_buffer <- st_transform(offshore_buffer, crs = st_crs(world_coastline))
  inland_buffer <- st_transform(inland_buffer, crs = st_crs(world_coastline))
  
  # Create simulated gridded data to test the buffer splitting
  # Define the bounding box for the region
  xmin <- 7.91755
  ymin <- -36.61979
  xmax <- 19.788742
  ymax <- -5.811113
  
  # Generate a grid of points within the bounding box
  lon <- seq(xmin, xmax, length.out = 200)
  lat <- seq(ymin, ymax, length.out = 400)
  grid <- expand.grid(lon = lon, lat = lat)
  
  # Add a simulated temperature column
  set.seed(42)  # For reproducibility
  grid$temperature <- runif(nrow(grid), min = 10, max = 30)
  
  # Convert the data frame to an sf object
  grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = st_crs(world_coastline))
  
  # Extract points within the inland buffer
  points_in_inland_buffer <- st_intersects(grid_sf, inland_buffer, sparse = FALSE)
  inland_data <- grid[apply(points_in_inland_buffer, 1, any), ]
  
  # Extract points within the offshore buffer
  points_in_offshore_buffer <- st_intersects(grid_sf, offshore_buffer, sparse = FALSE)
  offshore_data <- grid[apply(points_in_offshore_buffer, 1, any), ]
  
  
  data <- rbind(inland_data, offshore_data)
  
  return(data)
  
}

data <- buffer_func(100)



new_bind <- bind %>% 
  group_by(lon, lat) %>% 
  summarise()

points_in_buff <- function(df1, df2, df3){
  
  coordinates(df1) <- ~lon+lat
  
  coordinates(df2) <- ~lon+lat
  
  nn1 = get.knnx(coordinates(df2), coordinates(df1), 1)
  
  
  il = nn1$nn.dist[,1]
  
  df1$dist <- il
  
  new_data <- as.data.frame(df1)
  
  update <- new_data %>% 
    filter(dist <= 0.05)
  
  
  new_update <- update %>% 
    left_join(df3, by = c("lon", "lat")) %>% 
    na.omit() 
  
  return(new_update)
  
}


filt_data <- points_in_buff(new_bind, data, bind)
rm(data, bind,new_bind)
gc()

# First have to create 100 km buffer in order to stop edging effect on 50 km buffer.
# Then filter each month and save each monthly dataframe into .Rds file.


Dec <- filt_data %>% 
  filter(t == "12")

rm(Jun, Jul, Aug, Sep, Oct, Nov, Dec)
gc()

base::saveRDS(Dec, "C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/Dec.Rds")

# Calculate wind stress curl in python.
# Save dataframe with wind stress curl into csv.
# Read csv files into R, and combine the monthly dataframes into one dataframe.

Jan <- read.csv("C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/WSC_Jan.csv")

Jan$month <- "Jan"

bind <- rbind(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)
rm(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)
gc()

bind$month = factor(bind$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                           'Sep', 'Oct', 'Nov', 'Dec'))

base::saveRDS(filt_data, "C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/WSC.Rds")

# Then create a 50 km buffer, and perform a nearest neighbour analysis to filter
# out the coordinates that fall outside of the 50 km buffer zone.
# This eliminates the edging effect.

buffer_func <- function(buffer){
  
  # Download the coastline data for the entire world
  world_coastline <- ne_download(scale = "large", type = "coastline", category = "physical", returnclass = "sf")
  
  # Define the bounding box for Africa
  africa_bbox <- st_bbox(c(xmin = -20, ymin = -35, xmax = 52, ymax = 38), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box for Africa
  africa_coastline <- st_crop(world_coastline, africa_bbox)
  
  # The BCLME region:
  bclme_bbox <- st_bbox(c(xmin = 7.91755, ymin = -36.61979, xmax = 19.788742, ymax = -5.811113), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box
  bclme_coastline <- st_crop(africa_coastline, bclme_bbox)
  
  # Create a buffer of 50 nautical miles (1 nautical mile = 1852 meters)
  buffer_50nm <- st_make_valid(st_union(st_buffer(bclme_coastline, dist = buffer * 1852)))
  
  # Ensure coastline extends beyond buffer: Extend the coastline line if necessary
  # (I added/subtracted 3 degrees to the bbox to get a longer coastline)
  # Create a bbox slightly larger than the BCLME extent to get a longer coastline
  # that will insersect the buffer and thus can be used for splitting it lengthwise
  extended_bbox <- st_bbox(c(xmin = 7.91755 - 3,
                             ymin = -36.61979 - 3,
                             xmax = 19.788742 + 3,
                             ymax = -5.811113 + 3),
                           crs = st_crs(world_coastline))
  
  # Crop the coastline data to the outer bbox; this is the coastline that
  # will be used for splitting the buffer
  extended_coastline <- st_crop(africa_coastline, extended_bbox)
  
  # Convert the outer coastline to a single LINESTRING object
  extended_coastline_line <- st_union(st_cast(extended_coastline, "LINESTRING"))
  extended_coastline_line <- st_make_valid(extended_coastline_line)
  #extended_coastline_line_extended <- st_segmentize(extended_coastline_line, dfMaxLength = 1000000) # Adjust the max length as needed
  
  # Split the buffer using the extended coastline line
  split_buffers <- st_split(buffer_50nm, extended_coastline_line)
  
  # Extract the resulting polygons
  split_buffers_sf <- st_collection_extract(split_buffers)
  
  # Separate into two objects
  offshore_buffer <- split_buffers_sf[1, drop = FALSE]
  inland_buffer <- split_buffers_sf[2, drop = FALSE]
  
  # Ensure buffers are in the original CRS
  offshore_buffer <- st_transform(offshore_buffer, crs = st_crs(world_coastline))
  inland_buffer <- st_transform(inland_buffer, crs = st_crs(world_coastline))
  
  # Create simulated gridded data to test the buffer splitting
  # Define the bounding box for the region
  xmin <- 7.91755
  ymin <- -36.61979
  xmax <- 19.788742
  ymax <- -5.811113
  
  # Generate a grid of points within the bounding box
  lon <- seq(xmin, xmax, length.out = 200)
  lat <- seq(ymin, ymax, length.out = 400)
  grid <- expand.grid(lon = lon, lat = lat)
  
  # Add a simulated temperature column
  set.seed(42)  # For reproducibility
  grid$temperature <- runif(nrow(grid), min = 10, max = 30)
  
  # Convert the data frame to an sf object
  grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = st_crs(world_coastline))
  
  # Extract points within the inland buffer
  points_in_inland_buffer <- st_intersects(grid_sf, inland_buffer, sparse = FALSE)
  inland_data <- grid[apply(points_in_inland_buffer, 1, any), ]
  
  # Extract points within the offshore buffer
  points_in_offshore_buffer <- st_intersects(grid_sf, offshore_buffer, sparse = FALSE)
  offshore_data <- grid[apply(points_in_offshore_buffer, 1, any), ]
  
  
  data <- rbind(inland_data, offshore_data)
  
  return(data)
  
}

data <- buffer_func(50)



new_bind <- bind %>%
  rename(lon = longitude,
         lat = latitude) %>% 
  group_by(lon, lat) %>% 
  summarise()

points_in_buff <- function(df1, df2, df3){
  
  coordinates(df1) <- ~lon+lat
  
  coordinates(df2) <- ~lon+lat
  
  nn1 = get.knnx(coordinates(df2), coordinates(df1), 1)
  
  
  il = nn1$nn.dist[,1]
  
  df1$dist <- il
  
  new_data <- as.data.frame(df1)
  
  update <- new_data %>% 
    filter(dist <= 0.05)
  
  
  new_update <- update %>% 
    left_join(df3, by = c("lon", "lat")) %>% 
    na.omit() 
  
  return(new_update)
  
}


filt_data <- points_in_buff(new_bind, data, bind)



library(viridis)

ggplot(filt_data %>% 
         filter(month == "Jan"|
                  month == "Feb"|
                  month == "Dec") %>% 
         group_by(longitude, latitude) %>% 
         summarise(curl = mean(curl_tau)), aes(x = longitude, y = latitude)) +
  #coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  coord_fixed(xlim = c(14, 20), ylim = c(-35, -27)) +
  #coord_fixed(xlim = c(11, 17), ylim = c(-27, -17)) +
  geom_tile(aes(fill = curl)) +
  borders("world", regions = c("South Africa", "Namibia")) +
  scale_fill_viridis(option = "turbo") +
  #facet_wrap(~month) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Wind Stress Curl")


