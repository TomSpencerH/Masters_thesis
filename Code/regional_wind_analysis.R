library(tidyverse)
library(ncdf4)
library(tidync)
library(CFtime)
library(ggpubr)
library(coastR)
library(marmap)
library(lattice)
library(latticeExtra)
library(RColorBrewer)
library(viridis)
library(FNN)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lwgeom) # For st_split function
library(nngeo)
library(terra)

our_nc_data <- nc_open("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/Dist2Coast/Full.nc")


lat <- ncvar_get(our_nc_data, "latitude")

lon <- ncvar_get(our_nc_data, "longitude")

#get the variable in "matrix slices"
dist <- ncvar_get(our_nc_data) 

lonlattime <- as.matrix(expand.grid(lon,lat))

dist_vec <- as.vector(dist)


df1 <- data.frame(cbind(lonlattime, dist_vec))

df2 <- df1 %>% 
  filter_all(any_vars(. %in% c('1')))

colnames(df2) <- c("lon", "lat", "dist")


SA_west_coast <- df2 %>% 
  filter(between(lon, 16, 18.8),
         between(lat, -34.5, -33.9))

ggplot(SA_west_coast, aes(x = lon, y = lat)) +
  borders("world", region = "South Africa") +
  coord_fixed(xlim = c(14, 20), ylim = c(-35, -30)) +
  geom_point()

rm(data, df1, df2, dist, lonlattime, our_nc_data, SA_west_coast, dist_vec)
gc()

# Using points along the coastline to determine the coordinates we want to use.

data <- tidync("C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/Raw_data/Wind/ccam_8km_197901_202304.nc")


data$transforms$time$time <- as.POSIXct((data$transforms$time$time)*60, origin = "1979-01-01")

wind_load <- function(df, lon1, lon2, lat1, lat2, time_start, time_end) {
  wind_dat <- df %>%
    hyper_filter(lon = between(lon, lon1, lon2),
                 lat = between(lat, lat1, lat2),
                 time = time >= time_start & time <= time_end) %>%
    hyper_tibble(select_var = c("uas", "vas"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(t = time, u10 = uas, v10 = vas) %>% 
    dplyr::mutate(Hour = format(t, "%H"),
                  t = as.Date(t)) %>% 
    dplyr::filter(Hour == "14") %>% 
    dplyr::select(-Hour)
  
  return(wind_dat)
  
}

load <- wind_load(data, lon1 = 16, lon2 = 22, lat1 = -34.5, lat2 = -33.9, time_start = "2002-06-16", time_end = "2013-06-16")

wind_data <- rbind(load, load1)
rm(load, load1)
gc()


buffer_func <- function(buffer){
  
  # Download the coastline data for the entire world
  world_coastline <- ne_download(scale = "large", type = "coastline", category = "physical", returnclass = "sf")
  
  # Define the bounding box for Africa
  africa_bbox <- st_bbox(c(xmin = -20, ymin = -35, xmax = 52, ymax = 38), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box for Africa
  africa_coastline <- st_crop(world_coastline, africa_bbox)
  
  # The BCLME region:
  bclme_bbox <- st_bbox(c(xmin = 7.91755, ymin = -36.61979, xmax = 21.788742, ymax = -5.811113), crs = st_crs(world_coastline))
  
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
                             xmax = 21.788742 + 3,
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
  xmax <- 21.788742
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


new_bind <- wind_data %>% 
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


filt_data <- points_in_buff(new_bind, data, wind_data)

ggplot(filt_data %>% 
         filter(t == "2016-01-01"), aes(x = lon, y = lat)) +
  borders("world", region = "South Africa") +
  coord_fixed(xlim = c(14, 22), ylim = c(-36, -30)) +
  geom_tile(aes(fill = u10))



saveRDS(filt_data, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/wind_for_py.Rds"))

new_wind <- read.csv("WSC_cc.csv")

new_wind <- new_wind %>% 
  rename(lon = longitude,
         lat = latitude,
         u10 = u_component,
         v10 = v_component,
         t = time)

wind_data <- readRDS("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/wind_for_py.Rds")

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


new_bind <- new_wind %>% 
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


filt_data <- points_in_buff(new_bind, data, new_wind)

ggplot(filt_data %>% 
         filter(t == "11854"), aes(x = lon, y = lat)) +
  borders("world", region = "South Africa") +
  coord_fixed(xlim = c(16, 21), ylim = c(-35.5, -33)) +
  geom_tile(aes(fill = u10))

new_data <- filt_data %>% 
  left_join(wind_data, by = c("lon", "lat", "u10", "v10")) %>% 
  na.omit() %>% 
  select(-t.x, -dist.x, -dist.y) %>% 
  rename(t = t.y)

rm(new_bind, wind_data, new_wind, filt_data, data)
gc()

new <- new_data[1:300000,]

sbux_sf <- st_as_sf(new, 
                    coords = c("lon", "lat"),  
                    crs = 4326)

model <- new %>% 
  group_by(lon, lat, t) %>% 
  mutate(model = lm(curl_tau ~ t, data = .,))

install.packages("bmstdr")
library(bmstdr)

?Bspatial

model <- Bspatial(model="lm", formula=curl_tau~t, data=new_data, 
         coordtype="lonlat", coords=c("lon", "lat"))
  
new_data$fitted <- model[["fitteds"]]
new_data$residuals <- model[["residuals"]]

ggplot(new_data %>% 
         filter(t == "2003-04-09"), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = residuals)) +
  borders("world", regions = "South Africa") +
  coord_fixed(xlim = c(16, 20), ylim = c(-34.7, -33.8)) +
  scale_fill_viridis_c()


library(plyr)

models <- dlply(new_data, c("lon","lat"), function(df) 
  lm(curl_tau ~ t, data = df))

# Apply coef to each model and return a data frame
df <- ldply(models, coef)

if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE) 

colnames(df) <- c("lon", "lat", "intercept", 'slope')

ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = slope)) +
  borders("world", regions = "South Africa") +
  coord_fixed(xlim = c(16, 20), ylim = c(-34.7, -33.8)) +
  scale_fill_viridis_c()

slope_df <- new_data %>% 
  left_join(df, by = c("lon", "lat", "t")) 

slope <- slope_df %>% 
  group_by(group) %>% 
  filter(min(dist) <= 9) 

plot <- new_data %>% 
  mutate(month = month(t)) %>% 
  filter(month == "1"|
           month == "2"|
           month == "12") %>% 
  mutate(p5 = quantile(curl_tau, probs = 0.0025)) %>% 
  filter(curl_tau <= p5)

ggplot(plot, aes(x = lon, y = lat)) +
  borders("world", region = "South Africa") +
  coord_fixed(xlim = c(16, 21), ylim = c(-35.5, -33)) +
  geom_tile(aes(fill = curl_tau)) +
  scale_fill_viridis(option = "turbo")

sum_data <- plot %>% 
  left_join(new_data, by = c("lon", "lat")) %>% 
  na.omit() %>% 
  select(lon, lat, t.y, curl_tau.y) %>% 
  rename(curl = curl_tau.y,
         t = t.y) %>% 
  distinct(.keep_all = TRUE)

ggplot(sum_data %>% 
         filter(t == "2002-06-20"), aes(x = lon, y = lat)) +
  borders("world", region = "South Africa") +
  coord_fixed(xlim = c(16, 21), ylim = c(-35.5, -33)) +
  geom_tile(aes(fill = curl)) +
  scale_fill_viridis(option = "turbo")


ggplot(sum_data %>% 
         mutate(year = year(t)) %>%
         filter(between(year, 2003, 2022)) %>% 
         group_by(year) %>% 
         summarise(curl = mean(curl)), aes(x = year, y = curl)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_smooth()
