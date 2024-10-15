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
           t = format(t, "%m")) %>% 
    dplyr::filter(Hour == "02") %>% 
    group_by(t, lon, lat) %>% 
    summarise(u10 = mean(u10),
              v10 = mean(v10))
  
  return(wind_dat)

}

load4 <- wind_load(data, lon1 = 11, lon2 = 22, lat1 = -35, lat2 = -17, time_start = "2021-01-01", time_end = "2021-12-31")

wind_data <- rbind(load, load1, load2, load3, load4)
rm(load, load1, load2, load3, load4)
gc()

bind <- wind_data %>% 
  group_by(t, lon, lat) %>% 
  summarise(u10 = mean(u10),
            v10 = mean(v10))

rm(data, wind_data)
gc()

saveRDS(bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/Francois_wind_02h00.Rds"))

bind <- readRDS("C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/Processed_data/Wind/Francois_wind_14h00.Rds")

wind <- bind %>% 
  mutate(wind_spd = sqrt(u10^2 + v10^2),
         CD = ifelse(wind_spd <= 1, 0.00218,
                     ifelse(wind_spd > 1 & wind_spd < 3, (0.62 + 1.56/ wind_spd)*0.001,
                            ifelse(wind_spd > 3 & wind_spd < 10, 0.00114, (0.49 + 0.065*wind_spd)*0.001))), 
         taux = 1.293*CD*wind_spd*u10,
         tauy = 1.293*CD*wind_spd*v10,
         curl = -(lat/abs(lat))*(taux*cos(160-90) + tauy*sin(160-90)))

new_wind <- wind %>% 
  mutate(wind_dir_trig_to = atan2(u10/wind_spd, v10/wind_spd) ,
         wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi,
         wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180)


rm(bind)
gc()

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


new_bind <- wind %>% 
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

filt_data$t <- as.numeric(filt_data$t)

filt_data$t <- month.abb[filt_data$t]


filt_data$t = factor(filt_data$t, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                                   'Sep', 'Oct', 'Nov', 'Dec'))

ggplot(filt_data, aes(x = lon, y = lat)) +
  coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  geom_tile(aes(fill = curl)) +
  borders("world", regions = c("South Africa", "Namibia")) +
  scale_fill_viridis(option = "turbo") +
  facet_wrap(~t) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Wind Stress Curl")




# Seasonal (Summer vs Winter)

sum_chl <- filt_data %>% 
  filter(t == "Dec"|
           t == "Jan"|
           t == "Feb") %>% 
  group_by(lon, lat) %>% 
  summarise(curl = mean(curl),
            wind_spd = mean(wind_spd))


wint_chl <- filt_data %>% 
  filter(t == "Jun"|
           t == "Jul"|
           t == "Aug") %>% 
  group_by(lon, lat) %>% 
  summarise(curl = mean(curl),
            wind_spd = mean(wind_spd))

# Inter-annual

sbus <- filt_data %>% 
  filter(between(lon, 14, 20),
         between(lat, -35, -27))

nbus <- filt_data %>% 
  filter(between(lon, 11, 17),
         between(lat, -27, -17))

sbus$month = factor(sbus$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                         'Sep', 'Oct', 'Nov', 'Dec'))
nbus$month = factor(nbus$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                         'Sep', 'Oct', 'Nov', 'Dec'))

library(viridis)

# Seasonal

p1 <- ggplot(sum_chl, aes(x = lon, y = lat)) +
  coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  geom_tile(aes(fill = curl)) +
  borders("world", regions = c("South Africa", "Namibia")) +
  scale_fill_viridis(option = "turbo") +
  ggtitle("Summer") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill='Wind Stress Curl') 

ggarrange(p1, p2, ncol = 2)

# Inter-annual

ggplot(nbus, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = curl)) +
  borders("world", regions = c("South Africa", "Namibia")) +
  #coord_fixed(xlim = c(14, 20), ylim = c(-35, -27)) +
  coord_fixed(xlim = c(11, 16), ylim = c(-27, -17)) +
  #coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  scale_fill_viridis(option = "turbo") +
  ggtitle("SBUS") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill='Wind Stress Curl') +
  facet_wrap(~t)


hovmoller <- filt_data %>% 
  group_by(lat, t) %>% 
  summarise(curl = mean(curl))

# Hovmoller Plot:


myTheme <- custom.theme(region=rev(brewer.pal(n=10, 'RdBu')))

hovmoller$month1 <- match(hovmoller$t, month.abb)


levelplot(curl ~ month1*lat,
          data=hovmoller,
          xlab='Month', ylab='Lat',
          par.settings=myTheme,
          aspect = 2)

