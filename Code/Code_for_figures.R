# Load required libraries
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(cowplot)
library(marmap)

# Load the Africa map
africa <- ne_countries(continent = "Africa", returnclass = "sf")

# Define bounding box for Benguela upwelling region
bbox <- st_as_sfc(st_bbox(c(xmin = 11, xmax = 20, ymin = -35, ymax = -17), crs = st_crs(africa)))

SA_bath <- getNOAA.bathy(lon1 = 11, lon2 = 20,
                         lat1 = -35, lat2 = -17, resolution = 0.6)


SB_bath <- fortify.bathy(SA_bath)

SB_bath <- SB_bath %>%
  rename(lon = x,
         lat = y,
         bathy = z) %>% 
  filter(bathy <= 0) 



SB_bath$lon <- round(SB_bath$lon, 2)
SB_bath$lat <- round(SB_bath$lat, 2)

bath <- SB_bath %>%
  filter(bathy >= -250) %>% 
  group_by(lat, lon) %>% 
  summarise(bathy = mean(bathy))

# Create the main map for the Benguela upwelling region
main_map <- ggplot() +
  geom_contour(data = SB_bath %>% 
                 filter(bathy < -250), aes(x = lon, y = lat, z = bathy), col = "cyan4") +
  geom_contour(data = bath, aes(x = lon, y = lat, z = bathy), col = "black") +
  geom_sf(data = africa, fill = "tan", color = "black") +  # Background of Africa
  coord_sf(xlim = c(11, 20), ylim = c(-35, -17), expand = FALSE) +  # Zoom into the Benguela region
  theme_minimal() + 
  annotation_scale(location = "bl", width_hint = 0.3) +  # Scale bar for main map
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering()) +
  xlab("Longitude") +
  ylab("Longitude")

# Create the inset map for the full continent of Africa
base_map <- ggplot() +
  geom_sf(data = africa, fill = "lightgray", color = "cyan4") +
  geom_sf(data = bbox, fill = NA, color = "red", linewidth = 1) +
  #coord_sf(xlim = c(-20, 60), ylim = c(-40, 40)) +  # View of Africa
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = "black"),
        panel.grid.minor = element_line(color = "black")) 

 

# Combine the main and inset maps
combined_map <- ggdraw() +
  draw_plot(base_map, 0.56, 0.57, 0.3, 0.5, scale = 0.7) +
  draw_plot(main_map, 0, 0, 1, 1)

# Print the final combined map
print(combined_map)


# Shore-Normal vs Cross-Shore Transects -----------------------------------


our_nc_data <- nc_open("Raw_data/Full.nc")

transects_func <- function(df){
  
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
    filter(between(lon, 11, 20),
           between(lat, -36, -17))
  
  #SA_west_coast <- df2 %>% 
  #filter(between(lon, 8, 16),
  #between(lat, -26, -18))
  
  
  
  SA_bath <- transects(SA_west_coast, spread = 65)
  
  
  
  heading2 <- data.frame(geosphere::destPoint(p = select(SA_bath, lon, lat),  
                                              b = SA_bath$heading, d = 200000))
  
  # Add the new coordinates tot he site list
  site_list <- SA_bath %>% 
    mutate(lon_dest = heading2$lon,
           lat_dest = heading2$lat)
  
  return(site_list)
  
}

site_list <- transects_func(our_nc_data)

interp <- function(rng, n) {
  seq(rng[1], rng[2], length = n)
}


fix_func <- function(df, x){
  
  df <- df[x, ]
  
  
  df2 <- data.frame(
    
    lon = c(df$lon, df$lon_dest),
    lat = c(df$lat, df$lat_dest)
  )
  
  
  
  munched <- data.frame(
    lon = round(interp(df2$lon, 41), 2),
    lat = round(interp(df2$lat, 41), 2),
    dist = c(0:40),
    group = x
  )
  
  
  return(munched)
  
}

y <- NULL


for(i in 1:nrow(site_list)){
  
  dfs <- fix_func(site_list, i)
  y <- rbind(y, dfs)
  
  
}



bathy <- function(res){
  
  SA_bath <- getNOAA.bathy(lon1 = 11, lon2 = 20,
                           lat1 = -35, lat2 = -17, resolution = res)
  
  
  SB_bath <- fortify.bathy(SA_bath)
  
  SB_bath <- SB_bath %>%
    rename(lon = x,
           lat = y,
           bathy = z) %>% 
    filter(bathy <= 0) 
  
  
  
  SB_bath$lon <- round(SB_bath$lon, 2)
  SB_bath$lat <- round(SB_bath$lat, 2)
  
  bath <- SB_bath %>%
    filter(bathy >= -250) %>% 
    group_by(lat, lon) %>% 
    summarise(bathy = mean(bathy))
  
  return(bath)
  
}

bath <- bathy(res = 0.6)

final_transect <- function(df1, df2){
  
  new_bathy <- df1 %>% 
    group_by(lat) %>% 
    mutate(min_lon = min(lon))
  
  
  
  d1 <- data.frame(lon = df2$lon,
                    lat = df2$lat)
  
  d2 <- data.frame(lon = new_bathy$min_lon,
                    lat = new_bathy$lat)
  
  
  
  era <- dplyr::intersect(d1, d2)
  
  
  
  
  connect <- era %>% 
    left_join(df2, by = c("lon", "lat")) %>% 
    na.omit()
  
  
  y2 <- connect %>% 
    left_join(df2, by = c("group")) %>% 
    select(lon.y, lat.y, group, dist.y, dist.x) %>% 
    rename(lon = lon.y,
           lat = lat.y,
           dist = dist.y,
           max.dist = dist.x)
  
  return(y2)
  
  
}

transects <- final_transect(df1 = bath, df2 = y)
  
  transects <- transects %>% 
    group_by(group) %>% 
    filter(dist == 0:max.dist) %>% 
    filter(row_number() %% 2 != 0) %>% 
    slice_sample(n = 6) %>% 
    filter(min(dist) <= 9)
  






coastline <- ne_countries(scale = "medium", returnclass = "sf", continent = "Africa")

# Define the bounding box for the Southern African west coast
bounding_box <- st_bbox(c(xmin = 10, xmax = 20.5, ymin = -35, ymax = -15), crs = st_crs(4326))

sf_use_s2(FALSE)

# Crop the coastline to the bounding box
sa_coastline <- st_crop(coastline, bounding_box)


p1 <- ggplot() +
  geom_contour(data = bath, aes(x = lon, y = lat, z = bathy), col = "cyan4") +
  geom_line(data = transects, aes(x = lon, y = lat, group = group)) +
  geom_point(data = transects, aes(x = lon, y = lat), col = "red3") +
  geom_sf(data = sa_coastline, fill = "grey80", color = "black") +
  coord_sf(xlim = c(12, 20), ylim = c(-35, -27)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Shore-normal")

lon <- c(17.17, 12.5, 17.65, 12.8, 18.25, 13.5, 17.85, 12.9)
lat <- c(-30, -30.6, -31, -31.5, -32, -32, -33, -33.1)
point <- c("A", "A", "B", "B", "C", "C", "D", "D")  

cross <- data.frame(lon, lat, point)

p2 <- ggplot() +
  geom_contour(data = bath, aes(x = lon, y = lat, z = bathy), col = "cyan4") +
  geom_sf(data = sa_coastline, fill = "grey80", color = "black") +
  coord_sf(xlim = c(12, 20), ylim = c(-35, -27)) +
  geom_line(data = cross, aes(x = lon, y = lat, group = point)) +
  geom_point(data = cross, aes(x = lon, y = lat), size = 4, col = "red3") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Cross-shore")

ggarrange(p1, p2, ncol = 2, align = "h")



