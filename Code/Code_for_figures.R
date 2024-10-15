# Load required libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(cowplot)

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
  coord_sf(xlim = c(-20, 60), ylim = c(-40, 40)) +  # View of Africa
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = "black"),
        panel.grid.minor = element_line(color = "black"))

 

# Combine the main and inset maps
combined_map <- ggdraw() +
  draw_plot(base_map, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(main_map, 0, 0, 1, 1)

# Print the final combined map
print(combined_map)
