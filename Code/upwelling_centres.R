our_nc_data <- nc_open("Raw_data/Full.nc")

  
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
  
frio <- SA_west_coast %>% 
  filter(between(lat, -19.2, -18.8)) %>% 
  mutate(Area = "Cape Frio")

walvis <- SA_west_coast %>% 
  filter(between(lat, -23.2, -22.8)) %>% 
  mutate(Area = "Walvis Bay")

luderitz <- SA_west_coast %>% 
  filter(lat == -26.5) %>% 
  mutate(Area = "Luderitz")

namaqua <- SA_west_coast %>% 
  filter(between(lat, -30.2, -29.8)) %>% 
  mutate(Area = "Namaqua")

cc <- SA_west_coast %>% 
  filter(lat == -32) %>% 
  mutate(Area = "Cape Columbine")

cape_town <- SA_west_coast %>% 
  filter(between(lat, -34.2, -33.8)) %>% 
  mutate(Area = "Cape Town")

area_data <- rbind(frio, walvis, cc, cape_town, namaqua, luderitz)
rm(frio, walvis, cc, cape_town, namaqua, luderitz, df1, df2, lonlattime, dist, our_nc_data, SA_west_coast, lat, lon, dist_vec)
gc()

area <- area_data %>% 
  group_by(Area) %>% 
  summarise(lon = mean(lon),
            lat = mean(lat))

area1 <- area$Area

ggplot(area, aes(x=lon, y=lat)) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  coord_fixed(xlim = c(10, 20), ylim = c(-35, -17)) +
  geom_point() +
  geom_label(
    label=area1,
    nudge_x=0.45, nudge_y=0.1,
    label.padding=unit(0.2, "lines"),
    color="white",
    fill="darkcyan"
  )

