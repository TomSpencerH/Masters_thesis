library(tidyverse)
library(ncdf4)
library(marmap)
library(data.table)
library(coastR)
library(e1071) # C-means Fuzzy clustering for larger datasets
library(lattice)
library(grDevices)
library(ggpubr)
library(metR) 
library(zoo)
library(lattice)
library(latticeExtra)
library(RColorBrewer)

# Load Data: --------------------------------------------------------------

setwd("C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/")

data <- readRDS("Processed_data/SST/OSPO_full.Rds")

# SBUS

sbus <- data %>% 
  filter(between(lat, -35, -27))

nbus <- data %>% 
  filter(between(lat, -27, -17))




our_nc_data <- nc_open("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/Dist2Coast/NBUS.nc")

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
  
  
  #SA_west_coast <- df2 %>% 
    #filter(between(lon, 14, 20),
           #between(lat, -36, -27))
  
  SA_west_coast <- df2 %>% 
    filter(between(lon, 8, 16),
           between(lat, -26, -18))
  
  
  
  SA_bath <- transects(SA_west_coast, spread = 30)
  
  
  
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
    lon = round(interp(df2$lon, 201), 2),
    lat = round(interp(df2$lat, 201), 2),
    dist = c(0:200),
    group = x
  )
  
  
  return(munched)
  
}

fix_func(site_list, 1)

y <- NULL


for(i in 1:nrow(site_list)){
  
  dfs <- fix_func(site_list, i)
  y <- rbind(y, dfs)
  
  
}

sbus$lon <- round(sbus$lon, 2)
sbus$lat <- round(sbus$lat, 2)

nbus$lon <- round(nbus$lon, 2)
nbus$lat <- round(nbus$lat, 2)


SA_bath <- getNOAA.bathy(lon1 = 14, lon2 = 20,
                         lat1 = -35, lat2 = -27, resolution = 0.6)

SA_bath <- getNOAA.bathy(lon1 = 8, lon2 = 16,
                         lat1 = -26, lat2 = -18, resolution = 0.6)


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

new_bathy <- bath %>% 
  group_by(lat) %>% 
  mutate(min_lon = min(lon))



df1 <- data.frame(lon = y$lon,
                  lat = y$lat)

df2 <- data.frame(lon = new_bathy$min_lon,
                  lat = new_bathy$lat)

era <- dplyr::intersect(df1, df2)



connect <- era %>% 
  left_join(y, by = c("lon", "lat")) %>% 
  na.omit()

y2 <- connect %>% 
  left_join(y, by = c("group")) %>% 
  select(lon.y, lat.y, group, dist.y, dist.x) %>% 
  rename(lon = lon.y,
         lat = lat.y,
         dist = dist.y,
         max.dist = dist.x)


n2 <- y2 %>% 
  group_by(group) %>% 
  filter(dist == 0:max.dist)


n3 <- n2 %>% 
  left_join(nbus, by = c("lon", "lat")) %>%
  na.omit()

n3 <- n3 %>% 
  group_by(group) %>% 
  filter(n()>2)



library(plyr)

models <- dlply(n3, c("group", "month"), function(df) 
  lm(temp ~ dist, data = df))

# Apply coef to each model and return a data frame
df <- ldply(models, coef)

if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE) 

colnames(df) <- c("group", "month", "intercept", 'slope')

slope_df <- n3 %>% 
  left_join(df, by = c("group", "month")) 

slope <- slope_df %>% 
  group_by(group) %>% 
  filter(min(dist) <= 9)

#slope$month <- month.abb[slope$month]



slope$month = factor(slope$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                                   'Sep', 'Oct', 'Nov', 'Dec'))

#slope_mur <- slope %>% 
#group_by(group) %>% 
#slice_sample(n = 50)



rm(slope_df, n1, connect, df, df1, df2, dfs, dist, era, lonlattime, models, new, new_bathy, new_data, new_y, new.df, our_nc_data, y, y2, dist_vec, i, lat, lon)
gc()

library(viridis)

ggplot(slope, aes(x = lon, y = lat)) +
  geom_contour(data = bath, aes(z = bathy), colour = "black", alpha = 0.6) +
  geom_point(aes(col = slope)) +
  geom_line(aes(group = group, col = slope), linewidth = 0.1) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  #coord_fixed(xlim = c(14, 20), ylim = c(-35, -27)) +
  coord_fixed(xlim = c(11, 17), ylim = c(-27, -17)) +
  geom_text(data = area, aes(x = lon, y = lat, label = Area), size = 2) +
  scale_color_viridis(option = "turbo") +
  scale_fill_viridis_c() +
  ggtitle("OSPO") +
  facet_wrap(~month)

ggarrange(p1, p2, p3, ncol = 3, legend = "none")



hovmoller1 <- slope1 %>% 
  group_by(lat, month) %>% 
  summarise(slope = mean(slope))

# Hovmoller Plot:

myTheme <- custom.theme(region=rev(brewer.pal(n=10, 'RdBu')))

p2 <- levelplot(slope ~ month*lat,
          data=hovmoller1,
          xlab='Month', ylab='Lat',
          par.settings=myTheme)

library(gridExtra)

grid.arrange(p1,p4, p3, p2, ncol = 4)
