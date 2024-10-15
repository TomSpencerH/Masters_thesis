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

setwd("C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/")


bind <- readRDS("Processed_data/SST/OSPO_full.Rds")



our_nc_data <- nc_open("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/Dist2Coast/Full.nc")

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

y <- NULL


for(i in 1:nrow(site_list)){
  
  dfs <- fix_func(site_list, i)
  y <- rbind(y, dfs)
  
  
}

bind$lon <- round(bind$lon, 2)
bind$lat <- round(bind$lat, 2)



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
  left_join(bind, by = c("lon", "lat")) %>%
  na.omit()

n3 <- n3 %>% 
  group_by(group) %>% 
  filter(n()>2)

models <- n3 %>% 
  group_by(group, month) %>% 
  group_modify(~ broom::tidy(lm(temp ~ dist, data = .x))) %>%
  na.omit() %>% 
  filter(term == "dist") %>% 
  select(group, month, estimate) %>% 
  rename(slope = estimate)


slope_df <- n3 %>% 
  left_join(models, by = c("group", "month")) %>% 
  na.omit()

slope <- slope_df %>% 
  group_by(group) %>% 
  filter(min(dist) <= 9) 



 #slope$month <- month.abb[slope$month]



slope$month = factor(slope$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                           'Sep', 'Oct', 'Nov', 'Dec'))

sbus <- slope %>% 
  filter(between(lat, -35, -27))


 #slope <- slope %>% 
  #group_by(group) %>% 
  #slice_sample(n = 50)

rm(bind, connect, n3, df, df1, df2, dfs, era, models, new_bathy, our_nc_data, SB_bath, site_list, slope_df, y, y2, i, SA_bath)
gc()

library(viridis)
library(colorspace)

ggplot(slope %>% 
         filter(month == "Dec"|
                  month == "Jan"|
                  month == "Feb"), aes(x = lon, y = lat)) +
  geom_contour(data = bath, aes(z = bathy), colour = "black", alpha = 0.6) +
  geom_point(aes(col = slope)) +
  geom_line(aes(group = group, col = slope), linewidth = 1) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  #coord_quickmap(expand = F) +
  coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  scale_color_viridis(option = "turbo") +
  #geom_text(data = area, aes(x = lon, y = lat, label = Area)) +
  #facet_wrap(~month) +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(color='SST Gradient') +
  ggtitle("OSTIA")

ost_sbus <- ggplot(sbus, aes(x = lon, y = lat)) +
  geom_contour(data = bath, aes(z = bathy), colour = "black", alpha = 0.6) +
  geom_point(aes(col = slope)) +
  geom_line(aes(group = group, col = slope), linewidth = 1) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  #coord_quickmap(expand = F) +
  coord_fixed(xlim = c(14, 20), ylim = c(-35, -27)) +
  scale_color_viridis(option = "turbo") +
  geom_text(data = area, aes(x = lon, y = lat, label = Area)) +
  facet_wrap(~month) +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(color='SST Gradient') +
  ggtitle("OSTIA")

mur_wint <- ggplot(slope_mur2, aes(x = lon, y = lat)) +
  geom_contour(data = bath, aes(z = bathy), colour = "black", alpha = 0.6) +
  geom_point(aes(col = slope)) +
  geom_line(aes(group = group, col = slope), linewidth = 1) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  #coord_quickmap(expand = F) +
  coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  scale_color_viridis(option = "turbo") +
  geom_text(data = area, aes(x = lon, y = lat, label = Area)) +
  #facet_wrap(~month) +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(color='SST Gradient') +
  ggtitle("MUR")

ggarrange(osp_wint, ost_wint, dmi_wint, mur_wint, ncol = 4, common.legend = TRUE)

ost_hovmoller <- slope %>% 
  group_by(lat, month) %>% 
  summarise(slope = mean(slope))

# Hovmoller Plot:

myTheme <- custom.theme(region=rev(brewer.pal(n=10, 'RdBu')))


dmi_hovmoller$month1 <- match(dmi_hovmoller$month, month.abb)

mur_hov <- levelplot(slope ~ month1*lat,
                data=mur_hovmoller,
                xlab='Month', ylab='Latitude',
                par.settings=myTheme,
                main = "MUR")

ggarrange(ospo_hov, ostia_hov, dmi_hov, mur_hov, ncol = 4)


