library(ncdf4)
library(tidyverse)
library(ggpubr)
library(coastR)
library(marmap)
library(lattice)
library(latticeExtra)
library(RColorBrewer)



# Copernicus Marine Data Store:
# 
# Online ------------------------------------------------------------------
setwd('/Users/tomsp/Desktop/Masters/Practice_thesis/data/')

our_nc_data <- nc_open("CHL/chl_2017.nc")
print(our_nc_data)

extract_nc <- function(df){
  
  lat <- ncvar_get(our_nc_data, "latitude")
  
  lon <- ncvar_get(our_nc_data, "longitude")
  
  time <- ncvar_get(our_nc_data, "time")
  
  #get the variable in "matrix slices"
  chl <- ncvar_get(our_nc_data, "CHL") 
  
  
  fillvalue_chl <- ncatt_get(our_nc_data, "CHL", "_FillValue")
  
  
  #right away let's replace the nc FillValues with NAs
  chl[chl==fillvalue_chl$value] <- NA
  
  
  #convert the hours into date + hour
  #as_datetime() function of the lubridate package needs seconds
  timestamp <- as_datetime(c(time),origin="1970-01-01", tz = "GMT")
  
  new_timestamp <- format(timestamp, "%b")
  
  lonlattime <- as.matrix(expand.grid(lon,lat, new_timestamp))
  
  
  
  chl_vec_long <- as.vector(chl)
  
  
  df1 <- data.frame(cbind(lonlattime, chl_vec_long))
  
  
  
  
  colnames(df1) <- c("lon","lat","month","chl")
  
  df1$chl <- as.numeric(df1$chl)
  df1$lon <- as.numeric(df1$lon)
  df1$lat <- as.numeric(df1$lat)
  
  df1[df1 == "NaN"] <- NA
  
  
  
  df1 <- df1 %>% 
    na.omit() %>% 
    group_by(month, lon, lat) %>% 
    summarise(chl = mean(chl))
  
  return(df1)
  
}

chl_2021 <- extract_nc(our_nc_data)

rm(our_nc_data)
gc()



bind <- rbind(chl_2017, chl_2018, chl_2019, chl_2020, chl_2021)

rm(chl_2017, chl_2018, chl_2019, chl_2020, chl_2021)
gc()

new_bind <- bind %>% 
  group_by(month, lon, lat) %>% 
  summarise(chl = mean(chl))

rm(bind)
gc()


base::saveRDS(new_bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/chl_climatology.Rds"))


setwd("C:/Users/tomsp/Desktop/Masters/Practice_thesis/Masters/")

chl <- readRDS("Processed_data/Chl-a/chl_climatology.Rds")


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


chl$lon <- round(chl$lon, 2)
chl$lat <- round(chl$lat, 2)


final_transect <- function(df1, df2, df3) {
  
  
  new_bathy <- df3 %>% 
    group_by(lat) %>% 
    mutate(min_lon = min(lon))
  
  
  
  d1 <- data.frame(lon = df1$lon,
                   lat = df1$lat)
  
  d2 <- data.frame(lon = new_bathy$min_lon,
                   lat = new_bathy$lat)
  
  
  
  era <- dplyr::intersect(d1, d2)
  
  
  
  
  connect <- era %>% 
    left_join(y, by = c("lon", "lat")) %>% 
    na.omit()
  
  
  y2 <- connect %>% 
    left_join(df1, by = c("group")) %>% 
    select(lon.y, lat.y, group, dist.y, dist.x) %>% 
    rename(lon = lon.y,
           lat = lat.y,
           dist = dist.y,
           max.dist = dist.x)
  
  
  n2 <- y2 %>% 
    group_by(group) %>% 
    filter(dist == 0:max.dist)
  
  
  
  n3 <- n2 %>% 
    left_join(df2, by = c("lon", "lat")) %>%
    na.omit()
  
  
  models <- n3 %>% 
    group_by(group, month) %>% 
    group_modify(~ broom::tidy(lm(chl ~ dist, data = .x))) %>%
    na.omit() %>%  
    filter(term == "dist") %>% 
    select(group, month, estimate, p.value) %>% 
    rename(slope = estimate)
 
  
  slope_df <- n3 %>% 
    left_join(models, by = c("group", "month")) %>% 
    na.omit()
  
  slope <- slope_df %>% 
    group_by(group) %>% 
    mutate(p75 = max.dist*0.75) %>% 
    filter(min(dist) <= 9|
             max(dist) > p75) 
  
  return(slope)
  
}


slope <- final_transect(df1 = y, df2 = chl, df3 = bath)


slope$month = factor(slope$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                         'Sep', 'Oct', 'Nov', 'Dec'))

coastline <- ne_countries(scale = "medium", returnclass = "sf", continent = "Africa")

# Define the bounding box for the Southern African west coast
bounding_box <- st_bbox(c(xmin = 10, xmax = 20.5, ymin = -35, ymax = -15), crs = st_crs(4326))

sf_use_s2(FALSE)

# Crop the coastline to the bounding box
sa_coastline <- st_crop(coastline, bounding_box)

library(viridis)
library(colorspace)

sum_slope <- slope %>% 
  filter(month == "Dec"|
           month == "Jan"|
           month == "Feb")



wint_slope <- slope %>% 
  filter(month == "Jun"|
           month == "Jul"|
           month == "Aug")


range(p.value$slope)

ggplot() +
  geom_contour(data = bath, aes(x = lon, y = lat, z = bathy), colour = "black", alpha = 0.6) +
  geom_point(data = p.value, aes(x = lon, y = lat, col = slope)) +
  geom_line(data = p.value, aes(x = lon, y = lat, group = group, col = slope), linewidth = 1) +
  geom_sf(data = sa_coastline, fill = "grey80", color = "black") +
  #coord_sf(xlim = c(17, 20), ylim = c(-35, -32)) +
  scale_color_continuous_diverging(palette = "Green-brown", l1 = 20, l2 = 90, p1 = 0.7, p2 = 1,
                                   rev = FALSE,
                                   limits = c(-0.3, 0.3),
                                   # guide = "colourbar",
                                   guide = guide_legend(even.steps = FALSE,
                                                        show.limits = TRUE)) +
  #geom_text(data = area, aes(x = lon, y = lat, label = Area)) +
  #facet_wrap(~month) +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(color='CHL Gradient') 

  

hovmoller <- new_update %>% 
  group_by(lat, month) %>% 
  summarise(chl = mean(chl))

# Hovmoller Plot:


chlTheme <- custom.theme(region= brewer.pal(n=10, 'BrBG'))

myTheme <- custom.theme(region=rev(brewer.pal(n=10, 'RdBu')))

hovmoller$month1 <- match(hovmoller$month, month.abb)


 levelplot(chl ~ month1*lat,
          data=hovmoller,
          xlab='Month', ylab='Lat',
          par.settings=myTheme,
          aspect = 2)




