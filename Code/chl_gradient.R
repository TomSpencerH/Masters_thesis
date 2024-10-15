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

new_bind <- chl %>% 
  group_by(lon, lat) %>% 
  summarise()


SA_bath <- getNOAA.bathy(lon1 = 11, lon2 = 20,
                         lat1 = -35, lat2 = -17, resolution = 0.6)


SB_bath <- fortify.bathy(SA_bath)

SB_bath <- SB_bath %>%
  rename(lon = x,
         lat = y,
         bathy = z) %>% 
  filter(bathy <= 0) 


bath <- SB_bath %>%
  filter(bathy >= -250) %>% 
  group_by(lat, lon) %>% 
  summarise(bathy = mean(bathy))

library(FNN)


bathy <- bath

coordinates(new_bind) <- ~lon+lat

coordinates(bath) <- ~lon+lat

nn1 = get.knnx(coordinates(bath), coordinates(new_bind), 1)


il = nn1$nn.dist[,1]

new_bind$dist <- il

new_data <- as.data.frame(new_bind)

update <- new_data %>% 
  filter(dist <= 0.04)



new_update <- update %>% 
  left_join(chl, by = c("lon", "lat")) %>% 
  na.omit()

new_update$month = factor(new_update$month, levels=c('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                                         'Sep', 'Oct', 'Nov', 'Dec'))

ggplot(new_update, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = chl)) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  #coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  coord_fixed(xlim = c(14, 20), ylim = c(-35, -27)) +
  scale_fill_viridis(option = "turbo") +
  facet_wrap(~month) +
  ggtitle("Summer") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill='Chlorophyll') 

# Seasonal (Summer vs Winter)

sum_chl <- new_update %>% 
  filter(month == "Dec"|
           month == "Jan"|
           month == "Feb") %>% 
  group_by(lon, lat) %>% 
  summarise(chl = mean(chl))


wint_chl <- new_update %>% 
  filter(month == "Jun"|
           month == "Jul"|
           month == "Aug") %>% 
  group_by(lon, lat) %>% 
  summarise(chl = mean(chl))

# Inter-annual

sbus <- new_update %>% 
  filter(between(lon, 14, 20),
         between(lat, -35, -27))

nbus <- new_update %>% 
  filter(between(lon, 11, 17),
         between(lat, -27, -17))



library(viridis)

# Seasonal

p1 <- ggplot(sum_chl, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = chl)) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  coord_fixed(xlim = c(11, 20), ylim = c(-35, -17)) +
  scale_fill_viridis(option = "turbo") +
  ggtitle("Summer") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill='Chlorophyll') 

ggarrange(p1, p2, ncol = 2)

# Inter-annual

ggplot(nbus, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = chl)) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  #coord_fixed(xlim = c(14, 20), ylim = c(-35, -27)) +
  coord_fixed(xlim = c(11, 16), ylim = c(-27, -17)) +
  scale_fill_viridis(option = "turbo") +
  ggtitle("NBUS") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill='Chlorophyll') +
  facet_wrap(~month)
  

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




