library(ncdf4)
library(tidync) # tidync is unfortunately not loading on the computer with more RAM for some reason
library(tidyverse)
library(reshape2)
library(lubridate)
library(stringr)
library(circular)
library(doParallel); doParallel::registerDoParallel(cores = 2) 
library(heatwaveR)
library(ggthemes)
library(ggpubr)
library(data.table)



# Online ------------------------------------------------------------------
setwd('/Users/tomsp/ERA5/')

our_nc_data <- nc_open("C:/Users/tomsp/ERA5/download_2.nc")
print(our_nc_data)

tidy_nc <- function(our_nc_data){
  
  
  lat <- ncvar_get(our_nc_data, "latitude")
  
  lon <- ncvar_get(our_nc_data, "longitude")
  
  time <- ncvar_get(our_nc_data, "time")
  
  
  #get the variable in "matrix slices"
  u10 <- ncvar_get(our_nc_data, "u10") 
  v10 <- ncvar_get(our_nc_data, "v10")
  
  
  fillvalue_u <- ncatt_get(our_nc_data, "u10", "_FillValue")
  fillvalue_v <- ncatt_get(our_nc_data, "v10", "_FillValue")
  
  
  #right away let's replace the nc FillValues with NAs
  u10[u10==fillvalue_u$value] <- NA
  v10[v10==fillvalue_v$value] <- NA
  
  
  #convert the hours into date + hour
  #as_datetime() function of the lubridate package needs seconds
  timestamp <- as_datetime(c(time*60*60),origin="1900-01-01", tz = "GMT")
  new_timestamp <- format(timestamp, "%b")
  
  lonlattime <- as.matrix(expand.grid(lon,lat,new_timestamp))
  
  u10_vec_long <- as.vector(u10)
  v10_vec_long <- as.vector(v10)
  
  
  
  df1 <- data.frame(cbind(lonlattime, u10_vec_long, v10_vec_long))
  
  colnames(df1) <- c("lon","lat","month","u10", "v10")
  
  df1$u10 <- as.numeric(df1$u10)
  df1$v10 <- as.numeric(df1$v10)
  df1$lon <- as.numeric(df1$lon)
  df1$lat <- as.numeric(df1$lat)
  
  df1 <- df1 %>% 
    group_by(month, lon, lat) %>% 
    summarise(u10 = mean(u10),
              v10 = mean(v10)) %>% 
    na.omit()
  
  return(df1)
  
}

  second <- tidy_nc(our_nc_data)
  rm(our_nc_data)
  gc()
  

bind <- rbind(first, second)

rm(first, second)
gc()

bind <- bind %>% 
  na.omit() %>% 
  group_by(month, lon, lat) %>% 
  summarise(u10 = mean(u10),
            v10 = mean(v10))
gc()

base::saveRDS(bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/ERA5_Land_10km.Rds"))



