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



files <- list.files(path = "C:/Users/tomsp/Earthdata/2022_2/",pattern = "*.nc", full.names = TRUE)

O = lapply(files, function(x){
  
  ncFile <- nc_open(x)
  
    LonIdx <- which( ncFile$dim$lon$vals > 11 & ncFile$dim$lon$vals < 20)
    LatIdx <- which( ncFile$dim$lat$vals > -35 & ncFile$dim$lat$vals < -30)
    #LonIdx <- which( ncFile$dim$lon$vals > 11 & ncFile$dim$lon$vals < 20)
    #LatIdx <- which( ncFile$dim$lat$vals > -30 & ncFile$dim$lat$vals < -25)
  #LonIdx <- which( ncFile$dim$lon$vals > 11 & ncFile$dim$lon$vals < 20)
  #LatIdx <- which( ncFile$dim$lat$vals > -25 & ncFile$dim$lat$vals < -17)
  MyVariable <- ncvar_get( ncFile, "analysed_sst")[ LonIdx, LatIdx]
 
  
  lon <- ncFile$dim$lon$val[LonIdx]
  lat <- ncFile$dim$lat$val[LatIdx]
  time <- ncFile$dim$time$val
  
  
  timestamp <- as_datetime(c(time),origin="1981-01-01", tz = "GMT")
  
  lonlattime <- as.matrix(expand.grid(lon,lat,timestamp))
  
  sst <- as.vector(MyVariable)
  
  df1 <- data.frame(cbind(lonlattime, sst))
  
  return(df1)
  
}) 

df7 <- do.call(rbind, O)  
rm(O, files)  
gc()

clean_up <- function(df){
  
  df1 <- df %>% 
    na.omit()
  
  colnames(df1) <- c("lon","lat","date","sst")
  
  df1$lon <- as.numeric(df1$lon)
  df1$lat <- as.numeric(df1$lat)
  df1$sst <- as.numeric(df1$sst)
  df1$date <- as.Date(df1$date)
  
  df1 <- df1 %>% 
    mutate(month = month(date))
  
  df1$sst <- df1$sst - 272.15
  
  
  
  
  df1 <- df1 %>% 
    group_by(month, lon, lat) %>% 
    summarise(temp = mean(sst))
  
  return(df1)
  
}
  
df7 <- clean_up(df7)
  
gc()


new_bind <- rbind(df1, df2, df3, df4, df5, df6, df7)
rm(df1, df2, df3, df4, df5, df6, df7)
gc()

new_bind <- new_bind %>% 
  group_by(month, lon, lat) %>% 
  summarise(sst = mean(sst))

bind <- rbind(bind, new_bind)
rm(new_bind)

bind <- bind %>% 
  group_by(month, lon, lat) %>% 
  summarise(sst = mean(sst)) %>% 
  filter(lon >= 10)

base::saveRDS(bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/k10_NBUS.Rds"))

gc()

ggplot(bind, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = sst)) +
  xlim(8, 16) +
  ylim(-26, -17) +
  scale_fill_viridis_c() +
  facet_wrap(~month)


