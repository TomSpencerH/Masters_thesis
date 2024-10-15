library(tidyverse)
library(ncdf4)


files <- list.files(path = "C:/Users/tomsp/Earthdata/DMI/December/",pattern = "*.nc", full.names = TRUE)

  
  O = lapply(files, function(x){
  
  ncFile <- nc_open(x)
  
  LonIdx <- which( ncFile$dim$lon$vals > 11 & ncFile$dim$lon$vals < 20)
  LatIdx <- which( ncFile$dim$lat$vals > -35 & ncFile$dim$lat$vals < -17)
  
  MyVariable <- ncvar_get( ncFile, "analysed_sst")[ LonIdx, LatIdx]
  
  
  lon <- ncFile$dim$lon$val[LonIdx]
  lat <- ncFile$dim$lat$val[LatIdx]
  time <- "Dec"
  
  
  lonlattime <- as.matrix(expand.grid(lon,lat, time))
  
  sst <- as.vector(MyVariable)
  
  df1 <- data.frame(cbind(lonlattime, sst))
  
  return(df1)
  
}) 


df1 <- do.call(rbind, O)
rm(O, files)
gc()


clean_up <- function(df){
  
  df1 <- df %>% 
    na.omit()
  
  colnames(df1) <- c("lon","lat","month","temp")
  
  df1$lon <- as.numeric(df1$lon)
  df1$lat <- as.numeric(df1$lat)
  df1$temp <- as.numeric(df1$temp)
  
  
  df1$temp <- df1$temp - 272.15
  
  
  
  
  df1 <- df1 %>% 
    group_by(month, lon, lat) %>% 
    summarise(temp = mean(temp))
  
  return(df1)
  
}

Dec <- clean_up(df1)

rm(df1)
gc()

second <- rbind(Sep, Oct, Nov, Dec)
rm(Sep, Oct, Nov, Dec)
gc()

bind <- rbind(first, second)
rm(first, second)
gc()

base::saveRDS(bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/DMI_OI_Full.Rds"))

