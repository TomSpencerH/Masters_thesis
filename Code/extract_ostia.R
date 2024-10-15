library(tidyverse)
library(ncdf4)
library(tidync)



library(tidyverse)
library(tidync)

data <- tidync("Raw_data/OSTIA/ostia_2007.nc")

data$transforms$time$time <- as.POSIXct((data$transforms$time$time), origin = "1970-01-01")

sst_load <- function(df, lon1, lon2, lat1, lat2, time_start, time_end) {
  wind_dat <- df %>%
    hyper_filter(lon = between(lon, lon1, lon2),
                 lat = between(lat, lat1, lat2),
                 time = time >= time_start & time <= time_end) %>%
    hyper_tibble(select_var = c("analysed_sst"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(t = time, lon = longitude, lat = latitude, temp = analysed_sst) %>% 
    dplyr::group_by(year(t), month(t), lon, lat) %>% 
    dplyr:: summarise(temp = mean(temp)) 
  
  return(wind_dat)
  
}

load_sst <- sst_load(data, lon1 = 14, lon2 = 20, lat1 = -35, lat2 = -27, time_start = "2007-01-01", time_end = "2011-12-31")





base::saveRDS(new_bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/OSTIA_full.Rds"))
