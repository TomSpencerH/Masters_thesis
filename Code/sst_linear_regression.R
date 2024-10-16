# Linear regression for SST products:
# First product: OSTIA
# 

library(tidyverse)
library(tidync)
library(marmap)
library(FNN)
library(sp)
library(sf)
library(zoo)
library(metR)

data <- tidync("Raw_data/OSTIA/ostia_2024.nc")

data$transforms$time$time <- as.POSIXct((data$transforms$time$time), origin = "1970-01-01")

files <- list.files(path = "Raw_data/OSTIA/",pattern = "*.nc", full.names = TRUE)

O = lapply(files, function(x){
  
  data <- tidync(x)
  
  data$transforms$time$time <- as.POSIXct((data$transforms$time$time), origin = "1970-01-01")

  wind_dat <- data %>%
    hyper_tibble(select_var = c("analysed_sst"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(t = time, lon = longitude, lat = latitude, temp = analysed_sst) %>% 
    dplyr::group_by(year(t), month(t), lon, lat) %>% 
    dplyr:: summarise(temp = mean(temp) - 273.15)
    
  
  return(wind_dat)
  
})

bind <- do.call(rbind, O)  
rm(O, files)  
gc()


colnames(bind) <- c("year", "month", "lon", "lat", "temp")



saveRDS(bind, "Processed_data/ostia_data.Rds")

bind <- readRDS("Processed_data/ostia_data.Rds")



complete_data <- transform(bind, season = as.yearqtr(as.yearmon(paste(year, month, sep = "-")) + 1/12))


bind <- complete_data %>% 
  #filter(season != "2007 Q1") %>% 
  filter(str_detect(season, 'Q1')) %>% 
  group_by(season, lon, lat) %>% 
  summarise(temp = mean(temp))

rm(complete_data)
gc()

str(bind$season)

#bind <- readRDS("Processed_data/mur_sum_data.Rds")

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

new_sst <- bind %>% 
  group_by(lon, lat) %>% 
  summarise()

points_in_buff <- function(df1, df2, df3){
  
  coordinates(df1) <- ~lon+lat
  
  coordinates(df2) <- ~lon+lat
  
  nn1 = get.knnx(coordinates(df2), coordinates(df1), 1)
  
  
  il = nn1$nn.dist[,1]
  
  df1$dist <- il
  
  new_data <- as.data.frame(df1)
  
  update <- new_data %>% 
    filter(dist <= 0.05)
  
  
  new_update <- update %>% 
    left_join(df3, by = c("lon", "lat")) %>% 
    na.omit() %>% 
    select(-dist)
  
  return(new_update)
  
}

sst_bath <- points_in_buff(new_sst, bath, bind)


  
  models <- sst_bath %>% 
    group_by(lon, lat) %>% 
    group_modify(~ broom::tidy(lm(temp ~ season, data = .x))) %>% 
    filter(term != "(Intercept)") %>% 
    select(lon, lat, estimate, p.value) %>% 
    rename(slope = estimate)
  
  
  
  
  slope_data <- sst_bath %>% 
    left_join(models, by = c("lon", "lat"))
  

rm(models, sst_bath, new_sst, SA_bath)
gc()

coastline <- ne_countries(scale = "medium", returnclass = "sf", continent = "Africa")

# Define the bounding box for the Southern African west coast
bounding_box <- st_bbox(c(xmin = 10, xmax = 20.5, ymin = -35, ymax = -15), crs = st_crs(4326))

# Crop the coastline to the bounding box
sa_coastline <- st_crop(coastline, bounding_box)
  
  library(colorspace)
  p1 <- ggplot(slope_data, aes(x = lon, y = lat)) +
    metR::geom_contour_fill(aes(z = slope*10)) +
    geom_contour(data = slope_data %>% 
                         filter(p.value < 0.05), 
                 aes(x = lon, y = lat, z = slope), bins = 4, col = "black") +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 20, l2 = 90, p1 = 0.7, p2 = 1,
                                    rev = FALSE,
                                    limits = c(-0.4, 0.4),
                                    # guide = "colourbar",
                                    guide = guide_legend(title = "Decadal rate of change",
                                                         even.steps = FALSE,
                                                         show.limits = TRUE)) +
    geom_sf(data = sa_coastline, fill = "grey80", color = "black") +
    #scale_fill_viridis(option = "turbo") +
    ggtitle("Sea Surface Temperature (DJF: 2002 - 2024)") +
    xlab("Longitude") +
    ylab("Latitude")
  
  library(ggmagnify)
  p1 + geom_magnify(from = c(17.6, 19, -35, -33.5), 
               to = c(xmin = 11.5, xmax = 14, ymin = -33, ymax = -30.5)) # OSPO + DMI?
 

ggarrange(p1, p2, p3, ncol = 3)

# Data product: MUR

# The information for the NOAA OISST data
rerddap::info(datasetid = "jplMURSST41mday", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# This function downloads and prepares data based on user provided start and end dates
OISST_sub_dl <- function(time_df){
  OISST_dat <- rerddap::griddap(datasetx = "jplMURSST41mday",
                                url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                                time = c(time_df$start, time_df$end),
                                latitude = c(-35, -27),
                                longitude = c(14, 20),
                                fields = c("sst"))$data %>% 
    dplyr::mutate(time = base::as.Date(stringr::str_remove(time, "T12:00:00Z"))) %>% 
    dplyr::rename(t = time, temp = sst, lon = longitude, lat = latitude) %>% 
    dplyr::select(lon, lat, t, temp) %>% 
    stats::na.omit() %>% 
    dplyr::mutate(year = year(t),
           month = month(t)) %>% 
    dplyr::group_by(year, month, lon, lat) %>% 
    dplyr::summarise(temp = mean(temp))
}


# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:4,
                       start = as.Date(c("2022-06-16", "2022-12-16", 
                                         "2023-06-16", "2023-12-16")),
                       end = as.Date(c("2022-11-16", "2023-05-16", 
                                       "2023-11-16", "2024-02-16")))




doParallel::registerDoParallel(cores = 2)

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed
base::system.time(
  MUR_data6 <- dl_years %>% 
    dplyr::group_by(date_index) %>% 
    dplyr::group_modify(~OISST_sub_dl(.x)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(lon, lat, year, month, temp) 
) 

gc()



rm(one)
gc()



complete_data <- transform(bind, season = as.yearqtr(as.yearmon(paste(year, month, sep = "-")) + 1/12))

rm(MUR_data1, MUR_data2, MUR_data3, MUR_data4, MUR_data5, MUR_data6)
gc()

# Data Product: OSPO

files <- list.files(path = "Raw_data/OSPO/2024/",pattern = "*.nc", full.names = TRUE)

O = lapply(files, function(x){
  
  data <- tidync(x)
  
  data$transforms$time$time <- as.POSIXct((data$transforms$time$time), origin = "1981-01-01")
  
  wind_dat <- data %>%
    hyper_filter(lon = between(lon, 11, 20),
                 lat = between(lat, -35, -17)) %>% 
    hyper_tibble(select_var = c("analysed_sst"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(t = time, lon = lon, lat = lat, temp = analysed_sst) %>% 
    dplyr::group_by(year(t), month(t), lon, lat) %>% 
    dplyr:: summarise(temp = mean(temp) - 273.15) %>% 
    dplyr::distinct(.keep_all = TRUE) %>% 
    stats::na.omit()
  
  
  return(wind_dat)
  
})

bind_2024 <- do.call(rbind, O)  
rm(O, files)  
gc()


bind2 <- rbind(bind, bind_2020, bind_2021, bind_2022, bind_2023, bind_2024)
rm(bind, bind_2020, bind_2021, bind_2022, bind_2023, bind_2024)
gc()


colnames(bind2) <- c("year", "month", "lon", "lat", "temp")

bind2 <- bind2 %>% 
  group_by(year, month, lon, lat) %>% 
  summarise(temp = mean(temp))

bind <- rbind(bind1, bind2)
rm(bind1, bind2)
gc()

saveRDS(bind, "Processed_data/ospo_data.Rds")



# Data Product: DMI


files <- list.files(path = "Raw_data/DMI/2023_Q3/",pattern = "*.nc", full.names = TRUE)



O = lapply(files, function(x){
  
  data <- tidync(x)
  
  data$transforms$time$time <- as.numeric("2023")
  
  
  wind_dat <- data %>%
    hyper_filter(lon = between(lon, 11, 20),
                 lat = between(lat, -35, -17)) %>% 
    hyper_tibble(select_var = c("analysed_sst"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(season = time, lon = lon, lat = lat, temp = analysed_sst) %>% 
    dplyr::group_by(season, lon, lat) %>% 
    dplyr:: summarise(temp = mean(temp) - 273.15) %>% 
    dplyr::distinct(.keep_all = TRUE) %>% 
    stats::na.omit()
  
  
  return(wind_dat)
  
})

Q1_2023 <- do.call(rbind, O)  

rm(O, files)  
gc()


bind2 <- rbind(Q1_2019, Q1_2020, Q1_2021, Q1_2022, Q1_2023)
rm(Q1_2019, Q1_2020, Q1_2021, Q1_2022, Q1_2023)

bind <- rbind(bind1, bind2)
rm(bind1, bind2)
gc()

bind <- bind %>% 
  group_by(season, lon, lat) %>% 
  summarise(temp = mean(temp))


saveRDS(bind, "Processed_data/dmi_data_wint.Rds")

data <- readRDS("Processed_data/dmi_data_sum.Rds")


data <- data %>% 
  mutate(season = as.numeric(substring(t, 1, 4))) %>% 
  ungroup() %>% 
  select(-t)

saveRDS(data, "Processed_data/dmi_data_sum.Rds")
  