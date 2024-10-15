library(tidyverse)

# The information for the NOAA OISST data
rerddap::info(datasetid = "jplMURSST41mday", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# This function downloads and prepares data based on user provided start and end dates
OISST_sub_dl <- function(time_df){
  OISST_dat <- rerddap::griddap(datasetx = "jplMURSST41mday",
                                url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                                time = c(time_df$start, time_df$end),
                                latitude = c(-35, -17),
                                longitude = c(11, 20),
                                fields = c("sst"))$data %>% 
    dplyr::mutate(time = base::as.Date(stringr::str_remove(time, "T12:00:00Z"))) %>% 
    dplyr::rename(t = time, temp = sst, lon = longitude, lat = latitude) %>% 
    dplyr::select(lon, lat, t, temp) %>% 
    stats::na.omit()
}


# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:5,
                       start = as.Date(c("2021-01-01", "2021-02-01", 
                                         "2021-03-01", "2021-04-01", "2021-05-01")),
                       end = as.Date(c("2021-01-31", "2021-02-28", 
                                       "2021-03-31", "2021-04-30", "2021-06-30")))

# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:5,
                       start = as.Date(c("2021-07-01", "2021-08-01", 
                                         "2021-09-01", "2021-10-01", "2021-11-01")),
                       end = as.Date(c("2021-07-31", "2021-08-31", 
                                       "2021-09-30", "2021-10-31", "2021-12-31")))


doParallel::registerDoParallel(cores = 2)

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed
base::system.time(
  MUR_data <- dl_years %>% 
    dplyr::group_by(date_index) %>% 
    dplyr::group_modify(~OISST_sub_dl(.x)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(lon, lat, t, temp)
) 

ten <- MUR_data %>% 
  group_by(t, lon, lat) %>% 
  summarise(temp = mean(temp))


rm(MUR_data)
gc()

bind <- rbind(one, two, three, four, five, six, seven, eight, nine, ten)
rm(one, two, three, four, five, six, seven, eight, nine, ten)
gc()

new_bind <- bind %>% 
  mutate(month = month(t))

rm(bind)
gc()

bind <- new_bind %>% 
  group_by(month, lon, lat) %>% 
  summarise(temp = mean(temp))

rm(new_bind)
gc()

base::saveRDS(bind, file = ("C:/Users/tomsp/Desktop/Masters/Practice_thesis/data/MUR_full.Rds"))
