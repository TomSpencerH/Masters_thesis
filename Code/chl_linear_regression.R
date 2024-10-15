library(tidyverse)
library(tidync)
library(marmap)
library(FNN)
library(sp)
library(sf)
library(zoo)

data <- tidync("Raw_data/CHL/chl_2024.nc")

data$transforms$time$time <- as.POSIXct((data$transforms$time$time), origin = "1970-01-01")

  
  chl_2024 <- data %>%
    hyper_tibble(select_var = c("CHL"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(t = time, lon = longitude, lat = latitude, chl = CHL) %>% 
    dplyr::group_by(year(t), month(t), lon, lat) %>% 
    dplyr:: summarise(chl = mean(chl))
  

  
  bind2 <- rbind(chl_2017, chl_2018, chl_2019, chl_2020, chl_2021, chl_2022, chl_2023, chl_2024)
  rm(chl_2017, chl_2018, chl_2019, chl_2020, chl_2021, chl_2022, chl_2023, chl_2024)
  gc()
  
  rm(chl_dat, chl_dat1, chl_dat2, chl_dat3, chl_dat4, chl_dat5, chl_dat6, chl_dat7, chl_dat8)
  gc()

  

colnames(bind1) <- c("year", "month", "lon", "lat", "chl")

bind <- rbind(bind1, bind2)
rm(bind1, bind2)
gc()

complete_data <- transform(bind, season = as.yearqtr(as.yearmon(paste(year, month, sep = "-")) + 1/12))

saveRDS(complete_data, "Processed_data/chl_complete.Rds")


bind <- complete_data %>% 
  #filter(season != "2007 Q1" & season != "2024 Q1") %>% 
  filter(str_detect(season, 'Q3')) %>% 
  group_by(season, lon, lat) %>% 
  summarise(chl = mean(chl))

rm(complete_data)
gc()

saveRDS(bind, "Processed_data/chl_sum.Rds")
saveRDS(bind, "Processed_data/chl_wint.Rds")

bind <- readRDS("Processed_data/chl_sum.Rds")
bind <- readRDS("Processed_data/chl_wint.Rds")


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
  group_modify(~ broom::tidy(lm(chl ~ season, data = .x))) %>% 
  filter(term != "(Intercept)") %>% 
  select(lon, lat, estimate, p.value) %>% 
  rename(slope = estimate)




slope_data <- sst_bath %>% 
  left_join(models, by = c("lon", "lat"))


library(colorspace)

library(viridis)
 ggplot(slope_data, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = slope * 10)) +
  #geom_contour(data = bath, aes(z = bathy), colour = "black", alpha = 0.3, bins = 4) +
  metR::geom_contour_fill(data = slope_data %>% 
                            filter(p.value < 0.05),
                           aes(x = lon, y = lat,
                               z = slope), bins = 60) +
   scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 20, l2 = 100, p1 = 0.7, p2 = 1,
                                   rev = FALSE,
                                   n_interp = 31,
                                   # guide = "colourbar",
                                   guide = guide_legend(title = "Decadal rate of change",
                                                        even.steps = FALSE,
                                                        show.limits = TRUE)) +
  borders("world", regions = c("South Africa", "Namibia"), fill = "grey") +
  coord_fixed(xlim = c(11, 21), ylim = c(-35, -17)) +
   #scale_fill_gradient2(low = "blue3", high = "darkgreen") +
  ggtitle("Chlorophyll Concentration (DJF: 2000 - 2024)") +
  xlab("Longitude") +
  ylab("Latitude") +
   labs(fill = "Decadal rate of Change")
