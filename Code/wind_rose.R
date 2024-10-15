library(circular)


bind <- readRDS("Processed_data/wind_data_14h.Rds")

bind$wind_dir <- atan2(-bind$u, -bind$v) * (180 / pi)
bind$wind_dir <- (270 - bind$wind_dir) %% 360


# Bin wind speed data for categorization
bind <- bind %>%
  mutate(wind_spd_binned = cut(wind_spd, 
                               breaks = c(0, 2, 4, 6, 8, 10, 12, 14, Inf), 
                               labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12", "12-14", "14+")))

# Aggregate the data by wind direction and wind speed category
agg_data <- bind %>%
  group_by(month, wind_dir = round(wind_dir/10)*10, wind_spd_binned) %>%
  summarise(count = n()) %>%
  ungroup()

# Function to create a wind rose plot with facet_wrap
create_facet_wind_rose <- function(data) {
  ggplot(data, aes(x = factor(wind_dir), y = count, fill = wind_spd_binned)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(start = 0) +
    scale_x_discrete(drop = FALSE, breaks = seq(0, 360, by = 30)) +
    scale_fill_discrete_sequential(palette = "Batlow", rev = FALSE,
                                   name = "Wind Speed (m/s)") +
    facet_wrap(~ month, ncol = 4, scales = "free") +
    labs(x = "Wind Direction (degrees)", y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

library(colorspace)

facet_wind_rose <- create_facet_wind_rose(agg_data)
facet_wind_rose  



    