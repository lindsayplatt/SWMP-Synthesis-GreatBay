# Borrowed these code concepts for using `wsyn` from Danny Szydlowski

# Wilkinson, G. M., Walter, J., Fleck, R., & Pace, M. L. (2020). 
# Beyond the trends: The need to understand multiannual dynamics in aquatic ecosystems. 
# Limnology and Oceanography Letters, 5(4), 281–286. https://doi.org/10.1002/lol2.10153

# Run wavelets separately for each site (both will keep the precip data)
for(this_station in c('gb', 'sq')) {
  
  # Prepare nice names for plot titles
  this_station_name <- swmp_data_NAfilled %>% 
    filter(station == this_station) %>% 
    pull(stationf) %>% unique()
  weather_station_name <- swmp_data_NAfilled %>% 
    filter(station == 'gb-gl') %>% 
    pull(stationf) %>% unique()
  
  # Prepare data for calculating wavelets
  swmp_wavelet_data <- swmp_data_NAfilled %>%
    filter(station %in% c(this_station, 'gb-gl')) %>% 
    select(year, month, year_frac, param, value) %>% 
    pivot_wider(id_cols = c(year, month, year_frac), names_from = 'param', values_from = 'value') %>% 
    # Fill NAs for missing data
    mutate(temp_median = zoo::na.approx(temp_median, na.rm = FALSE),
           turb_median = zoo::na.approx(turb_median, na.rm = FALSE),
           no23f = zoo::na.approx(no23f, na.rm = FALSE)) %>% 
    # Filter out the first three values which cannot be linearly interpolated to fill NAs
    filter(!is.na(temp_median) & !is.na(turb_median) & !is.na(no23f))
  
  # Clean the data to fit wavelet assumptions
  no23f <- cleandat(swmp_wavelet_data$no23f, 1:length(swmp_wavelet_data$no23f), clev=5)$cdat
  precip <- cleandat(swmp_wavelet_data$totprcp_total, 1:length(swmp_wavelet_data$totprcp_total), clev=5)$cdat
  turb <- cleandat(swmp_wavelet_data$turb_median, 1:length(swmp_wavelet_data$turb_median), clev=5)$cdat
  temp <- cleandat(swmp_wavelet_data$temp_median, 1:length(swmp_wavelet_data$temp_median), clev=5)$cdat
  
  # we can then apply the wavelet transform to break the data down into oscillating
  # patterns at different timescales
  wt.no23f <- wt(no23f, 1:length(no23f))
  wt.precip <- wt(precip, 1:length(precip))
  wt.turb <- wt(turb, 1:length(turb))
  wt.temp <- wt(temp, 1:length(temp))
  
  # arrange them side by side
  par(mfrow=c(1,4), mai=c(0.3,0.3,0.3,0.3), omi=c(0.4,0.4,0.2,0.1))
  
  plotmag(wt.precip, title = sprintf("%s, Precipitation (mm)", weather_station_name))
  plotmag(wt.no23f, title = sprintf("%s, Nitrate/Nitrite (mg/L)", this_station_name))
  plotmag(wt.turb, title = sprintf("%s, Turbidity (NTU)", this_station_name))
  plotmag(wt.temp, title = sprintf("%s, Temperature (deg C)", this_station_name))
}

# Remove all objects from environment except for the ones from previous runs
rm(list=c(ls()[!grepl('swmp_|figure', ls())], "swmp_wavelet_data"))
