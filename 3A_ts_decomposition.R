
# Only pass in data for one site and one parameter at a time
# Needs the columns `station`, `param`, `value`, `date`
apply_seasonal_decomposition <- function(data) {
  
  # Prep the timeseries object
  ts_obj <- ts(data = data$value,  
               frequency=12)
  
  # Now decompose timeseries
  ts_decomposed_obj <- decompose(ts_obj)
  
  # Return a nice tibble with all of these
  tibble(station = unique(data$station),
         param = unique(data$param),
         observed = ts_decomposed_obj$x,
         trend = ts_decomposed_obj$trend,
         seasonal = ts_decomposed_obj$seasonal,
         random = ts_decomposed_obj$random) %>% 
    mutate(seasonalANDrandom = seasonal + random,
           year_frac = data$year_frac)
  
}

# For each station and parameter, apply time series decomposition
# and extract the resulting timeseries of the original (`observed`),
# the trend, the seasonal, and the random component.
decomposed_ts_all <- swmp_data_NAfilled %>% 
  split(.$station_param) %>% 
  map(~apply_seasonal_decomposition(.x)) %>% 
  bind_rows(.id = 'station_param') %>% 
  pivot_longer(cols = -c(station_param, station, param, year_frac), names_to = 'ts_type') %>% 
  # Join the ordered factor column `station_param_f` back in after reformatting from lists
  left_join(distinct(select(swmp_data_NAfilled, station_param, station_param_f)), by = "station_param")

# Plot everything
ts_components_figure_all <- decomposed_ts_all %>% 
  filter(!ts_type %in% c('seasonalANDrandom')) %>% 
  ggplot(aes(x = year_frac, y = value, color = ts_type)) +
  geom_line(linewidth=1) + 
  facet_wrap(~station_param_f, scales='free', ncol=2) +
  ylab("Parameter Value") + xlab("Year-Month") +
  ggtitle('Decomposed time series components') +
  scale_color_scico_d(name = "TS Component", end = 0.80, palette = 'glasgow') +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, face='bold'))

# Plot only the trend vs observed
ts_components_figure_obs_trend <- decomposed_ts_all %>% 
  filter(ts_type %in% c('observed', 'trend')) %>% 
  ggplot(aes(x = year_frac, y = value, 
             color = ts_type, alpha = ts_type, 
             linewidth = ts_type)) + 
  geom_line() + 
  scale_linewidth_manual(values = c(observed=0.5, trend=1.25), guide = "none") +
  scale_alpha_manual(values = c(observed=0.5, trend=1), guide = "none") +
  facet_wrap(~station_param_f, scales='free', ncol=2) +
  ylab("Parameter Value") + xlab("Year-Month") +
  ggtitle('Decomposed time series components: observed signal + trend') + 
  scale_color_scico_d(name = "TS Component", end = 0.80, palette = 'glasgow') +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, face='bold'))

# Plot only the seasonal
seasonal_col <- scico::scico(n=4, palette = 'glasgow', end = 0.80, categorical=TRUE)[3]
ts_components_figure_seasonal <- decomposed_ts_all %>% 
  filter(ts_type %in% c('seasonal')) %>% 
  ggplot(aes(x = year_frac, y = value)) +
  geom_line(linewidth = 1, color = seasonal_col) + 
  facet_wrap(~station_param_f, scales='free', ncol=2) +
  ylab("Parameter Value") + xlab("Year-Month") +
  ggtitle('Decomposed time series components: seasonal signal') + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, face='bold'))

# Plot seasonal + random to see how year-to-year changes over the ts
ts_components_figure_seasonalrandom <- decomposed_ts_all %>% 
  filter(ts_type %in% c('seasonal', 'seasonalANDrandom')) %>% 
  mutate(ts_typef = factor(ts_type, levels = c('seasonal', 'seasonalANDrandom'),
                           labels = c('Seasonal', 'Seasonal + Random'),
                           ordered = TRUE)) %>% 
  ggplot(aes(x = year_frac, y = value, color = ts_typef, linewidth = ts_typef)) +
  geom_line() + 
  facet_wrap(~station_param_f, scales='free', ncol=2) +
  scale_color_manual(values = c(Seasonal = seasonal_col, 
                                `Seasonal + Random` = 'grey30'),
                     name = 'TS Component') +
  scale_linewidth_manual(values = c(Seasonal = 1.25,  `Seasonal + Random` = 0.75),
                         guide = 'none') +
  ylab("Parameter Value") + xlab("Year-Month") +
  ggtitle('Decomposed time series components: seasonal signal') + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, face='bold'))

# Clean up environment since this is likely being sourced to load the data
rm(apply_seasonal_decomposition, decomposed_ts_all, seasonal_col)
