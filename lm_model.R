#### Input: adj_data.csv
#### Output: model_fits.rds

library(tidyverse)
library(purrr)
library(broom)
library(modelr)

adj_data <- read_csv('output/adj_data.csv')
com_adj_data <- read_csv('output/com_adj_data.csv')
env_adj_data <- read_csv('output/env_adj_data.csv')
spa_adj_data <- read_csv('output/spa_adj_data.csv')
troph_total <- read_csv('output/troph_total.csv')



##### Functions

lowk_model <- function(df) {
  lm(Low_final_k ~ abund, data = df)
}
vmax_model <- function(df) {
  lm(Vmax ~ abund, data = df)
}
# Create nested data set by asv
  # Calculate linear model for each asv between abundance and function
fit_models <- function(data, model) {
  data %>% 
    gather(asv, abund, starts_with('asv_')) %>% 
    group_by(asv) %>% 
    nest() %>% 
    mutate(model = map(data, model)) %>% 
    return()
}

saveRDS(fit_models(adj_data, lowk_model), file = "output/lowk_fits.rds")
saveRDS(fit_models(com_adj_data, lowk_model), file = "output/lowk_fits_com.rds")
saveRDS(fit_models(env_adj_data, lowk_model), file = "output/lowk_fits_env.rds")
saveRDS(fit_models(spa_adj_data, lowk_model), file = "output/lowk_fits_spa.rds")
saveRDS(fit_models(troph_total, lowk_model), file = "output/lowk_fits_raw.rds")
saveRDS(fit_models(adj_data, vmax_model), file = "output/vmax_fits.rds")
saveRDS(fit_models(com_adj_data, vmax_model), file = "output/vmax_fits_com.rds")
saveRDS(fit_models(env_adj_data, vmax_model), file = "output/vmax_fits_env.rds")
saveRDS(fit_models(spa_adj_data, vmax_model), file = "output/vmax_fits_spa.rds")
saveRDS(fit_models(troph_total, vmax_model), file = "output/vmax_fits_raw.rds")
