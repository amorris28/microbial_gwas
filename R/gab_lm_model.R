library(tidyverse)
library(purrr)
library(broom)
library(modelr)
library(amorris)
source('functions.R')


adj_data_c <- read.table('../output/gab_adj_com.tsv')
adj_data_e <- read.table('../output/gab_adj_env.tsv')
adj_data_s <- read.table('../output/gab_adj_spa.tsv')
adj_data_ce <- read.table('../output/gab_adj_com_env.tsv')
adj_data_cs <- read.table('../output/gab_adj_com_spa.tsv')
adj_data_es <- read.table('../output/gab_adj_env_spa.tsv')
adj_data_ces <- read.table('../output/gab_adj_com_env_spa.tsv')

adj_data <- list(c = adj_data_c, e = adj_data_e, s = adj_data_s, 
                 ce = adj_data_ce, cs = adj_data_cs, es = adj_data_es,
                 ces = adj_data_ces)
lapply(names(adj_data), function(x) saveRDS(fit_models(adj_data[[x]], lowk_model), 
                                     file = paste0("../output/lowk_fits_", x, 
                                                   ".rds")))
saveRDS(fit_models(adj_data_c, lowk_model), file = "../output/lowk_fits_c.rds")
saveRDS(fit_models(com_adj_data, lowk_model), file = "../output/lowk_fits_com.rds")
saveRDS(fit_models(env_adj_data, lowk_model), file = "../output/lowk_fits_env.rds")
saveRDS(fit_models(spa_adj_data, lowk_model), file = "../output/lowk_fits_spa.rds")
saveRDS(fit_models(troph_total, lowk_model), file = "../output/lowk_fits_raw.rds")
saveRDS(fit_models(adj_data, vmax_model), file = "../output/vmax_fits.rds")
saveRDS(fit_models(com_adj_data, vmax_model), file = "output/vmax_fits_com.rds")
saveRDS(fit_models(env_adj_data, vmax_model), file = "output/vmax_fits_env.rds")
saveRDS(fit_models(spa_adj_data, vmax_model), file = "output/vmax_fits_spa.rds")
saveRDS(fit_models(troph_total, vmax_model), file = "output/vmax_fits_raw.rds")
