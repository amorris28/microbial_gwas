library(tidyverse)
library(purrr)
library(broom)
library(modelr)
library(amorris)
source('functions.R')

# Import all model output
temp <- list.files(path="../output", pattern="gab_adj*", full.names = TRUE)
adj_data <- sapply(temp, read.table, simplify = FALSE, USE.NAMES = TRUE)
# Rename data subsets
names(adj_data) <- reverse_substr(substr(names(adj_data), 19, 25), 5, 7)
# Fit models and save model outputs
lapply(names(adj_data), function(x) saveRDS(fit_models(adj_data[[x]], lowk_model), 
                                     file = paste0("../output/fits_lowk_", x, 
                                                   ".rds")))
lapply(names(adj_data), function(x) saveRDS(fit_models(adj_data[[x]], vmax_model), 
                                     file = paste0("../output/fits_vmax_", x, 
                                                   ".rds")))
