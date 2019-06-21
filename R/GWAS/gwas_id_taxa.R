library(tidyverse)
library(xtable)
library(amorris)
source('functions.R')

# import model output
temp <- list.files(path="../output", pattern="fits*", full.names = TRUE)
all_models <- sapply(temp, readRDS, simplify = FALSE, USE.NAMES = TRUE)
names(all_models) <- reverse_substr(substr(names(all_models), 16, 27), 5, 12)

taxon_table <- read_csv('../output/taxon_table.csv')

# Pull out sig ASVs
asvs <- lapply(all_models, id_asvs)

# Number of significant ASVs
lapply(asvs, length)

# Pull out taxonomy for sig ASVs
sig_taxa <- lapply(asvs, id_taxonomy, taxon_table) 

# Save taxa
sapply(names(sig_taxa), function (x) write.table(sig_taxa[[x]], file = paste0('../tables/sig_taxa_', x, ".txt")))
