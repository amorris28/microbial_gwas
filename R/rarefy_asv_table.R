library(tidyverse)
library(amorris)

# Import ASV table
asv_table <- read.csv('../output/gab_asv_table.csv')

comm_mat <- as.matrix(asv_table[, -1])

hist(rowSums(comm_mat))

samples <- tibble(Sample = asv_table$Sample)

rare_mat <- rarefy_mean(comm_mat, n = 1000)

rare_asv <- cbind(samples, as_tibble(rare_mat))

write_csv(rare_asv, '../output/gab_rare_asv_table.csv')

