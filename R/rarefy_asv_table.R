library(tidyverse)
source('R/Rarefy_mean.R')
source('R/misc_functions.R')

# Import ASV table
asv_table <- read_tsv('data/16S_dada2_table_RDP_tax_edit.txt', n_max = 100000)

# Pull out 
taxon_table <- 
  asv_table %>% 
  select(OTU_ID, ConsensusLineage) %>% 
  mutate(OTU_ID = paste('asv', OTU_ID, sep = "_")) %>% 
  select(asv = OTU_ID, ConsensusLineage) %>% 
  separate(ConsensusLineage, c('Domain', 'Phylum', 'Class', 'Order', 'Family',
                               'Genus', 'Species'), ";") %>% 
  select(-Species)
write_csv(taxon_table, 'output/taxon_table.csv')

# Rename ASVs to remove backticks in column names, add Sample ID column
# Convert all abundances to doubles
wide_asv <- asv_table %>% 
  mutate(asv = paste('asv', OTU_ID, sep = "_")) %>% 
  select(asv, starts_with('dna')) %>% 
  select(asv, ends_with('MO')) %>% 
  gather(Sample, abund, 2:length(.)) %>% 
  spread(asv, abund) %>% 
  mutate(Sample = substring(Sample, 5))
wide_asv <- wide_asv[, colSums(wide_asv != 0, na.rm = T) > 0]

comm_mat <- as.matrix(wide_asv[, -1])

samples <- tibble(Sample = wide_asv$Sample)

rare_mat <- rarefy.mean(comm_mat, n = 1000)

rare_asv <- cbind(samples, as_tibble(rare_mat))

wide_asv <- cbind(samples, as_tibble(comm_mat))


write_csv(wide_asv, 'output/asv_table.csv')

write_csv(rare_asv, 'output/rare_asv_table.csv')

