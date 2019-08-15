library(tidyverse)
library(morris)

# Import ASV table
asv_table <- read_tsv('../raw_data/16S_dada2_table_RDP_tax_edit.txt', n_max = 100000)

# Pull out taxonomy
taxon_table <- 
  asv_table %>% 
  select(OTU_ID, ConsensusLineage) %>% 
  mutate(OTU_ID = paste('asv', OTU_ID, sep = "_")) %>% 
  select(asv = OTU_ID, ConsensusLineage) %>% 
  separate(ConsensusLineage, c('Domain', 'Phylum', 'Class', 'Order', 'Family',
                               'Genus', 'Species'), ";") %>% 
  select(-Species)
write_csv(taxon_table, '../output/taxon_table.csv')

# Rename ASVs to remove backticks in column names, add Sample ID column
# Convert all abundances to doubles
wide_asv <- asv_table %>% 
  mutate(asv = paste('asv', OTU_ID, sep = "_")) %>% 
  select(asv, starts_with('dna')) %>% 
  gather(Sample, abund, 2:length(.)) %>% 
  spread(asv, abund) %>% 
  mutate(Sample = substring(Sample, 5))

# Remove zero abundance taxa and NA taxa
wide_asv <- wide_asv[, colSums(wide_asv != 0, na.rm = T) > 0]
dim(wide_asv)

write_csv(wide_asv, '../output/gab_asv_table.csv')
