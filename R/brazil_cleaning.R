library(tidyverse)

## Cleaning up metadata
metadata <- read_tsv('data/brazil/Dimensions_meta_data.txt')
metadata <- as.data.frame(metadata)
## Converting below detection limit to the detection limt
## Controversial decision - consider other options (zeros?)
metadata <- sapply(metadata,function(x) 
               gsub("<","",x) )
metadata <- as.data.frame(metadata)
metadata <- sapply(metadata,function(x) gsub(",",".",x))
metadata <- as_tibble(metadata)
metadata <- mutate_at(metadata, vars(pH:Error_N2O), as.numeric)
as.data.frame(metadata)

## ASVs
asvs <- read_tsv('data/brazil/16S_dada2_table_RDP_tax.txt')
colnames(asvs)
taxa_abund <- 
  asvs %>% 
  mutate(OTU_ID = paste0("asv_", OTU_ID)) %>% 
  select(-ConsensusLineage) %>% 
  gather(Sample_ID, Abundance, 2:75) %>% 
  spread(OTU_ID, Abundance)
gps <- read.table('data/brazil/Dimensions_GPS_dec_deg.txt', header = FALSE)
colnames(gps) <- c('Sample_ID', 'Lat', 'Long')

all_data <-
  gps %>% 
  right_join(metadata, by = 'Sample_ID') %>% 
  right_join(taxa_abund, by = 'Sample_ID')

all_data$library_size <- rowSums(taxa_abund[, -1])
write_tsv(all_data, 'output/dim_cleaned_data.tsv')
