library(tidyverse)
library(DESeq2)
library(data.table)

## Cleaning up metadata
metadata <- read_tsv('../data/brazil/Dimensions_meta_data.txt')
metadata <- as.data.frame(metadata)
## Converting below detection limit to the detection limt
## Controversial decision - consider other options (zeros?)
metadata <- sapply(metadata,function(x) 
               gsub("<","",x) )
metadata <- as.data.frame(metadata)
metadata <- sapply(metadata,function(x) gsub(",",".",x))
metadata <- as_tibble(metadata)
metadata <- mutate_at(metadata, vars(pH:Error_N2O), as.numeric)
## ASVs
asvs <- as.data.frame(data.table::fread('../data/brazil/16S_dada2_table_RDP_tax.txt'))

# Pull out taxonomy
taxon_table <- 
  asvs %>% 
  select(OTU_ID, ConsensusLineage) %>% 
  mutate(OTU_ID = paste('asv', OTU_ID, sep = "_")) %>% 
  select(asv = OTU_ID, ConsensusLineage) %>% 
  separate(ConsensusLineage, c('Domain', 'Phylum', 'Class', 'Order', 'Family',
                               'Genus', 'Species', 'Subspecies'), ";")
fwrite(taxon_table, '../output/brazil_taxon_table.csv')

# Pull out abundances and do variance stabilizing transformation
asvs <-
  asvs %>% 
  select(-ConsensusLineage)
asv_mat <- asvs[,-1]
rownames(asv_mat) <- asvs$OTU_ID
col <- data.frame(site = colnames(asv_mat))
dds <- DESeqDataSetFromMatrix(countData = asv_mat, design = ~ 1, colData = col)
dds <- estimateSizeFactors(dds, type = 'poscounts')
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)
varstab_commat <- assay(vsd)
vsd_data <- 
	cbind(tibble(asv = rownames(varstab_commat)), as_tibble(varstab_commat)) %>%
  mutate(asv = paste('asv', asv, sep = "_")) 

taxa_abund <- 
  vsd_data %>% 
  gather(Sample_ID, Abundance, 2:75) %>% 
  spread(asv, Abundance)
taxa_abund[taxa_abund < 0] <- 0

gps <- read.table('../data/brazil/Dimensions_GPS_dec_deg.txt', header = FALSE)
colnames(gps) <- c('Sample_ID', 'Lat', 'Long')

all_data <-
  gps %>% 
  right_join(metadata, by = 'Sample_ID') %>% 
  right_join(taxa_abund, by = 'Sample_ID')

#all_data$library_size <- rowSums(taxa_abund[, -1])
write_tsv(all_data, '../output/brazil_cleaned_data.tsv')
