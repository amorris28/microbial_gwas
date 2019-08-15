library(tidyverse)
library(morris)
library(phyloseq)

# Methanotrophy metadata
sample_data <- data.table::fread('../output/gab_troph_attr_table.csv', data.table = FALSE)
sample_data$Site <- paste0(sample_data$Site, 
				reverse_substr(sample_data$Sample, 3, 3))
geodist <- read.table('../output/geodist.tsv', stringsAsFactors = FALSE)
sample_data$Experiment <- reverse_substr(sample_data$Sample, 1, 2)
sample_data <- left_join(sample_data, geodist, by = "Site")
rownames(sample_data) <- sample_data$Sample
sample_data <- sample_data[, -1]

sample_data[sample_data$Y > 30000, 'geocode'] <- 'A'
sample_data[sample_data$Y < 30000 & sample_data$Y > 10000, 'geocode'] <- 'B'
sample_data[sample_data$Y < 10000, 'geocode'] <- 'C'
sample_data$geocode <- factor(sample_data$geocode)

## DESeq2 variance stabilized asv table
otu_table <- data.table::fread('../output/vst_asv_table.csv', data.table = FALSE)
rownames(otu_table) <- otu_table[, 1]
otu_table <- otu_table[, -1]
colnames(otu_table) <- substr(colnames(otu_table), 5, nchar(colnames(otu_table)[1]))
otu_table <- as.matrix(otu_table)

## Raw asv table
raw_otu_table <- data.table::fread('../output/gab_asv_table.csv', data.table = FALSE)
rownames(raw_otu_table) <- raw_otu_table[, 1]
raw_otu_table <- raw_otu_table[, -1]
colnames(raw_otu_table) <- substr(colnames(raw_otu_table), 5, nchar(colnames(raw_otu_table)[1]))
raw_otu_table <- as.matrix(raw_otu_table)

## Taxonomy
taxon_table <- data.table::fread('../output/taxon_table.csv', data.table = FALSE)
taxon_table$Species <- ""
rownames(taxon_table) <- substr(taxon_table$asv, 5, nchar(taxon_table$asv[1]))
taxon_table <- taxon_table[, -1]
taxon_table <- as.matrix(taxon_table)

## Read tree
tree <- read_tree('../raw_data/16S_rep_seqs_tree_Gabon.nwk')

## Make phyloseq object
ASV <- otu_table(otu_table, taxa_are_rows = FALSE)
TAX <- tax_table(taxon_table)
sample_data <- sample_data(sample_data)

physeq <- phyloseq(ASV, TAX, sample_data, tree)
physeq <- 
physeq %>% 
  subset_samples(Experiment == 'MO')
#physeq <- corncob::clean_taxa_names(physeq)
otu_table(physeq) <- otu_table(physeq)[, colSums(otu_table(physeq)) > 0]
saveRDS(physeq, '../output/physeq_vst.rds')


ASV <- otu_table(raw_otu_table, taxa_are_rows = FALSE)
physeq <- phyloseq(ASV, TAX, sample_data, tree)
physeq <- 
physeq %>% 
  subset_samples(Experiment == 'MO')
#physeq <- corncob::clean_taxa_names(physeq)
otu_table(physeq) <- otu_table(physeq)[, colSums(otu_table(physeq)) > 0]
saveRDS(physeq, '../output/physeq_raw.rds')
