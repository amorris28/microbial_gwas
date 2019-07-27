library(tidyverse)
library(morris)
library(phyloseq)

# Methanotrophy metadata
troph_attr_table <- data.table::fread('../output/gab_troph_attr_table.csv', data.table = FALSE)
troph_attr_table$Site <- paste0(troph_attr_table$Site, 
				reverse_substr(troph_attr_table$Sample, 3, 3))
# Methanogenesis metadata
gen_attr_table <- data.table::fread('../output/gab_gen_attr_table.csv', data.table = FALSE)
gen_attr_table$Site <- paste0(gen_attr_table$Site, 
				reverse_substr(gen_attr_table$Sample, 3, 3))
geodist <- read.table('../output/geodist.tsv', stringsAsFactors = FALSE)

sample_data <- full_join(troph_attr_table, gen_attr_table)
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
class(otu_table)
colnames(otu_table)

## Taxonomy
taxon_table <- data.table::fread('../output/taxon_table.csv', data.table = FALSE)
taxon_table$Species <- ""
rownames(taxon_table) <- substr(taxon_table$asv, 5, nchar(taxon_table$asv[1]))
taxon_table[1:10, ]
taxon_table <- taxon_table[, -1]
taxon_table <- as.matrix(taxon_table)

## Read tree
tree <- read_tree('../data/gabon/16S_rep_seqs_tree_Gabon.nwk')
class(tree)

## Make phyloseq object
ASV <- otu_table(otu_table, taxa_are_rows = FALSE)
sample_names(ASV)
TAX <- tax_table(taxon_table)
sample_names(TAX)
sample_data <- sample_data(sample_data)
sample_names(sample_data)

physeq <- phyloseq(ASV, TAX, sample_data, tree)
saveRDS(physeq, '../output/physeq.rds')

