library(DESeq2)
library(tidyverse)

asv_data <- read.delim('../raw_data/16S_dada2_table_RDP_tax_edit.txt', header = TRUE)
asv_mat <- as.matrix(asv_data[,-c(1, ncol(asv_data))])
rownames(asv_mat) <- asv_data$OTU_ID
# colnames(asv_mat) <- colnames(asv_data[, -1])
col <- data.frame(site = colnames(asv_data[, -c(1, ncol(asv_data))]))

dds <- DESeqDataSetFromMatrix(countData = asv_mat, design = ~ 1, colData = col)

dds
dds <- estimateSizeFactors(dds, type = 'poscounts')
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)

varstab_commat <- assay(vsd)
vsd_data <- 
	cbind(tibble(asv = rownames(varstab_commat)), as_tibble(varstab_commat)) %>%
  mutate(asv = paste('asv', asv, sep = "_")) 


asv_table <- vsd_data %>%
  select(asv, starts_with('dna')) %>% 
  gather(Sample, abund, 2:length(.)) %>% 
  spread(asv, abund) %>% 
  mutate(Sample = substring(Sample, 5))
asv_table[asv_table < 0] <- 0
write_csv(asv_table, '../output/vst_asv_table.csv')
