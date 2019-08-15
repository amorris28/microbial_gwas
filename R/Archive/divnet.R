library(phyloseq)
library(DivNet)

physeq <- readRDS('../output/physeq_raw.rds')

# Identify references otu for divnet
a <- as(otu_table(physeq), 'matrix')
a[a > 1] <- 1
ref_otu <- colnames(a)[colSums(a) == max(colSums(a))]

asv_div <- divnet(physeq, ncores = 28, base = ref_otu)

saveRDS(asv_div, '../output/divnet_asv.rds')
