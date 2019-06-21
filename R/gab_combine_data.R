library(tidyverse)
library(morris)

# Methanogen metadata
gen_attr_table <- read.csv('../output/gab_gen_attr_table.csv', stringsAsFactors = FALSE)
gen_attr_table$Site <- paste0(gen_attr_table$Site, 
				reverse_substr(gen_attr_table$Sample, 3, 3))

# Methanotrophy metadata
troph_attr_table <- read.csv('../output/gab_troph_attr_table.csv', stringsAsFactors=F)
troph_attr_table$Site <- paste0(troph_attr_table$Site, 
				reverse_substr(troph_attr_table$Sample, 3, 3))
## Rarefied ASV Table
#rare_asv_table <- read.csv('../output/gab_rare_asv_table.csv')

## DESeq2 variance stabilized asv table
vst_asv_table <- read.csv('../output/vst_asv_table.csv', stringsAsFactors = FALSE)

## Raw sequence counts
pa_asv_table <- vst_asv_table #read_csv('../output/gab_asv_table.csv')
asvs <- pa_asv_table[, -1] 
## Presence Absence
asvs[asvs>0] <- 1
pa_asv_table[, -1] <- asvs

#left_join(troph_attr_table, troph_asv_table, by = "Sample") %>%
#  write_csv('output/troph_total.csv')
geodist <- read.table('../output/geodist.tsv', stringsAsFactors = FALSE)

# Export
troph_total_vst <- left_join(troph_attr_table, geodist, by = "Site")
troph_total_vst <- inner_join(troph_total_vst, vst_asv_table, by = "Sample")
#troph_total_rare <- troph_total_rare[troph_total_rare$Wetland == 'Upland', ]
write_csv(troph_total_vst, '../output/gab_all_vst_troph.csv')

gen_total_vst <- left_join(gen_attr_table, geodist, by = "Site")
gen_total_vst <- inner_join(gen_total_vst, vst_asv_table, by = "Sample")
write_csv(gen_total_vst, '../output/gab_all_vst_gen.csv')

troph_total_pa <- left_join(troph_attr_table, geodist, by = "Site")
troph_total_pa <- inner_join(troph_total_pa, pa_asv_table, by = "Sample")
write_csv(troph_total_pa, '../output/gab_all_pa_troph.csv')

gen_total_pa <- left_join(gen_attr_table, geodist, by = "Site")
gen_total_pa <- inner_join(gen_total_pa, pa_asv_table, by = "Sample")
write_csv(gen_total_pa, '../output/gab_all_pa_gen.csv')

