library(tidyverse)
library(amorris)

gen_attr_table <- read_csv('../output/gab_gen_attr_table.csv')

gen_attr_table$Site <- paste0(gen_attr_table$Site, 
				reverse_substr(gen_attr_table$Sample, 3, 3))
troph_attr_table <- read.csv('../output/gab_troph_attr_table.csv', stringsAsFactors=F)
troph_attr_table$Site <- paste0(troph_attr_table$Site, 
				reverse_substr(troph_attr_table$Sample, 3, 3))
## Rarefied ASV Table
rare_asv_table <- read.csv('../output/gab_rare_asv_table.csv')

## DESeq2 variance stabilized asv table
# troph_asv_table <- read_csv('../output/vst_asv_table.csv')

## Raw sequence counts
pa_asv_table <- read_csv('../output/gab_asv_table.csv')
asvs <- pa_asv_table[, -1] 
## Presence Absence
asvs[asvs>0] <- 1
pa_asv_table[, -1] <- asvs

#left_join(troph_attr_table, troph_asv_table, by = "Sample") %>%
#  write_csv('output/troph_total.csv')
geodist <- read.table('../output/geodist.tsv')

# Export
troph_total_rare <- left_join(troph_attr_table, geodist, by = "Site")
troph_total_rare <- inner_join(troph_total_rare, rare_asv_table, by = "Sample")
troph_total_rare <- troph_total_rare[troph_total_rare$Wetland == 'Upland', ]
write_csv(troph_total_rare, '../output/gab_all_rare_troph.csv')

gen_total_rare <- left_join(gen_attr_table, geodist, by = "Site")
gen_total_rare <- inner_join(gen_total_rare, rare_asv_table, by = "Sample")
write_csv(gen_total_rare, '../output/gab_all_rare_gen.csv')

troph_total_pa <- left_join(troph_attr_table, geodist, by = "Site")
troph_total_pa <- inner_join(troph_total_pa, pa_asv_table, by = "Sample")
write_csv(troph_total_pa, '../output/gab_all_pa_troph.csv')

gen_total_pa <- left_join(gen_attr_table, geodist, by = "Site")
gen_total_pa <- inner_join(gen_total_pa, pa_asv_table, by = "Sample")
write_csv(gen_total_pa, '../output/gab_all_pa_gen.csv')

