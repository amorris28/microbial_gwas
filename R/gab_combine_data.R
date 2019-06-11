#### Input: troph_attr_table and rare_asv_table.csv
#### Output: troph_total.csv

library(tidyverse)
library(amorris)

#gen_attr_table <- read_csv('output/gen_attr_table.csv')

troph_attr_table <- read_csv('output/troph_attr_table.csv')
troph_attr_table$Site <- paste0(troph_attr_table$Site, 
				reverse_substr(troph_attr_table$Sample, 3, 3))
## Rarefied ASV Table
troph_asv_table <- read_csv('output/rare_asv_table.csv')

## DESeq2 variance stabilized asv table
# troph_asv_table <- read_csv('output/vst_asv_table.csv')

## Raw sequence counts
#troph_asv_table <- read_csv('output/asv_table.csv')
#asvs <- troph_asv_table[, -1] 
### Presence Absence
#asvs[asvs>0] <- 1
#troph_asv_table[, -1] <- asvs

#left_join(troph_attr_table, troph_asv_table, by = "Sample") %>%
#  write_csv('output/troph_total.csv')
geodist <- read.table('output/geodist.tsv')
troph_attr_table <- left_join(troph_attr_table, geodist, by = "Site")
troph_attr_table <- left_join(troph_attr_table, troph_asv_table, by = "Sample")
  write_csv(troph_attr_table, '../output/gab_troph_total.csv')

# gen_attr_table %>% 
# left_join(asv_table, by = "Sample") %>%
# write_csv('output/gen_total.csv')
