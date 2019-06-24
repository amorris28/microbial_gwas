
library(tidyverse)
library(morris)

# Import ------------------------------------------
raw_troph_data <- read_tsv('../data/gabon/Gabon_methanotrophy_table_scaled.txt')
raw_gen_data <- read_tsv('../data/gabon/Gabon_methanogenesis_table_scaled.txt')

# Organize Methanotroph ------------------------------------------
troph_data <- raw_troph_data %>% 
  dplyr::select(-Low_r2, -High_r2, -High_final_k) %>% 
  dplyr::select(Vmax, Low_final_k, Bulk_dens:PC3_16S_rna) %>% 
  mutate(Rep = reverse_substr(Sample_names, 3, 3)) %>% 
  mutate_at(vars(Site, Description, Rep, Land_type, Sample_names, Experiment), funs(factor)) %>% 
  select(Sample = Sample_names, Site, Rep, Experiment, Description, Land_type, 
         Vmax:C_percent, pmoa_copy_num:PC3_16S_rna)

# Add Wetland variable
troph_data %<>% left_join(data_frame(Land_type = unique(troph_data$Land_type), 
           Wetland = factor(c('Wetland', 'Upland', 'Upland', 'Wetland', 'Upland',
                              'Wetland', 'Upland', 'Upland', 'Upland'))), by = 'Land_type') %>% 
  select(Sample:Experiment, Wetland, Description:PC3_16S_rna)

## Add WFPS variable
#Site_Pd <- data_frame(
#  Site = factor(unique(troph_data$Site)),
#  Pd = c(2.65, 2.65, 2.65, 2.65, 1.55, 2.65, 1.55, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65)
#)

#troph_data %<>% left_join(Site_Pd, by = 'Site') %>% 
#  mutate(WFPS = calc_wfps(Bulk_dens, Mois_cont, Pd)) %>% 
#  select(Sample:Low_final_k, WFPS, Bulk_dens:PC3_16S_rna, -Pd)

# Separate numeric data (remove factors)
troph_num_data <- troph_data %>% 
  select_if(function(col) !is.factor(col))

# Remove response variables
num_predictors <- troph_num_data %>% dplyr::select(-Vmax, -Low_final_k)

# Attribute table
troph_attr_table <- troph_data %>% 
  select(-Rep, -Experiment, -Description, -Land_type)

# Organizing attribute table for easy lmer
atr_tbl_long <- troph_attr_table %>% 
  gather('attribute', 'value', pmoa_copy_num:PC3_16S_rna)
rna_dna <- data_frame(
  attribute = unique(atr_tbl_long$attribute),
  molecule = factor(c('dna', 'rna', 'rna', 'dna', 
                      'dna', 'dna', 'dna', 'rna', 
                      'rna', 'rna', 'dna', 'dna', 
                      'dna', 'rna', 'rna', 'rna', 
                      'dna', 'dna', 'dna', 'rna',
                      'rna', 'rna')))
atr_tbl_long <- atr_tbl_long %>% 
  left_join(rna_dna, by = 'attribute')

new_attribute <- data_frame(
  attribute = unique(atr_tbl_long$attribute),
  new_attribute = factor(c('pmoa_num','pmoa_num', 'troph_rel_abun', 'troph_rel_abun', 
                           'troph_rich', 'troph_shan', 'troph_even', 'troph_rich', 
                           'troph_shan', 'troph_even', 'PC1', 'PC2', 
                           'PC3', 'PC1', 'PC2', 'PC3', 
                           'PC1_16S', 'PC2_16S', 'PC3_16S', 'PC1_16S', 
                           'PC2_16S', 'PC3_16S'))
)

attr_tbl_by_mlcl <- atr_tbl_long %>% 
  left_join(new_attribute, by = 'attribute') %>% 
  select(-attribute) %>% 
  rename(attribute = new_attribute) %>% 
  gather(process, rate, Vmax:Low_final_k) %>% 
  spread(attribute, value)

by_mlcl_proc <-
  attr_tbl_by_mlcl %>% 
  group_by(molecule, process) %>% 
  nest()




# Organize Methanogen ------------------------------------------
gen_data <-
  raw_gen_data %>% 
  dplyr::select(-H_CH4, -A_CH4, -A_CH4_percent) %>% 
  mutate(Rep = reverse_substr(Sample_names, 3, 3), 
         Experiment = reverse_substr(Sample_names, 1, 2)) %>% 
  select(Sample = Sample_names, Site, Rep, Experiment, Description, Land_type, 
         CH4:C_percent, Copies:PC3_16S_tr_cdna) %>% 
  mutate_at(vars(Sample, Site, Rep, Experiment, Description, Land_type), funs(factor))

# Add CO2/CH4 ratio
gen_data <- gen_data %>% 
  mutate(CO2_CH4 = CO2 / CH4) %>% 
  dplyr::select(Sample:CH4, -CO2, CO2_CH4, H_CH4_percent:PC3_16S_tr_cdna)

# Add WFPS variable
gen_data <- troph_data %>% 
  select(Bulk_dens, Site, Rep) %>% 
  right_join(gen_data) %>% 
  mutate(Site = factor(Site))

#Site_Pd <- data_frame(
#  Site = factor(unique(gen_data$Site)),
#  Pd = c(2.65, 2.65, 1.55, 1.55, 2.65, 2.65)
#)

gen_data <-
  gen_data %>% #left_join(Site_Pd, by = 'Site') %>%
#  mutate(WFPS = calc_wfps(Bulk_dens, Mois_cont, Pd)) %>%
  select(Site:H_CH4_percent, Bulk_dens, Mois_cont:PC3_16S_tr_cdna)

# Gen Numeric data
gen_num_data <- gen_data %>% 
  select_if(function(col) !is.factor(col))

num_predictors <- gen_num_data %>% dplyr::select(-CH4, -H_CH4_percent)


# Attribute table
gen_attr_table <- gen_data %>% 
  select(-Rep, -Experiment, -Description, -Land_type)

# Remove extraneous attributes
colnames(gen_attr_table)
gen_attr_table <- gen_attr_table %>% 
  select(Site:Copies)
colnames(troph_attr_table)
troph_attr_table <- troph_attr_table %>% 
  select(Sample:pmoa_transcript_num)

write_csv(gen_attr_table, '../output/gab_gen_attr_table.csv')

write_csv(troph_attr_table, '../output/gab_troph_attr_table.csv')
