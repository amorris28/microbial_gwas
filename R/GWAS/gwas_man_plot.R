library(tidyverse)
source('functions.R')
taxon_table <- read_csv('output/taxon_table.csv')
lowk_by_asv <- readRDS('output/lowk_fits.rds')
lowk_by_asv_com <- readRDS('../output/lowk_fits_c.rds')
lowk_by_asv_env <- readRDS('output/lowk_fits_env.rds')
lowk_by_asv_raw <- readRDS('output/lowk_fits_raw.rds')
lowk_by_asv_spa <- readRDS('output/lowk_fits_spa.rds')

vmax_by_asv <- readRDS('output/vmax_fits.rds')
vmax_by_asv_com <- readRDS('output/vmax_fits_com.rds')
vmax_by_asv_env <- readRDS('output/vmax_fits_env.rds')
vmax_by_asv_raw <- readRDS('output/vmax_fits_raw.rds')
vmax_by_asv_spa <- readRDS('output/vmax_fits_spa.rds')

# Create man plots for adjusted and un-adjusted abundances
ggsave('figures/man_rare_all.pdf', plot = man_plot(by_asv, taxon_table), 
       width = 4*1.2, height = 3*1.2)
ggsave('figures/man_rare_env.pdf', plot = man_plot(by_asv_env, taxon_table), 
       width = 4*1.2, height = 3*1.2)
ggsave('figures/man_rare_com.pdf', plot = man_plot(by_asv_com, taxon_table), 
       width = 4*1.2, height = 3*1.2)
ggsave('figures/man_rare_spa.pdf', plot = man_plot(by_asv_spa, taxon_table), 
       width = 4*1.2, height = 3*1.2)
ggsave('figures/man_rare_raw.pdf', plot = man_plot(by_asv_raw, taxon_table), 
       width = 4*1.2, height = 3*1.2)

## Bonferroni correction (quite conservative) would suggest p-vale of 1 * 10 -4
## citation: Genome-Wide Association Studies, William S. Bush,  Jason H. Moore
## in PLOS 
## alpha:
#alpha = 0.05/nrow(by_asv)
#
##### Label significant taxa
## Identify ASVs with -log(p value) > 
#lowk_asvs <- by_asv %>% 
#  mutate(tidy = map(model, broom::tidy)) %>% 
#  unnest(tidy) %>% 
#  filter(term == 'abund') %>% 
#  filter(-log10(p.value) > -log10(alpha)) %>% 
#  pull(asv)
## Grab taxonomy of significant ASVs
#imp_tax <- taxon_table %>% 
#  filter(asv %in% lowk_asvs)
## Pull out the finest scale of taxonomy for each ASV
#imp_tax$taxon <- apply(imp_tax, 1, function(r) { r[ 
#                      if( any(is.na(r))){ 
#                          (which(is.na(r))[1]-1) 
#                       }else{
#                          (length(r)-0)}
#                                      ] }
#       )
## Pull out just the finest taxonomic scale and the ASV id
#imp_tax  <- imp_tax %>% 
#	select(asv, taxon)
#taxon_table  <- taxon_table %>% 
#	left_join(imp_tax, by = 'asv')
