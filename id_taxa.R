library(tidyverse)
library(xtable)
# import model output
by_asv_all <- readRDS('output/model_fits.rds')
by_asv_raw <- readRDS('output/model_fits_raw.rds')
by_asv_com <- readRDS('output/model_fits_com.rds')
by_asv_env <- readRDS('output/model_fits_env.rds')
by_asv_spa <- readRDS('output/model_fits_spa.rds')

new_data <- read_csv('output/adj_data.csv')
taxon_table <- read_csv('output/taxon_table.csv')

alpha <- 0.05 / nrow(by_asv_all)
# Identify ASVs with -log(p value) > 
id_asvs <- function(by_asv, alpha, taxon_table, new_data, file_name) {
lowk_asvs <- by_asv %>% 
  mutate(tidy = map(model, broom::tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'abund') %>% 
  filter(-log10(p.value) > -log10(alpha)) %>% 
  pull(asv)

# Pull out those ASVs
important_asvs <- new_data %>% 
  select(Low_final_k, one_of(lowk_asvs))%>% 
  gather(asv, abund, starts_with('asv'))
taxon_table %>% 
  filter(asv %in% lowk_asvs) %>% 
  select(-asv) %>%
  xtable() %>% 
  print.xtable(file = file_name)
}

id_asvs(by_asv_raw, alpha, taxon_table, new_data, file_name = "tables/table_rare_raw.tex")
id_asvs(by_asv_env, alpha, taxon_table, new_data, file_name = "tables/table_rare_env.tex")
id_asvs(by_asv_spa, alpha, taxon_table, new_data, file_name = "tables/table_rare_spa.tex")
id_asvs(by_asv_com, alpha, taxon_table, new_data, file_name = "tables/table_rare_com.tex")
id_asvs(by_asv_all, alpha, taxon_table, new_data, file_name = "tables/table_rare_all.tex")


## Individual correlation
imp_tax <- taxon_table %>% 
  filter(asv %in% important_asvs$asv)
imp_tax$taxon <- apply(imp_tax, 1, function(r) { r[ 
                      if( any(is.na(r))){ 
                          (which(is.na(r))[1]-1) 
                       }else{
                          (length(r)-0)}
                                      ] }
       )
imp_tax <- imp_tax %>% select(asv, taxon)
imp_tax$index <- seq(from=1, to=nrow(imp_tax))
#imp_tax$level <- as.character(1:5)
imp_data <-
	new_data %>% 
	select(Low_final_k, one_of(imp_asv_list$asv))
p <- imp_data %>% 
	gather(asv, abund, 2:ncol(imp_data)) %>% 
	left_join(imp_tax) %>% 
	ggplot(aes(x = abund, y = Low_final_k)) + 
	       geom_point() +
	       facet_wrap(. ~ index)+
	       theme_bw() +
	       labs(x = 'ASV Abundance', y = 'Methane Oxidation')
p
       ggsave('figures/sig_taxa.pdf', p, width = 8, height = 8)

summary(lm(Low_final_k ~ ., imp_data))
imp_data %>% 
	select(Low_final_k, starts_with('asv_2')) %>% 
	lm(Low_final_k ~ ., data = .)  %>% 
	summary()
imp_data %>% 
	select(Low_final_k, starts_with('asv_7')) %>% 
	lm(Low_final_k ~ ., data = .)  %>% 
	summary()
imp_data %>% 
	select(Low_final_k, starts_with('asv_1')) %>% 
	lm(Low_final_k ~ ., data = .)  %>% 
	summary()
imp_data %>% 
	select(Low_final_k, starts_with('asv_935')) %>% 
	gather(asv, abund, 2) %>% 
	ggplot(aes(x = abund, y = Low_final_k)) + 
	       geom_point() +
	       theme_bw() +
	       labs(x = expression('Abundance of' ~ italic(Nitrososphaera)), y = expression('Methane Oxidation Rate'~(italic(k)))) +
	       stat_smooth(method = 'lm') + 
	       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file = 'figures/eg_taxon.pdf', width = 4, height = 4)


fit <- lm(Low_final_k ~ ., data = imp_data)
summary(fit)
#car::avPlots(fit)
fit <- lm(Low_final_k ~ asv_9359e9c663671afaac28202af4c712bd, data  = imp_data)
summary(fit)
## Uncorrected values




by_asv_raw <- readRDS('output/model_fits_raw.rds')
raw_data <- read_csv('output/troph_uncorrected.csv')
alpha <- 0.05 / nrow(by_asv_raw)
# Identify ASVs with -log(p value) > 
lowk_asvs_raw <- by_asv_raw %>% 
  mutate(tidy = map(model, broom::tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'abund') %>% 
  filter(-log10(p.value) > -log10(alpha)) %>% 
  arrange(p.value) %>% 
  pull(asv)

# Pull out those ASVs
important_asvs_raw <- raw_data %>% 
  select(Low_final_k, one_of(lowk_asvs_raw))

imp_asv_list_raw <- important_asvs_raw %>% 
  gather(asv, abund, starts_with('asv'))
taxon_table %>% 
  filter(asv %in% imp_asv_list_raw$asv) %>% 
  select(-asv) %>%
  xtable() %>% 
  print.xtable(file = "tables/table_rare_non.tex")


fit <- lm(Low_final_k ~ asv_9359e9c663671afaac28202af4c712bd, data  = raw_data)
summary(fit)
raw_data %>% 
	select(Low_final_k, asv_9359e9c663671afaac28202af4c712bd) %>% 
	gather(asv, abund, 2) %>% 
	ggplot(aes(x = abund, y = Low_final_k)) + 
	       geom_point() +
	       theme_bw() +
	       labs(x = 'Abundance of a taxon (ASV)', y = expression('Methane Oxidation Rate'~(italic(k)))) +
	       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file = 'figures/eg_taxon_raw.pdf', width = 4, height = 4)

