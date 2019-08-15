#' ---
#' title: "Royal Society paper"
#' subtitle: "Linking microbial communities to ecosystem functions: what we can learn from genotype-phenotype mapping in organisms"
#' author: "Andrew Morris"
#' date: "`r Sys.Date()`"
#' output: 
#'   html_document:
#'     number_sections: TRUE
#' ---
#' 

#+ r global_options, include=FALSE
knitr::opts_chunk$set(fig.width=8, fig.height=6, fig.path='Figs/',
                      echo=FALSE, warning=FALSE, message=FALSE)

#+ r load_packages, include=FALSE

library(tidyverse)
library(broom)
library(broman)
library(varComp)
library(vegan)
library(qvalue)
library(kableExtra)
library(phyloseq)
library(gap)
library(breakaway)
library(ape)
library(ggpubr)

source('functions.R')

#+ r options
options(knitr.kable.NA = '')
theme_set(theme_bw())

#+ r import_data
physeq <- readRDS('../output/physeq_vst.rds')
physeq_raw <- readRDS('../output/physeq_raw.rds')

# Calculate richness from untransformed read counts using breakaway
richness <- breakaway(physeq_raw)
plot(richness)
sample_data(physeq)$rich_breakaway <- summary(richness)$estimate

# Add variables
com_mat <- otu_table(physeq)
bray <- vegdist(otu_table(physeq))
pcoa <- pcoa(bray)

# Screeplot for bray-curtis PCoA
ggplot(mapping = aes(seq_len(length(pcoa$values[, 1])), pcoa$values[, 1])) +
       geom_col() +
       labs(y = "Eigenvalues", x = "Principal Coordinate")

sample_data(physeq)$pc1 <- pcoa$vectors[, 1]
sample_data(physeq)$pc2 <- pcoa$vectors[, 2]
sample_data(physeq)$pc3 <- pcoa$vectors[, 3]

# Convert methane oxidation from negative to positive for easier interpetation
sample_data(physeq)$pos_lowk <- -sample_data(physeq)$Low_final_k

sample_data(physeq)$X_km <- sample_data(physeq)$X / 1000
sample_data(physeq)$Y_km <- sample_data(physeq)$Y / 1000

#' # Traditional correlations

(pmoa_plot <- sample_data(physeq) %>% 
  ggplot(aes(x = pmoa_copy_num, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = expression(italic('pmoA') ~ 'copy number'), y = expression('Methane oxidation rate (-k)')))

(rich_plot <- sample_data(physeq) %>% 
  ggplot(mapping = aes(x = rich_breakaway, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Richness', y = expression('Methane oxidation rate (-k)')))

ggarrange(pmoa_plot, rich_plot,
          labels = "AUTO",
          ncol = 2, nrow = 1)

ggsave(file = '../figures/trad_cors.pdf', width = 4*2, height = 4)

pmoa_fit <- lm(Low_final_k ~ pmoa_copy_num, data = as(sample_data(physeq), 'data.frame'))
summary(pmoa_fit)
glance(pmoa_fit)

rich_fit <- lm(Low_final_k ~ rich_breakaway, data = as(sample_data(physeq), 'data.frame'))
summary(rich_fit)
glance(rich_fit)

rbind(tidy(pmoa_fit)[2, ], tidy(rich_fit)[2, ])

# Pull out important axes an plot pairwise with axes proportional to var expl

ggplot(sample_data(physeq), aes(X_km, Y_km)) +
  geom_point() +
  stat_smooth() +
  coord_fixed()


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

(G1 <- ggplot(sample_data(physeq), aes(x = pc1, y = pc2, color = geocode)) +
  geom_point(size = 3) +
  labs(x = paste0('PC 1 (', myround(pcoa$values[1, 2]*100, 1),'%)'), 
       y = paste0('PC 2 (', myround(pcoa$values[2, 2]*100, 1),'%)')) +
  coord_fixed(pcoa$values[2, 2]/pcoa$values[1, 2]) +
  scale_color_manual(name = NULL, labels = c('Site A', 'Site B', 'Site C'),
                     values = cbPalette[6:8]))

(G2 <- ggplot(sample_data(physeq), aes(x = pc1, y = pc3, color = geocode)) +
  geom_point(size = 3) +
  labs(x = paste0('PC 1 (', myround(pcoa$values[1, 2]*100, 1),'%)'), 
       y = paste0('PC 3 (', myround(pcoa$values[3, 2]*100, 1),'%)')) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[1, 2]) +
  scale_color_manual(name = NULL, labels = c('Site A', 'Site B', 'Site C'),
                     values = cbPalette[6:8]))

(G3 <- ggplot(sample_data(physeq), aes(x = pc2, y = pc3, color = geocode)) +
  geom_point(size = 3) +
  labs(x = paste0('PC 2 (', myround(pcoa$values[2, 2]*100, 1),'%)'), 
       y = paste0('PC 3 (', myround(pcoa$values[3, 2]*100, 1),'%)')) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[2, 2]) +
  scale_color_manual(name = NULL, labels = c('Site A', 'Site B', 'Site C'),
                     values = cbPalette[6:8]))

(W1 <- ggplot(sample_data(physeq), aes(x = pc1, y = pc2, color = Wetland)) +
  geom_point(size = 3) +
  labs(x = paste0('PC 1 (', myround(pcoa$values[1, 2]*100, 1),'%)'), 
       y = paste0('PC 2 (', myround(pcoa$values[2, 2]*100, 1),'%)')) +
  coord_fixed(pcoa$values[2, 2]/pcoa$values[1, 2]) +
  scale_color_manual(name = NULL, values = cbPalette[2:3]))

(W2 <- ggplot(sample_data(physeq), aes(x = pc1, y = pc3, color = Wetland)) +
  geom_point(size = 3) +
  labs(x = paste0('PC 1 (', myround(pcoa$values[1, 2]*100, 1),'%)'), 
       y = paste0('PC 3 (', myround(pcoa$values[3, 2]*100, 1),'%)')) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[1, 2]) +
  scale_color_manual(name = NULL, values = cbPalette[2:3]))

(W3 <- ggplot(sample_data(physeq), aes(x = pc2, y = pc3, color = Wetland)) +
  geom_point(size = 3) +
  labs(x = paste0('PC 2 (', myround(pcoa$values[2, 2]*100, 1),'%)'), 
       y = paste0('PC 3 (', myround(pcoa$values[3, 2]*100, 1),'%)')) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[2, 2]) +
  scale_color_manual(name = NULL, values = cbPalette[2:3]))

ggarrange(ggarrange(G1, G2, G3,
          labels = c('A', 'B', 'C'),
          ncol = 1, nrow = 3,
          #align = 'v', 
          common.legend = TRUE),

ggarrange(W1, W2, W3,
          labels = c('D', 'E', 'F'),
          ncol = 1, nrow = 3,
          #align = 'v', 
          common.legend = TRUE))
ggsave('../figures/pcoa_multi.pdf', width = 5*2, height = 4*3)

#' ## Collinearity

#+ collinear_all_data

#+ r sim_mat

calc_sim_matrix <- function(x) {
1/(1+as.matrix(dist(as.matrix(scale(x)))))
}

## Community similarity
com_sim <- 1 - as.matrix(vegdist(otu_table(physeq), method = "bray"))
com_dis <- as.matrix(vegdist(otu_table(physeq), method = "bray"))
com_sim_bin <- 1 - as.matrix(vegdist(otu_table(physeq), method = "jaccard", binary = TRUE))
com_dis_bin <- as.matrix(vegdist(otu_table(physeq), method = "jaccard", binary = TRUE))

# Geographic similarity
geo_sim <- calc_sim_matrix(sample_data(physeq)[, c('X', 'Y')])
geo_dis <- as.matrix(dist(as.matrix(scale(sample_data(physeq)[, c('X', 'Y')]))))

# environmental similarity'
env_sim <- calc_sim_matrix(sample_data(physeq)[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])
env_dis <- as.matrix(dist(as.matrix(scale(sample_data(physeq)[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')]))))

ggplot(mapping = aes(x = c(geo_dis), y = c(com_dis))) + 
  geom_point()
ggplot(mapping = aes(x = c(geo_dis), y = c(env_dis))) + 
  geom_point()

print(mantel(com_sim_bin, geo_sim))
print(mantel(com_sim, geo_sim))
print(mantel(com_sim_bin, env_sim))
print(mantel(com_sim, env_sim))
print(mantel(geo_sim, env_sim))

#' # Variance Components

##+ r variance components
#Low_final_k <- sample_data(physeq)$Low_final_k
#geocode <- sample_data(physeq)$geocode
#com <- cholRoot(com_sim)
#geo <- geocode
##geo <- cholRoot(geo_sim)
#env <- cholRoot(env_sim)
#
#var_comp_pvalues <- data.frame(model = 
#                               c('CH4 \\textasciitilde 1 + com', 
#                                 'CH4 \\textasciitilde 1 + geo', 
#                                 'CH4 \\textasciitilde 1 + env', 
#                                 'CH4 \\textasciitilde 1 + geo + com', 
#                                 'CH4 \\textasciitilde 1 + geo + com', 
#                                 'CH4 \\textasciitilde 1 + geo + env', 
#                                 'CH4 \\textasciitilde 1 + geo + env', 
#                                 'CH4 \\textasciitilde 1 + com + env', 
#                                 'CH4 \\textasciitilde 1 + com + env', 
#                                 'CH4 \\textasciitilde 1 + geo + com + env', 
#                                 'CH4 \\textasciitilde 1 + geo + com + env', 
#                                 'CH4 \\textasciitilde 1 + geo + com + env'),
#           component = c('com', 'geo', 'env', 'geo', 'com', 'env', 'geo', 'com', 'env', 'env', 'geo', 'com'),
#           varcomp = NA,
#           se = NA,
#           p.value = NA)
#
## Single
## Com
#fit0 <- varComp(Low_final_k ~ 1)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[1] <- comps$h2G
#var_comp_pvalues$se[1] <- sqrt(comps$Varh2G)
#
#var_comp_pvalues$p.value[1] = p.value(varComp.test(fit, fit0))
#
## Geo
#fit <- varComp(Low_final_k ~ 1, random =  ~ geo)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[2] <- comps$h2G
#var_comp_pvalues$se[2] <- sqrt(comps$Varh2G)
#var_comp_pvalues$p.value[2] = p.value(varComp.test(fit, fit0))
#
## Env
#fit <- varComp(Low_final_k ~ 1, random =  ~ env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[3] <- comps$h2G
#var_comp_pvalues$se[3] <- sqrt(comps$Varh2G)
#var_comp_pvalues$p.value[3] = p.value(varComp.test(fit, fit0))
#
## Geo Com
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ com)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[4] <- comps$h2GE
#var_comp_pvalues$se[4] <- sqrt(comps$Varh2GE)
#var_comp_pvalues$p.value[4] = p.value(varComp.test(fit, fit0))
#
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ geo)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[5] <- comps$h2G
#var_comp_pvalues$se[5] <- sqrt(comps$Varh2G)
#var_comp_pvalues$p.value[5] = p.value(varComp.test(fit, fit0))
#
## Geo Env
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ geo)
#fit <- varComp(Low_final_k ~ 1, random =  ~ geo + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[6] <- comps$h2GE
#var_comp_pvalues$se[6] <- sqrt(comps$Varh2GE)
#var_comp_pvalues$p.value[6] = p.value(varComp.test(fit, fit0))
#
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ env)
#fit <- varComp(Low_final_k ~ 1, random =  ~ geo + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[7] <- comps$h2G
#var_comp_pvalues$se[7] <- sqrt(comps$Varh2G)
#var_comp_pvalues$p.value[7] = p.value(varComp.test(fit, fit0))
#
## Com Env
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ env)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[8] <- comps$h2G
#var_comp_pvalues$se[8] <- sqrt(comps$Varh2G)
#var_comp_pvalues$p.value[8] = p.value(varComp.test(fit, fit0))
#
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ com)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[9] <- comps$h2GE
#var_comp_pvalues$se[9] <- sqrt(comps$Varh2GE)
#var_comp_pvalues$p.value[9] = p.value(varComp.test(fit, fit0))
#
## Geo Com Env
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ com + geo)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[10] <- comps$h2GE[2]
#var_comp_pvalues$se[10] <- sqrt(comps$Varh2GE[2])
#var_comp_pvalues$p.value[10] = p.value(varComp.test(fit, fit0))
#
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ com + env)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[11] <- comps$h2GE[1]
#var_comp_pvalues$se[11] <- sqrt(comps$Varh2GE[1])
#var_comp_pvalues$p.value[11] = p.value(varComp.test(fit, fit0))
#
#fit0 <- varComp(Low_final_k ~ 1, random =  ~ geo + env)
#fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo + env)
#comps <- varCompGE(fit)
#var_comp_pvalues$varcomp[12] <- comps$h2G
#var_comp_pvalues$se[12] <- sqrt(comps$Varh2G)
#var_comp_pvalues$p.value[12] = p.value(varComp.test(fit, fit0))
#
#var_comp_pvalues$varcomp <- paste0(myround(var_comp_pvalues$varcomp, 2), ' (',
#                                   myround(var_comp_pvalues$se, 2), ')')
#
#var_comp_pvalues$p.value <- 
#  c(formatC(var_comp_pvalues$p.value[c(1)], format = "e", digits = 2), 
#    myround(var_comp_pvalues$p.value[2:7], 3),
#    formatC(var_comp_pvalues$p.value[c(8)], format = "e", digits = 2),
#    myround(var_comp_pvalues$p.value[9:12], 3))
#var_comp_pvalues[c(1, 2, 3, 7, 8), 5] <- 
#  cell_spec(var_comp_pvalues[c(1, 2, 3, 7, 8), 5],  "latex", bold = T)
##var_comp_pvalues[12, 5] <- 
##  cell_spec(var_comp_pvalues[12,5],  "latex", italic = T)
#var_comp_pvalues <- var_comp_pvalues[, -4]
#
#var_comp_pvalues %>% 
#  kable('latex', booktabs = T, digits = 2, escape = F,
#        caption = paste0('Variance component estimates from intercept-only models predicting the rate of high-affinity methane oxidation from soil. p-values
#        determined by linear score test comparing nested models with and without
#        each variance component. Significant p-values are bolded (p \\textless{} 0.05). 
#        n = 44.'), 
#        col.names = c('Model', 'Component', 'Estimate (SE)', 'p-value')) %>% 
#  collapse_rows(columns = 1) %>% 
#  writeLines('../tables/var_comp_pvalues.tex')

#' # Analyze model output

#+ r analyze_model_output
all_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_abund.csv'))
head(all_model_output)

#' ## Variance components

#+ r variance_components
all_model_output %>% 
  group_by(comps) %>% 
  summarize_at(vars(geo, com, env, err), mean) #%>% 
#  kable(caption = "Variance components for all data",
#        col.names = c('Components', 'Geography', 'Community', 'Environment', 'Error'),
#        digits = 2) %>% 
#  kable_styling()

# Calculate q values for FDR

for (i in unique(all_model_output$comps)) {
  print(i)
my_qs <- qvalue(all_model_output$p.value[all_model_output$comps == i])
summary(my_qs)
all_model_output$qvalues[all_model_output$comps == i] <- my_qs$qvalues
}

ggplot(all_model_output, aes(x = p.value)) +
  facet_wrap(~ comps) +
  geom_histogram()

#' ## Number of taxa identified with each covariate
#+ r n_taxa_ided
taxon_table <- read.csv('../output/taxon_table.csv')
n_taxa <- 
all_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05)  %>% 
  summarize(count = length(unique(asv)))

#' ## Taxa identified with each covariate

#' ### All Data Model

#+ r all_taxa


sig_taxa_all <- 
all_model_output  %>% 
left_join(taxon_table) %>% 
  group_by(comps) %>%
  filter(qvalues < 0.05 ) %>% 
  ungroup() %>% 
  select(comps, asv, Phylum:Genus)

taxa <- 
sig_taxa_all %>% 
  select(comps, Phylum:Genus)

sig_taxa_all %>% 
  filter(comps == 'gce')  %>% 
  arrange(Phylum:Genus) %>% 
  select(Phylum:Genus) %>% 
  kable('latex', booktabs = T, caption = 'Taxa significantly correlated with
        high-affinity methane oxidation rate after controlling for geographic
        proximity, environmental similarity, and community structure. 
        Significance determined by controlling
        the false discovery rate (q \\textless{} 0.05).') %>% 
  writeLines('../tables/all_taxa_sig.tex')
sig_taxa_gce <- sig_taxa_all %>% 
  filter(comps == 'gce')

sig_asv_ids <- substr(sig_taxa_gce$asv, 5, nchar(sig_taxa_gce$asv[1]))
asvs <- otu_table(physeq)
sig_asvs <- asvs[, sig_asv_ids]
sig_asvs <- as.data.frame(as(sig_asvs, 'matrix'))
sample_data <- as(sample_data(physeq), 'data.frame')
long_asvs <- 
cbind(sig_asvs, sample_data) %>% 
  gather(asv, abund, 1:6)

ggplot(long_asvs, aes(abund, pos_lowk)) +
  facet_wrap(~ asv) + 
  geom_point()

ggplot(long_asvs, aes(abund)) +
  facet_wrap(~ asv) + 
  geom_histogram()

effect_sizes <- data.frame(asv = colnames(asvs),
           estimate = NA,
           se = NA)

com_mat <- otu_table(physeq)[, !colnames(otu_table(physeq)) %in% sig_asv_ids[1]]
com_sim <- 1 - as.matrix(vegdist(com_mat, method = "bray"))

model <- varComp(pos_lowk ~ sig_asvs[, 1],
                 data = sample_data,
                 ~ geocode + cholRoot(com_sim) + cholRoot(env_sim))
summary(model)
effect_sizes[1, 2] <- 0.5395
effect_sizes[1, 3] <- 0.0823

com_mat <- otu_table(physeq)[, !colnames(otu_table(physeq)) %in% sig_asv_ids[2]]
com_sim <- 1 - as.matrix(vegdist(com_mat, method = "bray"))

model <- varComp(pos_lowk ~ sig_asvs[, 2],
                 data = sample_data,
                 ~ geocode + cholRoot(com_sim) + cholRoot(env_sim))
summary(model)
effect_sizes[2, 2] <- 0.5698
effect_sizes[2, 3] <- 0.0828

com_mat <- otu_table(physeq)[, !colnames(otu_table(physeq)) %in% sig_asv_ids[3]]
com_sim <- 1 - as.matrix(vegdist(com_mat, method = "bray"))

model <- varComp(pos_lowk ~ sig_asvs[, 3],
                 data = sample_data,
                 ~ geocode + cholRoot(com_sim) + cholRoot(env_sim))
summary(model)
effect_sizes[3, 2] <- 0.6171
effect_sizes[3, 3] <- 0.0930

com_mat <- otu_table(physeq)[, !colnames(otu_table(physeq)) %in% sig_asv_ids[4]]
com_sim <- 1 - as.matrix(vegdist(com_mat, method = "bray"))

model <- varComp(pos_lowk ~ sig_asvs[, 4],
                 data = sample_data,
                 ~ geocode + cholRoot(com_sim) + cholRoot(env_sim))
summary(model)
effect_sizes[4, 2] <- 0.5700
effect_sizes[4, 3] <- 0.0828

com_mat <- otu_table(physeq)[, !colnames(otu_table(physeq)) %in% sig_asv_ids[5]]
com_sim <- 1 - as.matrix(vegdist(com_mat, method = "bray"))

model <- varComp(pos_lowk ~ sig_asvs[, 5],
                 data = sample_data,
                 ~ geocode + cholRoot(com_sim) + cholRoot(env_sim))
summary(model)
effect_sizes[5, 2] <- 0.7911
effect_sizes[5, 3] <- 0.1261

com_mat <- otu_table(physeq)[, !colnames(otu_table(physeq)) %in% sig_asv_ids[6]]
com_sim <- 1 - as.matrix(vegdist(com_mat, method = "bray"))

model <- varComp(pos_lowk ~ sig_asvs[, 6],
                 data = sample_data,
                 ~ geocode + cholRoot(com_sim) + cholRoot(env_sim))
summary(model)
effect_sizes[6, 2] <- 1.5986
effect_sizes[6, 3] <- 0.2920

effect_sizes$asv <- paste("asv", effect_sizes$asv, sep = "_")
effect_sizes <- 
effect_sizes %>% 
  left_join(taxon_table)
effect_sizes$Genus <- as.character(effect_sizes$Genus)
effect_sizes[3, 'Genus'] <- 'Acidobacteria Gp1'  
taxon_labels <- effect_sizes$Genus
mylabels <- c(bquote(italic(.(taxon_labels[1]))),
                bquote(italic(.(taxon_labels[2]))),
                taxon_labels[3],
                bquote(italic(.(taxon_labels[4]))),
                bquote(italic(.(taxon_labels[5]))),
                bquote(italic(.(taxon_labels[6])))) 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(effect_sizes, aes(y = asv, x = estimate, xmin = estimate - se,
                         xmax = estimate + se, color = Phylum)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 1/4) +
  geom_vline(xintercept = 0) +
  labs(x = 'Estimate', y = 'ASV') +
  scale_y_discrete(labels = mylabels)+
  scale_colour_manual(values=cbPalette) +
  theme(axis.title.y = element_blank())
ggsave(file = '../figures/effect_sizes.pdf', width = 8, height = 4)

effect_sizes
# Plot largest estimate asv
asv <- as(otu_table(physeq)[, 'eb30799fb82dcac24f5d8629d8903ef2'], 'vector')
class(asv)
hist(asv)

sample_data <- as(sample_data(physeq), 'data.frame')
class(lowk)
ggplot(sample_data, mapping = aes(asv, pos_lowk, color = Land_type)) + 
  geom_point()
hist(otu_table(physeq)[, 'eb30799fb82dcac24f5d8629d8903ef2'])


sum(pull(sig_taxa_all[sig_taxa_all$comps == 'e', ], asv) %in% pull(sig_taxa_all[sig_taxa_all$comps == 'n', ], asv))

n_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'n', ], asv)
g_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'g', ], asv)
c_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'c', ], asv)
e_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'e', ], asv)
gc_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'gc', ], asv)
ge_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'ge', ], asv)
ce_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'ce', ], asv)
gce_ <- pull(sig_taxa_all[sig_taxa_all$comps == 'gce', ], asv)

comp_adds_removes <- 
data.frame(
           comps = c('none', 'geo', 'com', 'env', 'geo + com', 'geo + env', 'com + env', 'geo + com + env'), 
           removed = c(0, sum(!n_ %in% g_), sum(!n_ %in% c_), sum(!n_ %in% e_), 
                       sum(!n_ %in% gc_), sum(!n_ %in% ge_), sum(!n_ %in% ce_), 
                       sum(!n_ %in% gce_)),
           added = c(0, sum(!g_ %in% n_), sum(!c_ %in% n_), sum(!e_ %in% n_), 
                       sum(!gc_ %in% n_), sum(!ge_ %in% n_), sum(!ce_ %in% n_), 
                       sum(!gce_ %in% n_)),
           n_sig = c(length(n_), length(g_), length(c_), length(e_), length(gc_),
                     length(ge_), length(ce_), length(gce_)))

comp_adds_removes %>% 
  kable('latex', booktabs = TRUE, 
        caption = paste0("Number of taxa correlated with high-affinity
                         methane oxidation. Tested ", 
                         length(unique(all_model_output$asv)), " taxa.
        Significant taxa determined by controlling the false-discovery rate
        (q \\textless{} 0.05). n = ", nrow(all_data)),
        col.names = c('Components', 'Removed', 'Added', 'Significant')) %>% 
  kable_styling() %>% 
  writeLines('../tables/all_taxa_n.tex')

# Shared

sum(!g_ %in% ce_)
sum(!ce_ %in% g_)
sum(!c_ %in% e_)
sum(!c_ %in% e_)

comp_shared <- 
data.frame(
           comps = c('c in g', 'c in e', 'c in gc', 'c in ge', 'c in ce', 'c in gce', 
                     'g in c', 'g in e', 'g in gc', 'g in ge', 'g in ce', 'g in gce',
                     'e in c', 'e in c', 'e in gc', 'e in ge', 'e in ce', 'e in gce'),
           comparisons = c(sum(!c_ %in% g_), sum(!c_ %in% e_), sum(!c_ %in% gc_), 
                       sum(!c_ %in% ge_), sum(!c_ %in% ce_), sum(!c_ %in% gce_), 
                       sum(!g_ %in% e_), sum(!g_ %in% ge_), sum(!g_ %in% gc_)),
           added = c(0, sum(!g_ %in% n_), sum(!c_ %in% n_), sum(!e_ %in% n_), 
                       sum(!gc_ %in% n_), sum(!ge_ %in% n_), sum(!ce_ %in% n_), 
                       sum(!gce_ %in% n_)),
           n_sig = c(length(n_), length(g_), length(c_), length(e_), length(ge_),
                     length(gc_), length(ce_), length(gce_)))

#' ### Wetland Model

#+ r wetland_taxa

wet_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / (length(unique(asv)))) %>% 
  select(comps, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  kable() %>% 
  kable_styling() #%>% 
#  collapse_rows(columns = 1)

sig_taxa_wet <- 
wet_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05 ) %>% 
  select(comps, asv, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus)

sum(pull(sig_taxa_wet[sig_taxa_wet$comps == 'e', ], asv) %in% pull(sig_taxa_wet[sig_taxa_wet$comps == 'n', ], asv))

n_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'n', ], asv)
g_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'g', ], asv)
c_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'c', ], asv)
e_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'e', ], asv)
gc_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'gc', ], asv)
ge_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'ge', ], asv)
ce_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'ce', ], asv)
gce_ <- pull(sig_taxa_wet[sig_taxa_wet$comps == 'gce', ], asv)

comp_adds_removes <- 
data.frame(
           comps = c('n', 'g', 'c', 'e', 'gc', 'ge', 'ce', 'gce'), 
           removed = c(0, sum(!n_ %in% g_), sum(!n_ %in% c_), sum(!n_ %in% e_), 
                       sum(!n_ %in% gc_), sum(!n_ %in% ge_), sum(!n_ %in% ce_), 
                       sum(!n_ %in% gce_)),
           added = c(0, sum(!g_ %in% n_), sum(!c_ %in% n_), sum(!e_ %in% n_), 
                       sum(!gc_ %in% n_), sum(!ge_ %in% n_), sum(!ce_ %in% n_), 
                       sum(!gce_ %in% n_)),
           n_sig = c(length(n_), length(g_), length(c_), length(e_), length(ge_),
                     length(gc_), length(ce_), length(gce_)))

#' ### Upland Model

#+ r upland_taxa

dry_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05) %>% 
  select(comps, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  kable() %>% 
  kable_styling() %>% 
#  collapse_rows(columns = 1) %>% 
  save_kable(file = '../tables/dry_taxa_qvalue.html')

sig_taxa_dry_geo <- 
dry_model_output_geo  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05 ) %>% 
  select(comps, asv, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus)

sig_taxa_dry <- 
dry_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05 ) %>% 
  select(comps, asv, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus)


n_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'n', ], asv)
g_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'g', ], asv)
c_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'c', ], asv)
e_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'e', ], asv)
gc_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'gc', ], asv)
ge_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'ge', ], asv)
ce_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'ce', ], asv)
gce_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'gce', ], asv)

comp_adds_removes <- 
data.frame(
           comps = c('n', 'g', 'c', 'e', 'gc', 'ge', 'ce', 'gce'), 
           removed = c(0, sum(!n_ %in% g_), sum(!n_ %in% c_), sum(!n_ %in% e_), 
                       sum(!n_ %in% gc_), sum(!n_ %in% ge_), sum(!n_ %in% ce_), 
                       sum(!n_ %in% gce_)),
           added = c(0, sum(!g_ %in% n_), sum(!c_ %in% n_), sum(!e_ %in% n_), 
                       sum(!gc_ %in% n_), sum(!ge_ %in% n_), sum(!ce_ %in% n_), 
                       sum(!gce_ %in% n_)),
           n_sig = c(length(n_), length(g_), length(c_), length(e_), length(gc_),
                     length(ge_), length(ce_), length(gce_)))
#comp_adds_removes$math <- 124 - comp_adds_removes$removed + comp_adds_removes$added
comp_adds_removes %>% 
  kable('latex', booktabs = T, caption = "N taxa identified with upland data",
        col.names = c('Components', 'Removed', 'Added', 'Significant')) %>% 
  kable_styling() %>% 
  writeLines('../tables/dry_taxa_n.tex')
#' ### Binary Upland Model

#+ r upland_taxa

dry_model_output_bin  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05) %>% 
  select(comps, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  kable() %>% 
  kable_styling() %>% 
#  collapse_rows(columns = 1) %>% 
  save_kable(file = '../tables/dry_taxa_qvalue_bin.html')

sig_taxa_dry <- 
dry_model_output_bin  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05 ) %>% 
  select(comps, asv, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus)


n_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'n', ], asv)
g_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'g', ], asv)
c_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'c', ], asv)
e_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'e', ], asv)
gc_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'gc', ], asv)
ge_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'ge', ], asv)
ce_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'ce', ], asv)
gce_ <- pull(sig_taxa_dry[sig_taxa_dry$comps == 'gce', ], asv)

comp_adds_removes <- 
data.frame(
           comps = c('n', 'g', 'c', 'e', 'gc', 'ge', 'ce', 'gce'), 
           removed = c(0, sum(!n_ %in% g_), sum(!n_ %in% c_), sum(!n_ %in% e_), 
                       sum(!n_ %in% gc_), sum(!n_ %in% ge_), sum(!n_ %in% ce_), 
                       sum(!n_ %in% gce_)),
           added = c(0, sum(!g_ %in% n_), sum(!c_ %in% n_), sum(!e_ %in% n_), 
                       sum(!gc_ %in% n_), sum(!ge_ %in% n_), sum(!ce_ %in% n_), 
                       sum(!gce_ %in% n_)),
           n_sig = c(length(n_), length(g_), length(c_), length(e_), length(gc_),
                     length(ge_), length(ce_), length(gce_)))
#comp_adds_removes$math <- 124 - comp_adds_removes$removed + comp_adds_removes$added
comp_adds_removes %>% 
  kable(caption = "N taxa identified with binary upland data",
        col.names = c('Components', 'Removed', 'Added', 'n Significant')) %>% 
  kable_styling() %>% 
  save_kable(file = '../tables/dry_taxa_n_qvalue_bin.html')

#cbind(added = rep('a', 50), sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'e', ], asv), ]),
#sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'c', ], asv), ],
#sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'g', ], asv), ],
#sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'gc', ], asv), ],
#sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'ge', ], asv), ],
#sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'ce', ], asv), ],
#sig_taxa[sig_taxa$comps == 'n', ][!pull(sig_taxa[sig_taxa$comps == 'n', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'gce', ], asv), ])
#
#sig_taxa[sig_taxa$comps == 'e', ][!pull(sig_taxa[sig_taxa$comps == 'e', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ],
#sig_taxa[sig_taxa$comps == 'c', ][!pull(sig_taxa[sig_taxa$comps == 'c', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ],
#sig_taxa[sig_taxa$comps == 'g', ][!pull(sig_taxa[sig_taxa$comps == 'g', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ],
#
#sig_taxa[sig_taxa$comps == 'gc', ][!pull(sig_taxa[sig_taxa$comps == 'gc', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ],
#
#sig_taxa[sig_taxa$comps == 'ge', ][!pull(sig_taxa[sig_taxa$comps == 'ge', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ],
#
#sig_taxa[sig_taxa$comps == 'ce', ][!pull(sig_taxa[sig_taxa$comps == 'ce', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ],
#
#sig_taxa[sig_taxa$comps == 'gce', ][!pull(sig_taxa[sig_taxa$comps == 'gce', ], asv) %in% pull(sig_taxa[sig_taxa$comps == 'n', ], asv), ])
#)
#
#
#



