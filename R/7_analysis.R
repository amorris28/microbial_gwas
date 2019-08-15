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
#library(knitr)
library(kableExtra)
library(phyloseq)
library(gap)
#library(lmtest)
source('functions.R')
#library(DivNet)
library(breakaway)
library(ape)

#+ r options
options(knitr.kable.NA = '')
theme_set(theme_bw())

#+ r import_data
physeq <- readRDS('../output/physeq_vst.rds')
physeq_raw <- readRDS('../output/physeq_raw.rds')
#divnet <- readRDS('../talapas-output/divnet_asv.rds')
richness <- breakaway(physeq_raw)
plot(richness)

sample_data(physeq)$rich_breakaway <- summary(richness)$estimate
# Add variables
com_mat <- otu_table(physeq)
sample_data(physeq)$rich <- rowSums(com_mat > 0)
sample_data(physeq)$shan <- diversity(com_mat, "shannon")
sample_data(physeq)$simp <- diversity(com_mat, "simpson")
bray <- vegdist(otu_table(physeq))
pcoa <- pcoa(bray)

plot(sample_data(physeq)$rich, sample_data(physeq)$rich_breakaway)

ggplot(mapping = aes(seq_len(length(pcoa$values[, 1])), pcoa$values[, 1])) +
       geom_col() +
       labs(y = "Eigenvalues", x = "Principal Coordinate")

sample_data(physeq)$pc1 <- pcoa$vectors[, 1]
sample_data(physeq)$pc2 <- pcoa$vectors[, 2]
sample_data(physeq)$pc3 <- pcoa$vectors[, 3]

sample_data(physeq)$pos_lowk <- -sample_data(physeq)$Low_final_k

sample_data(physeq)$X_km <- sample_data(physeq)$X / 1000
sample_data(physeq)$Y_km <- sample_data(physeq)$Y / 1000
#' # Overview
#' 
#' The goal of this project is to answer the question: is microbial community 
#' composition important for determining the rate of ecosystem-scale functions independent 
#' of the underlying environmental variation? To address this, I partitioned
#' variation between community and environment (while controlling for geographic
#' location using site as a factor). I also identified taxa that are significantly
#' correlated with function after accounting for the covariance of geography,
#' environment, and community similarity. Currently, significant taxa are identified
#' using the Bonferroni correction (with the number of non-independent tests equal
#' to the number of taxa included in each model). Eventually, I would like to
#' determine p.values using a permutation test.
#'
#' These data were collected from Gabon, Africa and include `r n_samp` samples.
#' Microbial community data include 16S rRNA gene sequences (only DNA not RNA/cDNA)
#' with `r n_asvs` ASVs from the DADA2 pipeline. The community matrix was transformed
#' using the `varianceStabilizingTransformation` from `DESeq2`.
#' This analysis only includes high affinity
#' methane oxidation measurements (Low Final K). The spatial data were converted
#' from latitude/longitude to distance in meters by converting to UTM. 
#' Environmental data include bulk density, soil moisture, percent nitrogen, and 
#' percent carbon.
#'
#' # Overall Summary
#'
#' The traditional approach to microbial ecosystem function studies as I argue in the
#' paper is to correlate some aspect of microbial biodiversity such as abundance,
#' richness, diversity, or composition with the rate of an ecosystem function.
#' This largely fails as demonstrated by the Rocca and Graham papers. Here I
#' perform some of these analyses (PC regression, diversity, copy number/transcript
#' number, and various environmental factors) and my results agree with them.
#' Basically, the effectiveness of these hypotheses varies by function, but generally
#' don't hold up looking across ecosystems. For Low_final_k, a phylogenetically
#' conserved function (Martiny et al.) composition seems important. Although,
#' I'm arguing methanotrophs might not be the important limiting factor so
#' conservatism might not be important. Using PC regression, I identify some
#' candidate taxa.
#' 
#' In the paper I argue that this problem is less analagous to traditional
#' BEF literature where some variation in biodiversity is manipulated or 
#' measured in the field. Generally the function of interest there is biomass or
#' productivity at the community-level. In our case, we're sampling the
#' natural variation in microbial community composition across a landscape that
#' varies in environmental conditions and geographic distance between sites.
#' In addition, the shared history and community interactions of a site
#' could like to species associations unrelated to ecosystem function.
#' Therefore, we need to take into account these covariates, just like in 
#' GWAS studies.
#'
#' To that end, I use linear mixed-effects models to estimate the variance
#' explained by community and environment and used location as a factor
#' to control for variation among the clusters of sites. After including
#' these covariates in the model, many of the taxa drop out. 
#'
#' # Traditional correlations


(pmoa_plot <- sample_data(physeq) %>% 
  ggplot(aes(x = pmoa_copy_num, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = expression(italic('pmoA') ~ 'copy number'), y = expression('Methane oxidation rate (-k)')))
#  annotate(geom = 'text', x = min(sample_data(physeq)$pmoa_copy_num, na.rm = TRUE), 
#            y = max(sample_data(physeq)$pos_lowk), label = 'A', size = 6)

(rich_plot <- sample_data(physeq) %>% 
  ggplot(mapping = aes(x = rich_breakaway, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Richness', y = expression('Methane oxidation rate (-k)')))
#  annotate(geom = 'text', x = min(sample_data(physeq)$shan),
#           y = max(sample_data(physeq)$pos_lowk), label = 'B', size = 6)

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

#shan_fit_no_outlier <- lm(Low_final_k ~ shan, data = as(sample_data(physeq)[sample_data(physeq)$shan > 5, ], 'data.frame'))
#summary(shan_fit_no_outlier)
#glance(shan_fit_no_outlier)
#
rbind(tidy(pmoa_fit)[2, ], tidy(rich_fit)[2, ])

# Pull out important axes an plot pairwise with axes proportional to var expl

Z <- predict(loess(Y_km ~ X_km, as(sample_data(physeq), 'data.frame'))) 

ggplot(sample_data(physeq), aes(X_km, Y_km)) +
  geom_point() +
  stat_smooth() +
  coord_fixed()

ggplot(sample_data(physeq), aes(x = X_km, y = Y_km, color = Mois_cont)) +
  geom_point(size = 3, alpha = 0.1) +
  scale_color_gradient(name = "Moisture") +
  scale_shape(name = 'Site') +
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

#ggplot(sample_data(physeq), aes(x = pc1, y = pc2, color = interaction(geocode, Wetland))) +
#  geom_point(size = 3) +
#  labs(x = paste0('PC 1 (', myround(pcoa$values[1, 2]*100, 1),'%)'), 
#       y = paste0('PC 2 (', myround(pcoa$values[2, 2]*100, 1),'%)')) +
#  coord_fixed(pcoa$values[2, 2]/pcoa$values[1, 2]) +
#  scale_color_manual(name = "Wetland by Site", labels = c('Upland, Site A',
#       'Upland, Site B', 'Upland, Site C', 'Wetland, Site B', 'Wetland, Site C'),
#  values = my_colors)
#ggsave(file = '../figures/pc1pc2.pdf', width = 5, height = 4)
#
#ggplot(sample_data(physeq), aes(x = pc2, y = pc3, color = interaction(geocode, Wetland))) +
#  geom_point(size = 3) +
#  labs(x = paste0('PC 2 (', myround(pcoa$values[2, 2]*100, 1),'%)'), 
#       y = paste0('PC 3 (', myround(pcoa$values[3, 2]*100, 1),'%)')) +
#  coord_fixed(pcoa$values[3, 2]/pcoa$values[2, 2]) +
#  scale_color_manual(name = "Wetland by Site", labels = c('Upland, Site A',
#       'Upland, Site B', 'Upland, Site C', 'Wetland, Site B', 'Wetland, Site C'),
#  values = my_colors)
# # values = c('firebrick', 'red', 'pink', 'blue', 'lightblue'))
#ggsave(file = '../figures/pc2pc3.pdf', width = 5, height = 4)
#
#ggplot(sample_data(physeq), aes(x = pc1, y = pc3, color = interaction(geocode, Wetland))) +
#  geom_point(size = 3) +
#  labs(x = paste0('PC 1 (', myround(pcoa$values[1, 2]*100, 1),'%)'), 
#       y = paste0('PC 3 (', myround(pcoa$values[3, 2]*100, 1),'%)')) +
#  coord_fixed(pcoa$values[3, 2]/pcoa$values[1, 2]) +
#  scale_color_manual(name = "Wetland by Site", labels = c('Upland, Site A',
#       'Upland, Site B', 'Upland, Site C', 'Wetland, Site B', 'Wetland, Site C'),
#  values = my_colors)
# # values = c('firebrick', 'red', 'pink', 'blue', 'lightblue'))
#ggsave(file = '../figures/pc1pc3.pdf', width = 5, height = 4)

#+ trad_correlations

sample_data <- as(sample_data(physeq), 'data.frame')

sample_data(physeq) %>% 
  ggplot(mapping = aes(x = pmoa_copy_num, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = expression(italic('pmoA') ~ 'copy number'), y = expression('Methane oxidation (k)')) + 
  annotate(geom = 'text', x = 0, y = 2.5, label = 'A', size = 6)

ggsave(file = '../figures/lowk_pmoa.pdf', width = 4, height = 4)

sample_data(physeq) %>% 
  ggplot(mapping = aes(x = pc1, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'PC1', y = expression('Methane oxidation (k)')) + 
  annotate(geom = 'text', x = min(sample_data(physeq)$pc1),
           y = 2.5, label = 'B', size = 6)

ggsave(file = '../figures/lowk_pc1.pdf', width = 4, height = 4)

ggplot(sample_data(physeq), aes(x = richness, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Richness', y = expression('Methane oxidation (k)')) + 
  annotate(geom = 'text', x = min(sample_data(physeq)$richness), y = 2.5, label = 'C', size = 6)

ggsave(file = '../figures/lowk_rich.pdf', width = 4, height = 4)

# Low K and Shannon, w/ and w/o outlier

ggplot(sample_data(physeq), aes(x = shannon, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Shannon Diversity', y = expression('Methane oxidation (k)')) + 
  annotate(geom = 'text', x = min(sample_data(physeq)$shannon), y = 2.5, label = 'D', size = 6)

ggsave(file = '../figures/lowk_shan.pdf', width = 4, height = 4)

ggplot(sample_data[sample_data$shannon > 5, ], aes(x = shannon, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Shannon Diversity', y = expression('Methane oxidation (k)')) #+ 
#  annotate(geom = 'text', x = min(sample_data(physeq)[sample_data$shannon > 5, ]$shannon), y = 2.5, label = 'D', size = 6)


ggsave(file = '../figures/lowk_shan_no_outliers.pdf', width = 4, height = 4)

# Low k and Simpson, w/ and w/o outlier

ggplot(sample_data, aes(x = simpson, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "Simpson Diversity", y = expression('Methane oxidation (k)')) + 
  annotate(geom = 'text', x = min(sample_data(physeq)$simpson), y = 2.5, label = 'E', size = 6)

ggsave(file = '../figures/lowk_simp.pdf', width = 4, height = 4)

ggplot(sample_data[sample_data$simpson > 0.998, ], aes(x = simpson, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "Simpson Diversity", y = expression('Methane oxidation (k)')) #+ 
#  annotate(geom = 'text', x = min(sample_data(physeq)$simpson), y = 2.5, label = 'E', size = 6)


ggsave(file = '../figures/lowk_simp_no_outliers.pdf', width = 4, height = 4)

trad_cor_table <- data.frame(predictor = c('\\textit{pmoA} abundance',
                                           'PC1', 'Richness', 'Shannon Diversity',
                                           'Simpson Diversity'),
estimate = NA,
SE = NA,
t.value = NA,
p.value = NA,
r.squared = NA)

fit <- lm(Low_final_k ~ pmoa_copy_num, data = sample_data)
trad_cor_table[1, 2] <- summary(fit)$coefficients[2, 1]
trad_cor_table[1, 3] <- summary(fit)$coefficients[2, 2]
trad_cor_table[1, 4] <- summary(fit)$coefficients[2, 3]
trad_cor_table[1, 5] <- summary(fit)$coefficients[2, 4]
trad_cor_table[1, 6] <- summary(fit)$r.squared

fit <- lm(Low_final_k ~ pc1, data = sample_data)
trad_cor_table[2, 2] <- summary(fit)$coefficients[2, 1]
trad_cor_table[2, 3] <- summary(fit)$coefficients[2, 2]
trad_cor_table[2, 4] <- summary(fit)$coefficients[2, 3]
trad_cor_table[2, 5] <- summary(fit)$coefficients[2, 4]
trad_cor_table[2, 6] <- summary(fit)$r.squared

fit <- lm(Low_final_k ~ richness, data = sample_data)
trad_cor_table[3, 2] <- summary(fit)$coefficients[2, 1]
trad_cor_table[3, 3] <- summary(fit)$coefficients[2, 2]
trad_cor_table[3, 4] <- summary(fit)$coefficients[2, 3]
trad_cor_table[3, 5] <- summary(fit)$coefficients[2, 4]
trad_cor_table[3, 6] <- summary(fit)$r.squared

fit <- lm(Low_final_k ~ shannon, data = sample_data)
trad_cor_table[4, 2] <- summary(fit)$coefficients[2, 1]
trad_cor_table[4, 3] <- summary(fit)$coefficients[2, 2]
trad_cor_table[4, 4] <- summary(fit)$coefficients[2, 3]
trad_cor_table[4, 5] <- summary(fit)$coefficients[2, 4]
trad_cor_table[4, 6] <- summary(fit)$r.squared

fit <- lm(Low_final_k ~ simpson, data = sample_data)
trad_cor_table[5, 2] <- summary(fit)$coefficients[2, 1]
trad_cor_table[5, 3] <- summary(fit)$coefficients[2, 2]
trad_cor_table[5, 4] <- summary(fit)$coefficients[2, 3]
trad_cor_table[5, 5] <- summary(fit)$coefficients[2, 4]
trad_cor_table[5, 6] <- summary(fit)$r.squared

trad_cor_table$p.value <- 
  c(myround(trad_cor_table$p.value[1], 3),
    formatC(trad_cor_table$p.value[c(2)], format = "e", digits = 2),
    myround(trad_cor_table$p.value[3:5], 3))
trad_cor_table[c(2), 5] <- 
  cell_spec(trad_cor_table[c(2), 5],  "latex", bold = T)

trad_cor_table %>% 
  mutate(estimate = paste0(myround(estimate, 3), " (", myround(SE, 3), ")")) %>% 
  mutate(t.value = myround(t.value, 3)) %>% 
  select(-SE) %>% 
  kable('latex', booktabs = TRUE, escape = FALSE, digits = 3,
        caption = 'Aspects of microbial community structure. Significant p-values
        are bolded (p \\textless{} 0.05). n = 44 except for \\textit{pmoA} where n = 42',
        col.names = c('Variable', 'Estimate (SE)', 't.value', 'p.value', '$R^2$')) %>% 
  writeLines('../tables/trad_cor_table.tex')


#'
#' ## Summary of correlations
#'
#' There is a significant, but marginal difference between wetland and upland sites.
#' Vmax appears more strongly correlated with edaphic properties (moisture, nitrogen,
#' bulk density, carbon) and abundance than Low_final_k. In contrast, Low_final_k is
#' more strongly correlated with composition in the form of PC axis 1. Neither process
#' is correlated with species richness or diversity. 
#' 
#'
#' # Variance Component Analysis
#' 
#' This analysis fits linear mixed-effects models using the `varComp` package
#' version `r packageVersion('varComp')`. The fixed effect portion of each model
#'  is `Low final K ~ taxon abundance + location` with random effect
#' variance-covariance matrices for community similarity (from Bray-Curtis
#' distance) and environmental similarity (Euclidean distance) with all variables
#' centered and scaled to unit variance before calculating similarity metrics.
#' 
#' For each set of analyses (all data, wetland data, and upland data) asvs were
#' subset to include only asvs present in more than one sample (so those only present
#' in zero samples or one sample were removed).


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

dim(geo_sim)
dim(com_sim)
ggplot(mapping = aes(x = c(geo_dis), y = c(com_dis))) + 
  geom_point()
ggplot(mapping = aes(x = c(geo_dis), y = c(env_dis))) + 
  geom_point()

print(mantel(com_sim_bin, geo_sim))
print(mantel(com_sim, geo_sim))
print(mantel(com_sim_bin, env_sim))
print(mantel(com_sim, env_sim))
print(mantel(geo_sim, env_sim))

#' ## Geography Correlated with Community Composition
#'
#' Community and geography are correlated as demonstrated by the Mantel test for
#' the community similarity matrix and the geographic similarity matrix (see Mantel test output). 
#' Even after variance stabilization, the community similarity metric is skewed
#' towards one end of the variance space while geography uses up more of the 
#' similarity matrix space (see histograms). This results in geography
#' explaining virtually all of the variation in ecosystem function.
#' Since sites are clustered around three locations,
#' I coded geography as three factors (one for each cluster of sites)
#'Also, wetland/upland covaries with those three sites, but moisture content is
#'included in the environmental similarity matrix and should capture that variation
#' To insure that the covariation between geographic site and wetland/upland
#' isn't having a major influence on the results, I ran the models with all of the data
#' and then again within wetland and within upland

#' # Variance Components


#+ r variance components
Low_final_k <- sample_data(physeq)$Low_final_k
geocode <- sample_data(physeq)$geocode
com <- cholRoot(com_sim)
#geo <- geocode
geo <- cholRoot(geo_sim)
env <- cholRoot(env_sim)

var_comp_pvalues <- data.frame(model = 
                               c('CH4 \\textasciitilde 1 + com', 
                                 'CH4 \\textasciitilde 1 + geo', 
                                 'CH4 \\textasciitilde 1 + env', 
                                 'CH4 \\textasciitilde 1 + geo + com', 
                                 'CH4 \\textasciitilde 1 + geo + com', 
                                 'CH4 \\textasciitilde 1 + geo + env', 
                                 'CH4 \\textasciitilde 1 + geo + env', 
                                 'CH4 \\textasciitilde 1 + com + env', 
                                 'CH4 \\textasciitilde 1 + com + env', 
                                 'CH4 \\textasciitilde 1 + geo + com + env', 
                                 'CH4 \\textasciitilde 1 + geo + com + env', 
                                 'CH4 \\textasciitilde 1 + geo + com + env'),
           component = c('com', 'geo', 'env', 'geo', 'com', 'env', 'geo', 'com', 'env', 'env', 'geo', 'com'),
           varcomp = NA,
           se = NA,
           p.value = NA)

# Single
# Com
fit0 <- varComp(Low_final_k ~ 1)
fit <- varComp(Low_final_k ~ 1, random =  ~ com)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[1] <- comps$h2G
var_comp_pvalues$se[1] <- sqrt(comps$Varh2G)

var_comp_pvalues$p.value[1] = p.value(varComp.test(fit, fit0))

# Geo
fit <- varComp(Low_final_k ~ 1, random =  ~ geo)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[2] <- comps$h2G
var_comp_pvalues$se[2] <- sqrt(comps$Varh2G)
var_comp_pvalues$p.value[2] = p.value(varComp.test(fit, fit0))

# Env
fit <- varComp(Low_final_k ~ 1, random =  ~ env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[3] <- comps$h2G
var_comp_pvalues$se[3] <- sqrt(comps$Varh2G)
var_comp_pvalues$p.value[3] = p.value(varComp.test(fit, fit0))

# Geo Com
fit0 <- varComp(Low_final_k ~ 1, random =  ~ com)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[4] <- comps$h2GE
var_comp_pvalues$se[4] <- sqrt(comps$Varh2GE)
var_comp_pvalues$p.value[4] = p.value(varComp.test(fit, fit0))

fit0 <- varComp(Low_final_k ~ 1, random =  ~ geo)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[5] <- comps$h2G
var_comp_pvalues$se[5] <- sqrt(comps$Varh2G)
var_comp_pvalues$p.value[5] = p.value(varComp.test(fit, fit0))

# Geo Env
fit0 <- varComp(Low_final_k ~ 1, random =  ~ geo)
fit <- varComp(Low_final_k ~ 1, random =  ~ geo + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[6] <- comps$h2GE
var_comp_pvalues$se[6] <- sqrt(comps$Varh2GE)
var_comp_pvalues$p.value[6] = p.value(varComp.test(fit, fit0))

fit0 <- varComp(Low_final_k ~ 1, random =  ~ env)
fit <- varComp(Low_final_k ~ 1, random =  ~ geo + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[7] <- comps$h2G
var_comp_pvalues$se[7] <- sqrt(comps$Varh2G)
var_comp_pvalues$p.value[7] = p.value(varComp.test(fit, fit0))

# Com Env
fit0 <- varComp(Low_final_k ~ 1, random =  ~ env)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[8] <- comps$h2G
var_comp_pvalues$se[8] <- sqrt(comps$Varh2G)
var_comp_pvalues$p.value[8] = p.value(varComp.test(fit, fit0))

fit0 <- varComp(Low_final_k ~ 1, random =  ~ com)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[9] <- comps$h2GE
var_comp_pvalues$se[9] <- sqrt(comps$Varh2GE)
var_comp_pvalues$p.value[9] = p.value(varComp.test(fit, fit0))

# Geo Com Env
fit0 <- varComp(Low_final_k ~ 1, random =  ~ com + geo)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[10] <- comps$h2GE[2]
var_comp_pvalues$se[10] <- sqrt(comps$Varh2GE[2])
var_comp_pvalues$p.value[10] = p.value(varComp.test(fit, fit0))

fit0 <- varComp(Low_final_k ~ 1, random =  ~ com + env)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[11] <- comps$h2GE[1]
var_comp_pvalues$se[11] <- sqrt(comps$Varh2GE[1])
var_comp_pvalues$p.value[11] = p.value(varComp.test(fit, fit0))

fit0 <- varComp(Low_final_k ~ 1, random =  ~ geo + env)
fit <- varComp(Low_final_k ~ 1, random =  ~ com + geo + env)
comps <- varCompGE(fit)
var_comp_pvalues$varcomp[12] <- comps$h2G
var_comp_pvalues$se[12] <- sqrt(comps$Varh2G)
var_comp_pvalues$p.value[12] = p.value(varComp.test(fit, fit0))

var_comp_pvalues$varcomp <- paste0(myround(var_comp_pvalues$varcomp, 2), ' (',
                                   myround(var_comp_pvalues$se, 2), ')')

var_comp_pvalues$p.value <- 
  c(formatC(var_comp_pvalues$p.value[c(1)], format = "e", digits = 2), 
    myround(var_comp_pvalues$p.value[2:7], 3),
    formatC(var_comp_pvalues$p.value[c(8)], format = "e", digits = 2),
    myround(var_comp_pvalues$p.value[9:12], 3))
var_comp_pvalues[c(1, 2, 3, 7, 8), 5] <- 
  cell_spec(var_comp_pvalues[c(1, 2, 3, 7, 8), 5],  "latex", bold = T)
#var_comp_pvalues[12, 5] <- 
#  cell_spec(var_comp_pvalues[12,5],  "latex", italic = T)
var_comp_pvalues <- var_comp_pvalues[, -4]

var_comp_pvalues %>% 
  kable('latex', booktabs = T, digits = 2, escape = F,
        caption = paste0('Variance component estimates from intercept-only models predicting the rate of high-affinity methane oxidation from soil. p-values
        determined by linear score test comparing nested models with and without
        each variance component. Significant p-values are bolded (p \\textless{} 0.05). 
        n = 44.'), 
        col.names = c('Model', 'Component', 'Estimate (SE)', 'p-value')) %>% 
  collapse_rows(columns = 1) %>% 
  writeLines('../tables/var_comp_pvalues.tex')

#' # Analyze model output

#+ r analyze_model_output
all_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_abund.csv'))
#all_model_output <- as.data.frame(data.table::fread('../output/var_comp_fam.csv'))
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

#' ## Number of taxa identified with each covariate

for (i in unique(all_model_output$comps)) {
  print(i)
my_qs <- qvalue(all_model_output$p.value[all_model_output$comps == i])
summary(my_qs)
all_model_output$qvalues[all_model_output$comps == i] <- my_qs$qvalues
}
ggplot(all_model_output, aes(x = p.value)) +
  facet_wrap(~ comps) +
  geom_histogram()

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



