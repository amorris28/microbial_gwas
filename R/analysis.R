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
#library(broom)
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
library(DivNet)
library(ape)

#+ r options
options(knitr.kable.NA = '')
theme_set(theme_classic())

#+ r import_data
physeq <- readRDS('../output/physeq.rds')
physeq <- 
physeq %>% 
  subset_samples(Experiment == 'MO')
# Remove 0 abundance taxa
#otu_table(physeq) <- otu_table(physeq)[, colSums(otu_table(physeq)) > 0]

# Add variables
com_mat <- otu_table(physeq)
sample_data(physeq)$richness <- rowSums(com_mat > 0)
sample_data(physeq)$shannon <- diversity(com_mat, "shannon")
sample_data(physeq)$simpson <- diversity(com_mat, "simpson")

pca <- rda(otu_table(physeq), scale = TRUE)
head(summary(pca))
sample_data(physeq)$pc1 <- scores(pca, choices = 1, display = 'sites', scaling = 3)
sample_data(physeq)$pc2 <- scores(pca, choices = 2, display = 'sites', scaling = 3)
sample_data(physeq)$pc3 <- scores(pca, choices = 3, display = 'sites', scaling = 3)
sample_data(physeq)$pc4 <- scores(pca, choices = 4, display = 'sites', scaling = 3)
sample_data(physeq)$pc5 <- scores(pca, choices = 5, display = 'sites', scaling = 3)
sample_data(physeq)$pc6 <- scores(pca, choices = 6, display = 'sites', scaling = 3)
sample_data(physeq)$pc7 <- scores(pca, choices = 7, display = 'sites', scaling = 3)

sample_data(physeq) %>% 
  filter(Experiment == 'MO') %>% 
  ggplot(aes(Low_final_k)) +
  geom_density()

sample_data(physeq) %>% 
  filter(Experiment == 'MO') %>% 
  ggplot(aes(Low_final_k)) +
  facet_grid(Wetland ~ geocode) +
  geom_density()


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
all_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_out.csv'))

#' ## Variance components

#+ r variance_components
all_model_output %>% 
  group_by(comps) %>% 
  summarize_at(vars(geo, com, env, err), mean) %>% 
  kable(caption = "Variance components for all data",
        col.names = c('Components', 'Geography', 'Community', 'Environment', 'Error'),
        digits = 2) %>% 
  kable_styling()

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

n_taxa <- 
sample_data(physeq) %>% 
  left_join(tax_table(physeq)) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05)  %>% 
  summarize(count = length(unique(asv)))

my_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / length(unique(asv)))  %>% 
  summarize(count = length(unique(asv))) %>% 
  kable(caption = "N taxa identified with all data",
        col.names = c('Components', 'Count')) %>% 
  kable_styling()

wet_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / length(unique(asv)))  %>% 
  summarize(count = length(unique(asv))) %>% 
  kable(caption = "N taxa identified with wetland data",
        col.names = c('Components', 'Count')) %>% 
  kable_styling() 

n_taxa <- 
wet_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / length(unique(asv)))  %>% 
  summarize(count = length(unique(asv)))

dry_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05)  %>% 
  summarize(count = length(unique(asv))) %>%
  kable('latex', booktabs = TRUE, caption = "N taxa identified with upland data",
        col.names = c('Components', 'Count'))


n_taxa <- 
dry_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05)  %>% 
  summarize(count = length(unique(asv)))

n_taxa <- 
dry_model_output_geo %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05)  %>% 
  summarize(count = length(unique(asv)))

taxon_table <- data.table::fread('../output/taxon_table.csv', data.table = FALSE, na.strings = "")
model_output_abund %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05)  %>% 
  summarize(count = length(unique(asv)))

#' ## Taxa identified with each covariate



#+ r plot_shared_taxa, eval=FALSE
plot_asv_euler <- function(model_output) {
asvs <- 
model_output %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / length(unique(asv)))  %>% 
  select(comps, asv) %>% 
  as.data.frame()
asvs <- list(
                 'g' = asvs[asvs$comps == 'g', 2],
                 'c' = asvs[asvs$comps == 'c', 2],
                 'e' = asvs[asvs$comps == 'e', 2],
                 'n' = asvs[asvs$comps == 'n', 2],
                 'gc' = asvs[asvs$comps == 'gc', 2],
                 'ge' = asvs[asvs$comps == 'ge', 2],
                 'ce' = asvs[asvs$comps == 'ce', 2],
                 'gce' = asvs[asvs$comps == 'gce', 2]
                 )
plot(euler(asvs))
}

#' ### All Data Model

#+ r all_taxa

model_output_abund %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05 ) %>% 
  select(comps, asv, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  ungroup()

my_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / (length(unique(asv)))) %>% 
  select(comps, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  kable() %>% 
  kable_styling() %>% 
  collapse_rows(columns = 1)

sig_taxa_all <- 
all_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(qvalues < 0.05 ) %>% 
  select(comps, asv, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  ungroup()

  sig_taxa_gce <- 
all_data %>% 
  select(Low_final_k, Land_type, filter(sig_taxa_all, comps == 'gce')$asv)

head(sig_taxa_gce)


ggplot(sig_taxa_gce, aes(


sig_pca <- 
  sig_taxa_gce %>% 
  select(starts_with('asv')) %>% 
  as.matrix() %>% 
  rda(scale = TRUE)
biplot(sig_pca)

sig_taxa_gce$sig_pc1 <- 
  scores(sig_pca, choices = 1, display = 'sites', scaling = 3)
sig_taxa_gce$sig_pc2 <- 
  scores(sig_pca, choices = 2, display = 'sites', scaling = 3)

ggplot(sig_taxa_gce, aes(sig_pc1, Low_final_k)) +
  geom_point()
ggplot(sig_taxa_gce, aes(sig_pc2, Low_final_k)) +
  geom_point()

  scores(sig_pca, choices = 1, display = 'species', scaling = 3)
  scores(sig_pca, choices = 2, display = 'species', scaling = 3)

summary(lm(Low_final_k ~ sig_pc1, all_data))

fit <- varComp(Low_final_k ~ as.matrix(sig_taxa), data = all_data, varcov = list(com_sim, geo_sim, env_sim))
summary(fit)


colnames(sig_taxa_all)
taxa <- 
sig_taxa_all %>% 
  select(comps, Phylum:Genus)

sig_taxa_all %>% 
  filter(comps == 'gce')  %>% 
  select(Phylum:Genus) %>% 
  kable('latex', booktabs = T, caption = 'Taxa significantly correlated with
        high-affinity methane oxidation rate after controlling for geographic
        proximity, environmental similarity, and community structure. 
        Significance determined by controlling
        the false discovery rate (q \\textless{} 0.05).') %>% 
  writeLines('../tables/all_taxa_sig.tex')


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



