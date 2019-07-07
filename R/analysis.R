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
#library(lubridate)
library(broom)
#library(broman)
#library(morris)
library(varComp)
library(vegan)
library(parallel)
library(qvalue)
library(knitr)
library(kableExtra)
#library(eulerr)
#library(GGally)
library(phyloseq)
library(doParallel)
library(gap)
library(lmtest)

num.cores <- detectCores()
registerDoParallel(cores = 3)

#+ r tree_analysis
tree <- read_tree('../data/gabon/16S_rep_seqs_tree_Gabon.nwk')
#+ r import_data

taxon_table <- data.table::fread('../output/taxon_table.csv', data.table = FALSE)
all_data <- data.table::fread('../output/gab_all_vst_troph.csv', data.table = FALSE)
n_samp <- nrow(all_data)
n_asvs <- ncol(all_data[, grepl("^[asv]", colnames(all_data))])

## Remove ASVs only present in 1 sample
asvs <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
vars <- all_data[, !grepl("^[asv]", colnames(all_data))]
asvs <- asvs[, colSums(asvs) > 0]
all_data <- cbind(vars, asvs)

# Add location factor
all_data[all_data$Y > 30000, 'geocode'] <- 'A'
all_data[all_data$Y < 30000 & all_data$Y > 10000, 'geocode'] <- 'B'
all_data[all_data$Y < 10000, 'geocode'] <- 'C'
all_data$geocode <- factor(all_data$geocode)

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

colnames(all_data[, 1:13])

print(ggplot(all_data, aes(x = X, y = Y)) + 
  geom_point(aes(color = all_data$Wetland)) +
  labs(title = "Distribution of sites highlighted by wetland or upland"))

ggplot(all_data, aes(x = Vmax, y = Low_final_k)) + geom_point()

ggplot(all_data, aes(x = pmoa_copy_num, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Low_final_k ~ pmoa_copy_num, data = all_data))

ggplot(all_data, aes(x = pmoa_copy_num, y = Vmax)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Vmax ~ pmoa_copy_num, data = all_data))

adonis(asvs ~ all_data$Wetland)

pca <- rda(as.matrix(all_data[, grepl("^[asv]", colnames(all_data))]))
screeplot(pca)
head(summary(pca))
all_data$pc1 <- scores(pca, choices = 1, display = 'sites', scaling = 3)
all_data$pc2 <- scores(pca, choices = 2, display = 'sites', scaling = 3)
all_data$pc3 <- scores(pca, choices = 3, display = 'sites', scaling = 3)

ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = pc1, y = pc2))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = pc1, y = pc3))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = pc2, y = pc3))

ggplot(all_data, aes(color = geocode)) +
  geom_point(aes(x = pc1, y = pc2))
ggplot(all_data, aes(color = geocode)) +
  geom_point(aes(x = pc1, y = pc3))
ggplot(all_data, aes(color = geocode)) +
  geom_point(aes(x = pc2, y = pc3))

ggplot(all_data[all_data$pc1 < 0,], aes(x = pc1, y = Low_final_k)) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(Low_final_k ~ pc1, data = all_data))

ggplot(all_data, aes(x = pc1, y = Vmax)) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(Vmax ~ pc1, data = all_data))

#my_asvs <- tibble(asv = names(scores(pca)$species[, 'PC1']),
#          score = scores(pca)$species[, 'PC1'])
#taxon_table  %>% 
#  left_join(my_asvs) %>% 
#  filter(abs(score) > 0.5) %>% 
#  kable(digits = 2) %>% 
#  kable_styling()

ggplot(all_data, aes(x = pc2, y = Low_final_k)) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(Low_final_k ~ pc2, data = all_data))

env_pca <- rda(as.matrix(all_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')]))
screeplot(env_pca)
head(summary(env_pca))
all_data$env_pc1 <- scores(env_pca, choices = 1, display = 'sites', scaling = 3)
all_data$env_pc2 <- scores(env_pca, choices = 2, display = 'sites', scaling = 3)
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = env_pc1, y = env_pc2))
ggplot(all_data, aes(color = geocode)) +
  geom_point(aes(x = env_pc1, y = env_pc2))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = env_pc1, y = Mois_cont))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = env_pc1, y = C_percent))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = env_pc1, y = N_percent))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = env_pc1, y = Bulk_dens))
ggplot(all_data, aes(x = Wetland, y = env_pc1)) +
  geom_boxplot() + 
  geom_jitter()
summary(lm(env_pc1 ~ Wetland, all_data))
ggplot(all_data, aes(x = Wetland, y = Mois_cont)) +
  geom_boxplot() + 
  geom_jitter(aes(color = Site))
summary(lm(Mois_cont ~ Wetland, all_data))
## Spatial PCA

geo_pca <- rda(as.matrix(all_data[, c('X', 'Y')]))
screeplot(geo_pca)
summary(geo_pca)
all_data$geo_pc1 <- scores(geo_pca, choices = 1, display = 'sites', scaling = 3)

ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = geo_pc1, y = Low_final_k))
ggplot(all_data, aes(color = Wetland)) +
  geom_point(aes(x = geo_pc1, y = Mois_cont))

ggplot(all_data, aes(x = Wetland, y = Low_final_k)) +
       geom_boxplot() +
       geom_jitter()
summary(lm(Low_final_k ~ Wetland, data = all_data))

ggplot(all_data, aes(x = Mois_cont, y = Low_final_k)) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(Vmax ~ Mois_cont, data = all_data))

com_matrix <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
all_data$richness <- rowSums(com_matrix > 0)
all_data$shannon <- diversity(as.matrix(all_data[, grepl("^[asv]", colnames(all_data))]), "shannon")
all_data$simpson <- diversity(as.matrix(all_data[, grepl("^[asv]", colnames(all_data))]), "simpson")

ggplot(data = all_data, aes(x = richness, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Low_final_k ~ richness, data = all_data))

ggplot(data = all_data, aes(x = richness, y = Vmax)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Vmax ~ richness, data = all_data))

# Low K and Shannon, w/ and w/o outlier

ggplot(data = all_data, aes(x = shannon, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Low_final_k ~ shannon, data = all_data))

ggplot(data = all_data[all_data$shannon > 5, ], aes(x = shannon, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  labs(title = "Without sample 26")
summary(lm(Low_final_k ~ shannon, data = all_data[all_data$shannon > 5, ]))

# V max and Shannon, w/ and w/o outlier

ggplot(data = all_data, aes(x = shannon, y = Vmax)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Vmax ~ shannon, data = all_data))

ggplot(data = all_data[all_data$shannon > 6, ], aes(x = shannon, y = Vmax)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  labs(title = "Without sample 1 and 26")
summary(lm(Vmax ~ shannon, data = all_data[all_data$shannon > 6, ]))

# Low k and Simpson, w/ and w/o outlier

ggplot(data = all_data, aes(x = simpson, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Low_final_k ~ simpson, data = all_data))

ggplot(data = all_data[all_data$simpson > 0.998, ], aes(x = simpson, y = Low_final_k)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  labs(title = "Without sample 1 and 26")
summary(lm(Low_final_k ~ simpson, data = all_data[all_data$simpson > 0.998, ]))

# Vmax and simpson, w/ and w/o outlier

ggplot(data = all_data, aes(x = simpson, y = Vmax)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(Vmax ~ simpson, data = all_data))

ggplot(data = all_data[all_data$simpson > 0.998, ], aes(x = simpson, y = Vmax)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  labs(title = "Without sample 1 and 26")
summary(lm(Vmax ~ simpson, data = all_data[all_data$simpson > 0.998, ]))



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

#+ r sim_mat

calc_sim_matrix <- function(x) {
1/(1+as.matrix(dist(as.matrix(scale(x)))))
}


# Community similarity

x <- 3
asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
com_sim <- 1 - as.matrix(vegdist(asvs[, -c(x)], method = "bray"))
com_sim_bin <- 1 - as.matrix(vegdist(asvs[, -c(x)], method = "jaccard", binary = TRUE))

# Geographic similarity
geo_sim <- calc_sim_matrix(all_data[, c('X', 'Y')])
with(all_data, qplot(X, Y))
qplot(c(geo_sim))

# environmental similarity'
env_sim <- calc_sim_matrix(all_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])

#' ## Collinearity

#+ collinear_all_data
with(all_data, cor(pc1, env_pc1))
with(all_data, cor(pc3, env_pc1))
anova(lm(pc1 ~ geocode, all_data))
anova(lm(env_pc1 ~ geocode, all_data))
anova(lm(env_pc1 ~ Wetland, all_data))
anova(lm(pc1 ~ Wetland, all_data))
print(mantel(com_sim_bin, geo_sim))
print(mantel(com_sim, geo_sim))
print(mantel(com_sim_bin, env_sim))
print(mantel(com_sim, env_sim))
print(mantel(geo_sim, env_sim))

#' ## Model diagnostics for all data

#+ r varcomp_diagnostics

varcomp_diag <- function(data) {
print(mantel(com_sim, geo_sim))
print(mantel(com_sim, env_sim))
print(mantel(geo_sim, env_sim))


print(ggplot(data, aes(x = X, y = Y)) + 
  geom_point(aes(color = data$Wetland)) +
  labs(title = "Distribution of sites highlighted by wetland or upland"))
print(ggplot(data, aes(x = X, y = Y)) + 
  geom_point(aes(color = data$geocode)) +
  labs(title = "Distribution of sites highlighted by site code"))

print(ggplot(data, aes(Low_final_k)) + 
  geom_histogram(bins = 20) + 
  labs(x = "Low Final K", y = "Number of Samples"))
print(ggplot(data, aes(asv)) + 
  geom_histogram(bins = 20) + 
  labs(x = "Taxon Abundance", y = "Number of Samples",
  title = paste0('Abundance of ', asv_tax)))


print(ggplot(mapping = aes(c(com_sim))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Community Similarity (Bray-Curtis)", y = "Number of Samples"))
print(ggplot(mapping = aes(c(com_sim_bin))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Community Similarity (Jaccard)", y = "Number of Samples"))
print(ggplot(mapping = aes(c(geo_sim))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Spatial Similarity (Euclidean)", y = "Number of Samples"))
print(ggplot(mapping = aes(c(env_sim))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Environment Similarity (Euclidean)", y = "Number of Samples"))

model <- varComp(Low_final_k ~ asv, data = data,
                 varcov = list(com = com_sim, env = env_sim, geo = geo_sim) )
summary(model)
}
varcomp_diag(all_data)

#+ r com_var_comp_test

fit <- varComp(Low_final_k ~ 1, all_data, varcov = list(com = com_sim))
fit1 <- varComp(Low_final_k ~ 1, all_data)

varComp.test(fit, fit1)
lrtest(fit, fit1)
summary(fit)
h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

h2 <- h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
var_comps_se <- data.frame(comp = c('com'), h2 = c(h2$h2G), 
           SE = c(sqrt(h2$Varh2G)))
ggplot(var_comps_se, aes(x = comp, y = h2, ymax = h2 + SE, ymin = h2 - SE)) +
  geom_pointrange() 


fit <- varComp(Low_final_k ~ 1, all_data, varcov = list(com = com_sim_bin))
fit1 <- varComp(Low_final_k ~ 1, all_data)

varComp.test(fit, fit1)
lrtest(fit, fit1)
summary(fit)
h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

h2 <- h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
var_comps_se <- data.frame(comp = c('com'), h2 = c(h2$h2G), 
           SE = c(sqrt(h2$Varh2G)))
ggplot(var_comps_se, aes(x = comp, y = h2, ymax = h2 + SE, ymin = h2 - SE)) +
  geom_pointrange() 

fit <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim, env = env_sim))
fit1 <- varComp(Low_final_k ~ geocode, all_data, varcov = list(env = env_sim))

varComp.test(fit, fit1)
lrtest(fit, fit1)
summary(fit)
h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

h2 <- h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
var_comps_se <- data.frame(comp = c('com', 'env'), h2 = c(h2$h2G, h2$h2GE), 
           SE = c(sqrt(h2$Varh2G), sqrt(h2$Varh2GE)))
ggplot(var_comps_se, aes(x = comp, y = h2, ymax = h2 + SE, ymin = h2 - SE)) +
  geom_pointrange() 


fit <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim_bin, env = env_sim))
fit1 <- varComp(Low_final_k ~ geocode, all_data, varcov = list(env = env_sim))

varComp.test(fit, fit1)
lrtest(fit, fit1)
summary(fit)
h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

h2 <- h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
var_comps_se <- data.frame(comp = c('com', 'env'), h2 = c(h2$h2G, h2$h2GE), 
           SE = c(sqrt(h2$Varh2G), sqrt(h2$Varh2GE)))
ggplot(var_comps_se, aes(x = comp, y = h2, ymax = h2 + SE, ymin = h2 - SE)) +
  geom_pointrange() 


fit <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim, env = env_sim))
fit1 <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim))

varComp.test(fit, fit1)


fit <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim_bin, env = env_sim))
fit1 <- varComp(Low_final_k ~ 1, all_data, varcov = list(com = com_sim, env = env_sim))

varComp.test(fit, fit1)
h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
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

#+ r h2ge_attempt

#
#h2GE(all_data$Low_final_k, com_sim)
#h2GE
#h2G(all_data$Low_final_k, com_sim)
#h2G(all_data$Low_final_k, env_sim)
#h2G(all_data$Low_final_k, geo_sim)
#h2GE(c(model$varComps, model$sigma2))
#h2GE
#cov(c(model$varComps, model$sigma2))
#
#model <- varComp(Low_final_k ~ taxon + geocode, data = all_data,
#                 varcov = list(com = com_sim, env = env_sim) )
#anova(model)$`Pr(>F)`[2]
#calc_varcomp_per(model)
#
#model <- varComp(Low_final_k ~ taxon + geocode, data = all_data, 
#                 varcov = list(com = com_sim, env = env_sim) )
#summary(model)
#
#levels(all_data$geocode)
#all_data$geocode <- factor(all_data$geocode, levels(all_data$geocode)[c(2,3,1)])
#model1 <- varComp(lowk ~ taxon + all_data$geocode, varcov = list(com = com_sim, env = env_sim) )
#summary(model1)
#
#
#
#anova(model, model2)
#anova(model)


# ## Model diagnostics for Wetland data

#+ cor_using_varcomp, eval=FALSE
is.na(all_data$pmoa_copy_num)
fake_data <- all_data
fake_data$pmoa_copy_num[is.na(all_data$pmoa_copy_num )] <- mean(all_data$pmoa_copy_num, na.rm = T)
model <- varComp(Low_final_k ~ pmoa_copy_num + geocode,
                 data = fake_data, na.action = na.exclude,
                 varcov = list(com = com_sim, env = env_sim) )
summary(model)
#+ r varcomp_diagnostics_wet, eval=FALSE

wet_data <- filter(all_data, Wetland == 'Wetland')

## Remove ASVs only present in 1 sample
asvs <- as.matrix(wet_data[, grepl("^[asv]", colnames(wet_data))])
vars <- wet_data[, !grepl("^[asv]", colnames(wet_data))]
asvs <- asvs[, colSums(asvs) > 0]
wet_data <- cbind(vars, asvs)

# Geographic similarity
geo_sim <- calc_sim_matrix(wet_data[, c('X', 'Y')])

# environmental similarity'
env_sim <- calc_sim_matrix(wet_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])

asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))

varcomp_diag(wet_data)

# ## Model diagnostics for Upland data

#+ r varcomp_diagnostics_dry, eval=FALSE

dry_data <- filter(all_data, Wetland == 'Upland')

## Remove ASVs only present in 1 sample
asvs <- as.matrix(dry_data[, grepl("^[asv]", colnames(dry_data))])
vars <- dry_data[, !grepl("^[asv]", colnames(dry_data))]
asvs <- asvs[, colSums(asvs) > 0]
dry_data <- cbind(vars, asvs)

# Geographic similarity
geo_sim <- calc_sim_matrix(dry_data[, c('X', 'Y')])

# environmental similarity'
env_sim <- calc_sim_matrix(dry_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])

asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
com_sim_bin <- 1 - as.matrix(vegdist(asvs[, -x], method = 'jaccard', binary = TRUE))

pca <- rda(as.matrix(dry_data[, grepl("^[asv]", colnames(dry_data))]))
screeplot(pca)

dry_data$pc1 <- scores(pca, choices = 1, display = 'sites', scaling = 3)
dry_data$pc2 <- scores(pca, choices = 2, display = 'sites', scaling = 3)
dry_data$pc3 <- scores(pca, choices = 3, display = 'sites', scaling = 3)

env_pca <- rda(as.matrix(dry_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')]))
screeplot(env_pca)
summary(env_pca)
dry_data$env_pc1 <- scores(env_pca, choices = 1, display = 'sites', scaling = 3)
dry_data$env_pc2 <- scores(env_pca, choices = 2, display = 'sites', scaling = 3)

#+ collinear_dry_data, eval=FALSE
with(dry_data, cor(pc1, env_pc1))
with(dry_data, cor(pc2, env_pc1))
with(dry_data, cor(pc3, env_pc1))
anova(lm(pc1 ~ geocode, dry_data))
anova(lm(env_pc1 ~ geocode, dry_data))
print(mantel(com_sim_bin, geo_sim))
print(mantel(com_sim, geo_sim))
print(mantel(com_sim_bin, env_sim))
print(mantel(com_sim, env_sim))
print(mantel(geo_sim, env_sim))

varcomp_diag(dry_data)

#' # Analyze model output


#+ r analyze_model_output




my_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_out.csv'))
wet_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_out_wet.csv'))
dry_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_out_dry.csv'))


taxon_table$taxon <- apply(taxon_table, 1, function(r) { r[ 
                      if( any(is.na(r))){ 
                          (which(is.na(r))[1]-1) 
                       }else{
                          (length(r)-0)}
                                      ] }
       )
hist(my_model_output[my_model_output$comps == 'e', 'p.value'])
head(my_model_output)
myqs <- qvalue(my_model_output[my_model_output$comps == 'gce', 'p.value'], fdr.level = 0.1)
myqs$pi0
summary(myqs)
hist(myqs)
plot(myqs)
myqs
filter(taxon_table, asv %in% my_model_output[my_model_output$comps == 'e' & myqs$significant == TRUE, 'asv'])

#' ## Variance components

#+ r variance_components
my_model_output %>% 
  group_by(comps) %>% 
  summarize_at(vars(com, env, err), mean) %>% 
  kable(caption = "Variance components for all data",
        col.names = c('Components', 'Community', 'Environment', 'Error'),
        digits = 2) %>% 
  kable_styling()

wet_model_output %>% 
  group_by(comps) %>% 
  summarize_at(vars(com, env, err), mean) %>% 
  kable(caption = "Variance components for wetland data",
        col.names = c('Components', 'Community', 'Environment', 'Error'),
        digits = 2) %>% 
  kable_styling()

dry_model_output %>% 
  group_by(comps) %>% 
  summarize_at(vars(com, env, err), mean) %>% 
  kable(caption = "Variance components for upland data",
        col.names = c('Components', 'Community', 'Environment', 'Error'),
        digits = 2) %>% 
  kable_styling()

#' ## Number of taxa identified with each covariate

#+ r n_taxa_ided

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

dry_model_output %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / length(unique(asv)))  %>% 
  summarize(count = length(unique(asv))) %>% 
  kable(caption = "N taxa identified with upland data",
        col.names = c('Components', 'Count')) %>% 
  kable_styling()

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
my_model_output  %>% 
  left_join(taxon_table) %>% 
group_by(comps) %>% 
  filter(p.value <= 0.05 / (length(unique(asv)))) %>% 
  select(comps, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  kable() %>% 
  kable_styling() %>% 
  collapse_rows(columns = 1)


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


#' ### Upland Model

#+ r upland_taxa

dry_model_output  %>% 
  left_join(taxon_table) %>% 
  group_by(comps) %>% 
  filter(p.value <= 0.05 / (length(unique(asv)))) %>% 
  select(comps, Phylum:Genus) %>% 
  arrange(comps, Phylum, Class, Order, Family, Genus) %>% 
  kable() %>% 
  kable_styling() %>% 
  collapse_rows(columns = 1)



#+ phylo_tree, eval=FALSE

x <- 3
asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
com_sim <- 1 - as.matrix(vegdist(asvs[, -c(x)]))

# Geographic similarity
geo_sim <- calc_sim_matrix(all_data[, c('X', 'Y')])

# environmental similarity'
env_sim <- calc_sim_matrix(all_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])
asvs <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
asv <- asvs[, x]
colnames(asvs) <- substr(colnames(asvs), 5, nchar(colnames(asvs)[1]))
otu_table <- otu_table(asvs, taxa_are_rows = FALSE)
phylo_object <- phyloseq(otu_table, tree)
unifrac <- UniFrac(phylo_object, TRUE)
uni_sim <- 1/(1+as.matrix(unifrac))
ggplot(mapping = aes(c(uni_sim))) +
  geom_histogram(bins = 20)


ggplot(mapping = aes(x = c(com_sim), y = c(uni_sim))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
model <- varComp(Low_final_k ~ 1, data = all_data,
                 varcov = list(com = uni_sim) )
summary(model)

h2G(c(model$varComps['com'], err = model$sigma2), vcov(model, what = 'varComp'))

model <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = uni_sim) )
summary(model)

h2G(c(model$varComps['com'], err = model$sigma2), vcov(model, what = 'varComp'))

model <- varComp(Low_final_k ~ asv + geocode, data = all_data,
                 varcov = list(com = uni_sim, env = env_sim) )
all_data$geocode
class(uni_sim)
summary(model)
vcov(model, what = 'varComp')
h2 <- h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp'))

var_comps_se <- data.frame(comp = c('uni', 'env'), h2 = c(h2$h2G, h2$h2GE), 
           SE = c(sqrt(h2$Varh2G), sqrt(h2$Varh2GE)))
ggplot(var_comps_se, aes(x = comp, y = h2, ymax = h2 + SE, ymin = h2 - SE)) +
  geom_pointrange() 
sessionInfo()

#+ exploratory_stuff, eval=FALSE

with(all_data, cor(pc1, env_pc1))
with(all_data, cor(pc3, env_pc1))
anova(lm(pc1 ~ geocode, all_data))
anova(lm(env_pc1 ~ geocode, all_data))
anova(lm(env_pc1 ~ Wetland, all_data))
anova(lm(pc1 ~ Wetland, all_data))
fit <- varComp(Low_final_k ~ asv, dry_data, ~ pc1 + env_pc1 + geocode)
summary(fit)
h2G(c(fit$varComps['pc1'], err = fit$sigma2), vcov(fit, what = 'varComp')[c(1,4), c(1,4)])
V <- c(fit$varComps, err = fit$sigma2)
VCOV <- vcov(fit, what = 'varComp', drop = FALSE)
h2GE(V, VCOV)

fit2 <- varComp(Low_final_k ~ asv, all_data, ~ env_pc1 + geocode)
fit3 <- varComp(Low_final_k ~ asv, all_data, ~ pc1 + geocode)
fit4 <- varComp(Low_final_k ~ asv, all_data, ~ env_pc1 + pc1 )
fit5 <- varComp(Low_final_k ~ asv, all_data, ~  pc1 )
fit6 <- varComp(Low_final_k ~ asv, all_data, ~  pc2 )
fit7 <- varComp(Low_final_k ~ asv, all_data)
fit8 <- varComp(Low_final_k ~ asv, all_data, ~ geocode)
summary(fit)
lrtest(fit, fit2)
lrtest(fit, fit3)
lrtest(fit, fit4)
lrtest(fit4, fit5)
lrtest(fit7, fit8)


model <- varComp(Low_final_k ~ asv + geocode, data = all_data,
                 varcov = list(com = com_sim, env = env_sim) )

h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp'))

model <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim, geo = geo_sim, env = env_sim) )
model1 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(geo = geo_sim, env = env_sim) )
model2 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim, env = env_sim) )
model3 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim, geo = geo_sim) )

lrtest(model, model1)
lrtest(model, model2)
lrtest(model, model3)

model <- varComp(Low_final_k ~ asv, data = all_data)
model1 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(geo = geo_sim) )
model2 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(env = env_sim) )
model3 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim) )

lrtest(model, model1)
lrtest(model, model2)
lrtest(model, model3)

model <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim) )
model1 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(geo = geo_sim, env = env_sim) )
model2 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(geo = geo_sim, com = com_sim) )
model3 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim) )

lrtest(model, model1)
lrtest(model2, model)
lrtest(model, model3)

model3 <- varComp(Low_final_k ~ asv, data = all_data)

lrtest(model2, model3)


model2 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim_bin) )

model3 <- varComp(Low_final_k ~ asv, data = all_data)

lrtest(model2, model3)

fit <- varComp(Low_final_k ~ asv + geocode, all_data, ~ pc1 + env_pc1 + geocode/Wetland)
summary(fit)
fit1 <- varComp(Low_final_k ~ asv, all_data, ~ env_pc1 + geo_pc1)
summary(fit1)
lrtest(fit, fit1)

h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

fit <- varComp(Low_final_k ~ asv, all_data, ~ pc1 + geocode + env_pc1)
summary(fit)
fit1 <- varComp(Low_final_k ~ asv, all_data, ~ geocode + env_pc1)
summary(fit1)
lrtest(fit, fit1)

h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

fit <- varComp(Low_final_k ~ asv + geocode, all_data, ~ pc1 + Wetland + env_pc1)
summary(fit)
fit1 <- varComp(Low_final_k  ~ asv + geocode, all_data, ~ Wetland + env_pc1)
summary(fit1)
lrtest(fit, fit1)
ps <- foreach(asv = seq_len(ncol(asvs)), .combine = rbind) %dopar% {
cbind(asv = colnames(asvs)[asv], tidy(lm(Low_final_k ~ asvs[, asv], all_data))[2, 5])
}
nrow(ps)
ps[ps$p.value < 0.000000001,]
tidy(lm(Low_final_k ~ asv, all_data))[2, 5]

com_sim <- 1 - as.matrix(vegdist(asvs[, -1], method = "bray"))
com_sim_bin <- 1 - as.matrix(vegdist(asvs[, -1], method = "jaccard", binary = TRUE))
fit <- varComp(Low_final_k ~ asvs[, 1], all_data, 
varcov = list(com = com_sim, env = env_sim))
summary(fit)
h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

com_sim <- 1 - as.matrix(vegdist(asvs[, -6203], method = "bray"))
com_sim_bin <- 1 - as.matrix(vegdist(asvs[, -6203], method = "jaccard", binary = TRUE))
fit <- varComp(Low_final_k ~ asvs[, 6203], all_data, 
varcov = list(com = com_sim, env = env_sim))
summary(fit)
h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

asv <- 1
com_sim <- 1 - as.matrix(vegdist(asvs[, -asv], method = "bray"))
model <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim))
summary(model)
h2G(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))
model_n <- varComp(Low_final_k ~ geocode, all_data)
summary(model_n)
lrtest(model, model_n)

model1 <- varComp(Low_final_k ~ asvs[, asv] + geocode, all_data, varcov = list(jom = com_sim))
summary(model1)
h2G(c(model1$varComps, err = model1$sigma2), vcov(model1, what = 'varComp', drop = FALSE))
model1_n <- varComp(Low_final_k ~ asvs[, asv] + geocode, all_data)
summary(model1_n)
lrtest(model1, model1_n)



ggplot(all_data, aes(x = Mois_cont, y = Low_final_k)) +
#  facet_wrap(~ Wetland, scales = 'free') +
  geom_point() + 
  stat_smooth(method = 'lm')
summary(lm(Low_final_k ~ Mois_cont + Wetland, all_data))
summary(lm(Low_final_k ~ Mois_cont, all_data))

h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
fit1 <- varComp(Low_final_k ~ asv , all_data, ~ pc1)
summary(fit1)
lrtest(fit, fit1)
levels(all_data$geocode) <- levels(all_data$geocode)[c(2, 3, 1)]
fit <- varComp(Low_final_k ~ geocode, all_data, varcov = list(com = com_sim, env = env_sim))
fit1 <- varComp(Low_final_k ~ geocode, all_data, varcov = list(env = env_sim))

varComp.test(fit, fit1)
lrtest(fit, fit1)
summary(fit)
h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))

fit <- varComp(Low_final_k ~ asv, all_data, ~  pc1 + env_pc1)
summary(fit)
fit1 <- varComp(Low_final_k ~ asv, all_data, ~  pc1)
summary(fit1)
lrtest(fit, fit1)

fit <- varComp(Low_final_k ~ asv, all_data, ~  pc1 + env_pc1)
summary(fit)
fit1 <- varComp(Low_final_k ~ asv, all_data, ~  env_pc1)
summary(fit1)
lrtest(fit, fit1)


with(all_data, cor(pc1, env_pc1))
with(all_data, cor(pc1, geo_pc1))
with(all_data, cor(env_pc1, geo_pc1))

model <- varComp(Low_final_k ~ asv + geocode + Wetland, data = all_data,
                 varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
summary(model)

model1 <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
summary(model1)
lrtest(model, model1)

h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))
all_data <- all_data[all_data$Wetland == 'Upland', ]


model <- varComp(Low_final_k ~ asv , data = dry_data, ~ pc1  + env_pc1 + geocode)
summary(model)

h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(Low_final_k ~ geocode, data = all_data,
                 varcov = list(com = com_sim, env = env_sim) )
summary(model)

model1 <- varComp(Low_final_k ~ asv, data = dry_data,
                 varcov = list(com_sim,  env_sim) )
lrtest(model, model1)

model <- varComp(Low_final_k ~ asv, data = all_data, ~ geo_pc1 + pc1 + 
pc2 + pc3 + env_pc1 + env_pc2)
summary(model)
