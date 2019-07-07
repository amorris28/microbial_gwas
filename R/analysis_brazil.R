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
#library(gap)
library(parallel)
#library(qvalue)
library(knitr)
library(kableExtra)
#library(eulerr)
#library(GGally)
library(gap)
library(lmtest)

num.cores <- detectCores() - 1
#+ r import_data

taxon_table <- data.table::fread('../output/brazil_taxon_table.csv', data.table = FALSE)
all_data <- data.table::fread('../output/brazil_cleaned_data.tsv', data.table = FALSE)
all_data <- all_data[!is.na(all_data$Long), ]
n_samp <- nrow(all_data)
n_asvs <- ncol(all_data[, grepl("^[asv]", colnames(all_data))])

## Remove ASVs only present in 1 sample
asvs <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
vars <- all_data[, !grepl("^[asv]", colnames(all_data))]
asvs <- asvs[, colSums(asvs) > 0]
all_data <- cbind(vars, asvs)

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
#' These data were collected from two locations in Brazil and include `r n_samp` samples.
#' Microbial community data include 16S rRNA gene sequences 
#' with `r n_asvs` ASVs from the DADA2 pipeline. The community matrix was transformed
#' using the `varianceStabilizingTransformation` from `DESeq2`.
#' This analysis includes the flux of methane, carbon dioxide, and nitrous oxide.
#' Due to strong right-skewness and the presence of negative fluxes, I transformed
#' CH4 by adding a constant `r abs(min(all_data$CH4)) + 1` and taking the naural log.
#' CO2 was also skewed, but had no negative values so I took the natural log.
#' N2O was roughly normally distributed and so received on transformation.
#' The spatial data are latitude and longitude. 
#' Environmental data include a number of soil physicochemical factors including
#' pH, soil moisture, macro- and micro-nutrients, structure, and texture. 
#'
#' # Overall Summary
#'
#' The traditional approach to microbial ecosystem function studies as I argue in the
#' paper is to correlate some aspect of microbial biodiversity such as abundance,
#' richness, diversity, or composition with the rate of an ecosystem function.
#' This largely fails as demonstrated by the Rocca and Graham papers. Here I
#' perform some of these analyses, notably excluding marker genes (but including 
#' PC regression, diversity,
#' richness, and various environmental factors). None of the geographic
#' or environmental variables have a correlation coefficient (r) greater
#' than about ~0.4. CH4 does not vary by region, but is significantly increased
#' in the pasture, particularly in FNV. CH4 is positively correlated with PC1
#' and CO2 is negatively correlated with PC1. CH4 is negatively correlated
#' with PC2. Although,
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

#+ trad_correlations, fig.width=10, fig.height=10

colnames(all_data[, 1:40])
#ggcorr(all_data[, c(2:3, 8:34, 36, 38)], nbreaks = 4, label = TRUE)

#+ other_correlations
all_data$log_CH4_c <- with(all_data, log(CH4 + abs(min(CH4)) + 1))
pca <- rda(as.matrix(all_data[, grepl("^[asv]", colnames(all_data))]))
screeplot(pca)

all_data$pc1 <- scores(pca, choices = 1, display = 'sites', scaling = 3)
all_data$pc2 <- scores(pca, choices = 2, display = 'sites', scaling = 3)
all_data$pc3 <- scores(pca, choices = 3, display = 'sites', scaling = 3)

ggplot(all_data, aes(x = pc1, y = CH4)) +
       geom_point() +
       stat_smooth(method = 'lm')
ggplot(all_data, aes(x = pc1, y = rank(CH4))) +
       geom_point() +
       stat_smooth(method = 'lm')
ggplot(all_data, aes(x = pc1, y = log(CH4 + abs(min(CH4)) + 1))) +
       geom_point() +
       stat_smooth(method = 'lm')
summary(lm(log(CH4 + abs(min(CH4)) + 1) ~ pc1, data = all_data))

ggplot(all_data, aes(x = pc1, y = log(CO2))) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(log(CO2) ~ pc1, data = all_data))

# Looking at this regression, I think cauchy regression would work
# better for N2O
ggplot(all_data, aes(x = pc1, y = N2O)) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(N2O ~ pc1, data = all_data))

ggplot(all_data, aes(x = pc2, y = log(CH4 + abs(min(CH4)) + 1))) +
       geom_point() +
       stat_smooth(method = 'lm')
summary(lm(log(CH4 + abs(min(CH4)) + 1) ~ pc2, data = all_data))

ggplot(all_data, aes(x = pc2, y = log(CO2))) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(log(CO2) ~ pc2, data = all_data))

# Looking at this regression, I think cauchy regression would work
# better for N2O
ggplot(all_data, aes(x = pc2, y = N2O)) +
       geom_point() +
       geom_smooth(method = 'lm')
summary(lm(N2O ~ pc2, data = all_data))

#my_asvs <- tibble(asv = names(scores(pca)$species[, 'PC1']),
#          score = scores(pca)$species[, 'PC1'])
#taxon_table  %>% 
#  left_join(my_asvs) %>% 
#  filter(abs(score) > 0.6) %>% 
#  kable(digits = 2) %>% 
#  kable_styling()
#
#my_asvs <- tibble(asv = names(scores(pca)$species[, 'PC2']),
#          score = scores(pca)$species[, 'PC2'])
#taxon_table  %>% 
#  left_join(my_asvs) %>% 
#  filter(abs(score) > 0.5) %>% 
#  kable(digits = 2) %>% 
#  kable_styling()
env_data <- select(all_data, pH:N_total)
is.na(env_data)
env_pca <- rda(env_data, na.action = na.exclude, scale = TRUE)
screeplot(env_pca)
env_pca
all_data$env_pc1 <- scores(env_pca, choices = 1, display = 'sites', scaling = 3)
all_data$env_pc2 <- scores(env_pca, choices = 2, display = 'sites', scaling = 3)
all_data$env_pc3 <- scores(env_pca, choices = 3, display = 'sites', scaling = 3)

ggplot(all_data, aes(x = env_pc1, y = log_CH4_c)) +
  geom_point()

ggplot(all_data, aes(x = env_pc2, y = log_CH4_c)) +
  geom_point()

ggplot(all_data, aes(x = env_pc3, y = log_CH4_c)) +
  geom_point()

ggplot(all_data, aes(x = Total_moisture, y = log_CH4_c)) +
  geom_point()
summary(lm(log_CH4_c ~ Total_moisture, data = all_data))

ggplot(all_data, aes(x = Region, y = log_CH4_c)) +
       geom_boxplot() +
       geom_jitter()
summary(lm(log_CH4_c ~ Region, data = all_data))
ggplot(all_data, aes(x = Land_type, y = log_CH4_c)) +
       geom_boxplot() +
       geom_jitter()
summary(lm(log_CH4_c ~ Land_type, data = all_data))
ggplot(all_data, aes(x = interaction(Land_type, Region), y = log_CH4_c)) +
       geom_boxplot() +
       geom_jitter()
summary(lm(log_CH4_c ~ interaction(Land_type, Region), data = all_data))

com_matrix <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
richness <- rowSums(com_matrix > 0)
shannon <- diversity(com_matrix, "shannon")
simpson <- diversity(com_matrix, "simpson")
ggplot(data = all_data, aes(x = richness, y = log_CH4_c)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(log_CH4_c ~ richness, data = all_data))
ggplot(data = all_data, aes(x = richness, y = log(CO2))) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(log(CO2) ~ richness, data = all_data))
ggplot(data = all_data, aes(x = richness, y = N2O)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(N2O ~ richness, data = all_data))

ggplot(data = all_data, aes(x = shannon, y = log_CH4_c)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(log_CH4_c ~ shannon, data = all_data))
ggplot(data = all_data, aes(x = shannon, y = log(CO2))) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(log(CO2) ~ shannon, data = all_data))
ggplot(data = all_data, aes(x = shannon, y = N2O)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(N2O ~ shannon, data = all_data))

ggplot(data = all_data, aes(x = simpson, y = log_CH4_c)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(log_CH4_c ~ simpson, data = all_data))
ggplot(data = all_data, aes(x = simpson, y = log(CO2))) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(log(CO2) ~ simpson, data = all_data))
ggplot(data = all_data, aes(x = simpson, y = N2O)) +
  geom_point() +
  geom_smooth(method = 'lm')
summary(lm(N2O ~ simpson, data = all_data))
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
com_sim <- 1 - as.matrix(vegdist(asvs[, -c(x)]))
com_sim_bin <- 1 - as.matrix(vegdist(asvs[, -c(x)], method = 'jaccard', binary = TRUE))

# Geographic similarity
geo_sim <- calc_sim_matrix(all_data[, c('Lat', 'Long')])

# environmental similarity'
env_sim <- calc_sim_matrix(select(all_data, pH:Total_moisture))

#' ## Model diagnostics for all data

# Diagnostics
print(mantel(com_sim, geo_sim))
print(mantel(env_sim, geo_sim))
print(mantel(com_sim, env_sim))

print(ggplot(all_data, aes(x = Lat, y = Long)) + 
  geom_point(aes(color = all_data$Land_type)) +
  labs(title = "Distribution of sites highlighted by land use"))
print(ggplot(all_data, aes(x = Lat, y = Long)) + 
  geom_point(aes(color = all_data$Region)) +
  labs(title = "Distribution of sites highlighted by region"))

print(ggplot(all_data, aes(log_CH4_c)) + 
  geom_histogram(bins = 20) + 
  labs(x = "log_CH4_c", y = "Number of Samples"))
print(ggplot(mapping = aes(asv)) + 
  geom_histogram(bins = 20) + 
  labs(x = "Taxon Abundance", y = "Number of Samples",
  title = paste0('Abundance of ', asv_tax)))

print(ggplot(mapping = aes(c(com_sim))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Community Similarity (Bray-Curtis)", y = "Number of Samples"))
print(ggplot(mapping = aes(c(geo_sim))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Spatial Similarity (Euclidean)", y = "Number of Samples"))
print(ggplot(mapping = aes(c(env_sim))) + 
  geom_histogram(binwidth = 0.05) + 
  labs(x = "Environment Similarity (Euclidean)", y = "Number of Samples"))

model <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
model1 <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(geo = geo_sim, env = env_sim) )
varComp.test(model, model1)
lrtest(model, model1)
summary(model)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(com = com_sim) )
model1 <- varComp(log_CH4_c ~ 1,
                   data = all_data)
varComp.test(model, model1)
lrtest(model, model1)
summary(model)
h2G(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(com = com_sim_bin) )
model1 <- varComp(log_CH4_c ~ 1,
                   data = all_data)
varComp.test(model, model1)
lrtest(model, model1)
summary(model)
h2G(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim_bin, env = env_sim) )
model1 <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(geo = geo_sim, env = env_sim) )
varComp.test(model, model1)
lrtest(model, model1)
summary(model)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
model1 <- varComp(log_CH4_c ~ interaction(Region, Land_type),
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
lrtest(model, model1)
varComp.test(model, model1)
summary(model1)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model_null <- varComp(log_CH4_c ~ interaction(Region, Land_type),
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim) )
model <- varComp(log_CH4_c ~ interaction(Region, Land_type),
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
lrtest(model_null, model)
summary(model)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))


model <- varComp(log_CH4_c ~  Region + Land_type,
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim, env = env_sim) )
model1 <- varComp(log_CH4_c ~  Region + Land_type,
                   data = all_data,
                   varcov = list( com = com_sim, env = env_sim) )
varComp.test(model, model1)

model <- varComp(log_CH4_c ~  Region + Land_type,
                   data = all_data,
                   varcov = list(com = com_sim) )
model1 <- varComp(log_CH4_c ~ 1,
                   data = all_data,
                   varcov = list( com = com_sim) )

varComp.test(model1, model)
lrtest(model1, model)

model <- varComp(log_CH4_c ~ Region + Land_type,
                   data = all_data,
                   varcov = list(com = com_sim) )
model1 <- varComp(log_CH4_c ~ Region + Land_type,
                   data = all_data)
varComp.test(model, model1)
lrtest(model, model1)
h2G(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(log_CH4_c ~ Region + Land_type,
                   data = all_data,
                   varcov = list(com = com_sim, env = env_sim) )
model1 <- varComp(log_CH4_c ~ Region + Land_type,
                   data = all_data, varcov = list(env = env_sim))
varComp.test(model, model1)
lrtest(model, model1)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

model <- varComp(log_CH4_c ~ Region,
                   data = all_data,
                   varcov = list(com = com_sim, env = env_sim) )
model1 <- varComp(log_CH4_c ~ Region,
                   data = all_data, varcov = list(env = env_sim))
varComp.test(model, model1)
lrtest(model, model1)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

#+ exploratory_stuff, eval=FALSE

model_geo <- varComp(log_CH4_c ~ asv + Region + Land_type,
                   data = all_data,
                   varcov = list(com = com_sim, env = env_sim) )
model_env <- varComp(log_CH4_c ~ asv + Region + Land_type,
                   data = all_data,
                   varcov = list(geo = geo_sim, com = com_sim) )
lrtest(model_full, model_land, model_reg, model_com, model_geo, model_env)
lrtest(model_full, model_land)
lrtest(model_full, model_reg)
krtest(model_full, model_com)
lrtest(model_full, model_geo)
lrtest(model_full, model_env)

fit <- varComp(log_CH4_c ~ 1, all_data, varcov = list(geo = geo_sim, env = env_sim, com = com_sim))
fit1 <- varComp(log_CH4_c ~ 1, all_data, varcov = list(geo = geo_sim, env = env_sim))

varComp.test(fit, fit1)
lrtest(fit, fit1)
summary(fit)
summary(fit1)
h2GE(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
h2G(c(fit$varComps, err = fit$sigma2), vcov(fit, what = 'varComp', drop = FALSE))
model <- varComp(log_CH4_c ~ asv,
                   data = all_data,
                   varcov = list(com = com_sim) )
summary(model)

h2G(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

h2G(c(model$varComps['com'], err = model$sigma2), vcov(model, what = 'varComp', drop = FALSE))

all_data$locale <- interaction(all_data$Region, all_data$Land_type)
unique(all_data$locale)
model <- varComp(log_CH4_c ~ asv + locale,
                   data = all_data,
                   varcov = list(com = com_sim, env = env_sim) )
summary(model)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'random', drop = FALSE))

model <- varComp(log_CH4_c ~ asv ,
                   data = all_data, ~ pc1 + env_pc2 + Land_type)
summary(model)
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'random', drop = FALSE))

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


#' # Analyze model output

#+ r analyze_model_output



my_model_output <- as.data.frame(data.table::fread('../talapas-output/var_comp_out_brazil.csv'))


taxon_table$taxon <- apply(taxon_table, 1, function(r) { r[ 
                      if( any(is.na(r))){ 
                          (which(is.na(r))[1]-1) 
                       }else{
                          (length(r)-0)}
                                      ] }
       )

#' ## Variance components

#+ r variance_components
my_model_output %>% 
  group_by(comps) %>% 
  summarize_at(vars(com, env, err), mean) %>% 
  kable(caption = "Variance components for all data",
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


#' ## Taxa identified with each covariate


#+ r plot_shared_taxa
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

#plot_asv_euler(my_model_output)

model <- varComp(log_CH4_c ~ asv + Region + Land_type, data = all_data,
                 varcov = list(com = com_sim, env = env_sim) )
summary(model)
vcov(model, what = 'varComp')
h2 <- h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp'))

var_comps_se <- data.frame(comp = c('com', 'env'), h2 = c(h2$h2G, h2$h2GE), 
           SE = c(sqrt(h2$Varh2G), sqrt(h2$Varh2GE)))
ggplot(var_comps_se, aes(x = comp, y = h2, ymax = h2 + SE, ymin = h2 - SE)) +
  geom_pointrange() 

sessionInfo()
