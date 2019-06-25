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
library(eulerr)
num.cores <- detectCores()
#+ r import_data

all_data <- as.data.frame(data.table::fread('../output/gab_all_vst_troph.csv'))
n_samp <- nrow(all_data)
n_asvs <- ncol(all_data[, grepl("^[asv]", colnames(all_data))])

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
#' These data were collected from Gabon, Africa and include `r n_samp` samples
#' and `r n_asvs` ASVs from the DADA2 pipeline. The community matrix was transformed
#' using the `varianceStabilizingTransformation` from `DESeq2`.
#' This analysis only includes high affinity
#' methane oxidation measurements (Low Final K). The spatial data were converted
#' from latitude/longitude to distance in meters by converting to UTM. 
#' Environmental data include bulk density, soil moisture, percent nitrogen, and 
#' percent carbon.
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

## Remove ASVs only present in 1 sample
asvs <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
vars <- all_data[, !grepl("^[asv]", colnames(all_data))]
asvs <- asvs[, colSums(asvs) > 0]
all_data <- cbind(vars, asvs)

# Community similarity

x <- 3
asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
com_sim <- 1 - as.matrix(vegdist(asvs[, -c(x)]))

# Geographic similarity
geo_sim <- calc_sim_matrix(all_data[, c('X', 'Y')])

# environmental similarity'
env_sim <- calc_sim_matrix(all_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])

# 
all_data[all_data$Y > 30000, 'geocode'] <- 'A'
all_data[all_data$Y < 30000 & all_data$Y > 10000, 'geocode'] <- 'B'
all_data[all_data$Y < 10000, 'geocode'] <- 'C'
all_data$geocode <- factor(all_data$geocode)


#' ## Model diagnostics for all data

#+ r varcomp_diagnostics

varcomp_diag <- function(data) {
print(mantel(com_sim, geo_sim))

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


#' ## Model diagnostics for Wetland data

#+ r varcomp_diagnostics_wet

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

#' ## Model diagnostics for Upland data

#+ r varcomp_diagnostics_dry

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


varcomp_diag(dry_data)

#' # Analyze model output

#+ r analyze_model_output



taxon_table <- as.data.frame(data.table::fread('../output/taxon_table.csv'))
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

plot_asv_euler(my_model_output)

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

plot_asv_euler(wet_model_output)

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


plot_asv_euler(dry_model_output)

sessionInfo()
