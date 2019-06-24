#' ---
#' title: "Royal Society paper"
#' subtitle: "Linking microbial communities to ecosystem functions: what we can learn from genotype-phenotype mapping in organisms"
#' author: "Andrew Morris"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---
#' 
#' # Overview
#' 
#' This is the primary analysis script for the analysis associated with the Royal Society
#' paper. 

#+ r global_options, include=FALSE
knitr::opts_chunk$set(fig.width=8, fig.height=6, fig.path='Figs/',
                      echo=FALSE, warning=FALSE, message=FALSE)

#+ r load_packages, include=FALSE

library(tidyverse)
library(lubridate)
library(broom)
#library(broman)
#library(morris)
library(varComp)
library(vegan)
library(gap)
library(parallel)
#library(qvalue)
num.cores <- detectCores()
#+ r import_data

all_data <- as.data.frame(data.table::fread('../output/gab_all_vst_troph.csv'))

#' # Variance Component Analysis

#+ r sim_mat

calc_sim_matrix <- function(x) {
1/(1+as.matrix(dist(as.matrix(scale(x)))))
}

## Remove ASVs only present in 1 sample
asvs <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
vars <- all_data[, !grepl("^[asv]", colnames(all_data))]
asvs <- asvs[, colSums(asvs) > 0]
all_data <- cbind(vars, asvs)

## Community similarity
#asv <- 3
#taxon <- asvs[, asv]
#colnames(asvs)[asv]
#com_sim <- 1 - as.matrix(vegdist(asvs[, -c(asv)]))
#
## Geographic similarity
#geo_sim <- calc_sim_matrix(all_data[, c('X', 'Y')])
#
## environmental similarity'
#env_sim <- calc_sim_matrix(all_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])

# 
all_data[all_data$Y > 30000, 'geocode'] <- 'A'
all_data[all_data$Y < 30000 & all_data$Y > 10000, 'geocode'] <- 'B'
all_data[all_data$Y < 10000, 'geocode'] <- 'C'
all_data$geocode <- factor(all_data$geocode)

#+ r varcomp_diagnostics
#
#mantel(com_sim, geo_sim)
#qplot(x = xy$X, y = xy$Y) + geom_point(aes(color = all_data$Wetland))
#qplot(x = xy$X, y = xy$Y) + geom_point(aes(color = all_data$geocode))
#
#ggplot(all_data, aes(Low_final_k)) + 
#  geom_histogram(bins = 20) + 
#  labs(x = "Low Final K", y = "Number of Samples")
#ggplot(all_data, aes(taxon)) + 
#  geom_histogram(bins = 20) + 
#  labs(x = "Taxon Abundance", y = "Number of Samples")
#
#
#ggplot(mapping = aes(c(com_sim))) + 
#  geom_histogram(binwidth = 0.05) + 
#  labs(x = "Community Similarity (Bray-Curtis)", y = "Number of Samples")
#ggplot(mapping = aes(c(geo_sim))) + 
#  geom_histogram(binwidth = 0.05) + 
#  labs(x = "Spatial Similarity (Euclidean)", y = "Number of Samples")
#ggplot(mapping = aes(c(env_sim))) + 
#  geom_histogram(binwidth = 0.05) + 
#  labs(x = "Environment Similarity (Euclidean)", y = "Number of Samples")

#'
#'Community and geography are correlated and sites cluster around three locations
#'Also, geography uses up more of the similarity matrix space (see histograms)
#'Therefore, I coded geography as three factors (one for each cluster of sites)
#'Also, wetland/upland covaries with those three sites, but moisture content is
#'included in the environmental similarity matrix and should capture that variation

#+ r fit_var_comp


#model <- varComp(Low_final_k ~ taxon, data = all_data,
#                 varcov = list(com = com_sim, env = env_sim, geo = geo_sim) )
#summary(model)
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



calc_varcomp_per <- function(model) {
  as_tibble(as.list(c(p.value = anova(model)$`Pr(>F)`[2],
c(model$varComps, err = model$sigma2)/sum(c(model$varComps, model$sigma2)))))
}

fit_varComps <- function(x, all_data) {
  print(x/ncol(asvs))
asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
parts <- NULL
model <- varComp(Low_final_k ~ asv + geocode,
                 data = all_data,
                 varcov = list(com = com_sim, env = env_sim) )
parts <- bind_cols(tibble(asv = asv_tax, comps = 'gce', geo = T), calc_varcomp_per(model))

model <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim, env = env_sim) )
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'ce', geo = F), calc_varcomp_per(model)))

model <- varComp(Low_final_k ~ asv + geocode, data = all_data,
                 varcov = list(com = com_sim))

parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'gc', geo = T), calc_varcomp_per(model)))
model <- varComp(Low_final_k ~ asv + geocode, data = all_data,
                 varcov = list(env = env_sim))

parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'ge', geo = T), calc_varcomp_per(model)))
model <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(com = com_sim))

parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'c', geo = F), calc_varcomp_per(model)))
model <- varComp(Low_final_k ~ asv, data = all_data,
                 varcov = list(env = env_sim))

parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'e', geo = F), calc_varcomp_per(model)))
model <- varComp(Low_final_k ~ asv + geocode, data = all_data)

parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'g', geo = T), calc_varcomp_per(model)))
model <- varComp(Low_final_k ~ asv, data = all_data)

parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, comps = 'n', geo = F), calc_varcomp_per(model)))
parts
}


#system.time(my_model_output <- suppressMessages(suppressWarnings(do.call(rbind, mclapply(seq(ncol(asvs)), fit_varComps, mc.cores = num.cores)))))
#write.csv(my_model_output, "../output/var_comp_out.csv")

#+ fit_var_comp_wet

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


system.time(wet_model_output <- suppressMessages(suppressWarnings(do.call(rbind, mclapply(seq(100), fit_varComps, all_data = wet_data, mc.cores = 4)))))
system.time(wet_model_output <- suppressMessages(suppressWarnings(do.call(rbind, mclapply(seq(100), fit_varComps, all_data = wet_data, mc.cores = num.cores)))))
#system.time(wet_model_output <- suppressMessages(suppressWarnings(do.call(rbind, mclapply(seq(ncol(asvs)), fit_varComps, all_data = wet_data, mc.cores = num.cores)))))
#write.csv(wet_model_output, "../output/var_comp_out_wet.csv")

#+ fit_var_comp_dry

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

#system.time(dry_model_output <- suppressMessages(suppressWarnings(do.call(rbind, mclapply(seq(ncol(asvs)), fit_varComps, all_data = dry_data, mc.cores = num.cores)))))
#write.csv(dry_model_output, "../output/var_comp_out_dry.csv")

#+ analyze_model_output
#
#my_model_output <- as.data.frame(data.table::fread('../output/var_comp_out.csv'))
#
#head(my_model_output)
#colnames(my_model_output)
#my_model_output %>% 
#  group_by(comps) %>% 
#  summarize_at(vars(com, env, err), mean)
#
#qvalue(my_model_output$p.value)
