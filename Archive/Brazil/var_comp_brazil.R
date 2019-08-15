library(tidyverse)
library(broom)
library(varComp)
library(vegan)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(microbenchmark)

registerDoParallel(cores = 3)
num.cores <- detectCores() - 1

all_data <- as.data.frame(fread('../output/brazil_cleaned_data.tsv'))
all_data <- all_data[!is.na(all_data$Long), ]

# Variance Component Analysis

# sim_mat
calc_sim_matrix <- function(x) {
1/(1+as.matrix(dist(as.matrix(scale(x)))))
}

## Remove ASVs only present in 1 sample
asvs <- as.matrix(all_data[, grepl("^[asv]", colnames(all_data))])
vars <- all_data[, !grepl("^[asv]", colnames(all_data))]
asvs <- asvs[, colSums(asvs) > 0]
all_data <- cbind(vars, asvs)

## Functions
#hist(all_data[, 'CH4'])
#hist(log(all_data$CH4 + abs(min(all_data$CH4)) + 1))
#hist(all_data[, 'CO2'])
#hist(log(all_data[, 'CO2']))
#hist(all_data[, 'N2O'])
all_data$log_CH4_c <- log(all_data$CH4 + abs(min(all_data$CH4)) + 1)
all_data$log_CO2 <- log(all_data$CO2)
functions <- all_data[, c('log_CH4_c', 'log_CO2', 'N2O')]
# Community similarity
asv <- 3
taxon <- asvs[, asv]
asv_tax <- colnames(asvs)[asv]
com_sim <- 1 - as.matrix(vegdist(asvs[, -c(asv)]))

# Geographic similarity
geo_sim <- calc_sim_matrix(all_data[, c('Lat', 'Long')])

# environmental similarity'
env_sim <- calc_sim_matrix(select(all_data, pH:Total_moisture))

#plot(all_data$Lat, all_data$Long)


## Permute response
#permute_no_replace_for <- function(perm_vec, nperm = 1000) {
#samples <- matrix(ncol = nperm, nrow = length(perm_vec))
#samples <- cbind(samples, truth = perm_vec)
#for (i in seq_len(nperm)) {
#  sample_i <- sample(perm_vec)
#  while (any(apply(samples, 2, function(x) {all(sample_i == x)}), na.rm = T)) {
#    sample_i <- sample(perm_vec)
#  }
#  samples[, i] <- sample_i
#  colnames(samples)[i] <- paste0('perm.', i)
#}
#return(samples)
#}
#all_perm_lowk <- permute_no_replace_for(all_data$Low_final_k, nperm = 1000)


## Diagnostics
#print(mantel(com_sim, geo_sim))
#
#print(ggplot(all_data, aes(x = Lat, y = Long)) + 
#  geom_point(aes(color = all_data$Land_type)) +
#  labs(title = "Distribution of sites highlighted by land use"))
#print(ggplot(all_data, aes(x = Lat, y = Long)) + 
#  geom_point(aes(color = all_data$Region)) +
#  labs(title = "Distribution of sites highlighted by region"))
#
#print(ggplot(all_data, aes(CH4)) + 
#  geom_histogram(bins = 20) + 
#  labs(x = "CH4", y = "Number of Samples"))
#print(ggplot(mapping = aes(taxon)) + 
#  geom_histogram(bins = 20) + 
#  labs(x = "Taxon Abundance", y = "Number of Samples",
#  title = paste0('Abundance of ', asv_tax)))
#
#print(ggplot(mapping = aes(c(com_sim))) + 
#  geom_histogram(binwidth = 0.05) + 
#  labs(x = "Community Similarity (Bray-Curtis)", y = "Number of Samples"))
#print(ggplot(mapping = aes(c(geo_sim))) + 
#  geom_histogram(binwidth = 0.05) + 
#  labs(x = "Spatial Similarity (Euclidean)", y = "Number of Samples"))
#print(ggplot(mapping = aes(c(env_sim))) + 
#  geom_histogram(binwidth = 0.05) + 
#  labs(x = "Environment Similarity (Euclidean)", y = "Number of Samples"))
#
#model <- varComp(log10(CH4 + abs(min(CH4)) + 1) ~ taxon + Region + Land_type,
#                   data = all_data,
#                   varcov = list(com = com_sim, env = env_sim) )
#summary(model)
#
#ggplot(all_data, aes(y = log10(CH4 + abs(min(CH4)) + 1), x = interaction(Region, Land_type))) + 
#  geom_boxplot() + 
#  geom_jitter()
#
# Fit models

calc_varcomp_per <- function(model) {
  as_tibble(as.list(c(p.value = anova(model)$`Pr(>F)`[2], F.value = anova(model)$`F value`[2],
c(model$varComps, err = model$sigma2)/sum(c(model$varComps, model$sigma2)))))
}

fit_varComps <- function(x, y, all_data) {
  asv <- asvs[, x]
  asv_tax <- colnames(asvs)[x]
  com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
  func_rate <- functions[, y]
  which_function <- colnames(functions)[y]
  parts <- NULL
  model <- varComp(func_rate ~ asv + Land_type + Region,
                   data = all_data,
                   varcov = list(com = com_sim, env = env_sim) )
  parts <- bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'gce'), calc_varcomp_per(model))
  
  model <- varComp(func_rate ~ asv, data = all_data,
                   varcov = list(com = com_sim, env = env_sim) )
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'ce'), calc_varcomp_per(model)))
  
  model <- varComp(func_rate ~ asv + Land_type + Region, data = all_data,
                   varcov = list(com = com_sim))
  
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'gc'), calc_varcomp_per(model)))
  model <- varComp(func_rate ~ asv + Land_type + Region, data = all_data,
                   varcov = list(env = env_sim))
  
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'ge'), calc_varcomp_per(model)))
  model <- varComp(func_rate ~ asv, data = all_data,
                   varcov = list(com = com_sim))
  
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'c'), calc_varcomp_per(model)))
  model <- varComp(func_rate ~ asv, data = all_data,
                   varcov = list(env = env_sim))
  
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'e'), calc_varcomp_per(model)))
  model <- varComp(func_rate ~ asv + Land_type + Region, data = all_data)
  
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'g'), calc_varcomp_per(model)))
  model <- varComp(func_rate ~ asv, data = all_data)
  
  parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, func = which_function, comps = 'n'), calc_varcomp_per(model)))
  parts
}

#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind', 
#                               .packages = c('vegan', 'varComp', 'tidyverse')) %:% 
#  foreach(j = seq_len(ncol(functions)), .combine = 'rbind') %dopar% {
#            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = all_data)))
#                 })
#fwrite(results, "../output/var_comp_out_brazil.csv")
microbenchmark(
system.time(results <- foreach(i = seq_len(10), .combine = 'rbind', 
                               .packages = c('vegan', 'varComp', 'tidyverse')) %:% 
  foreach(j = seq_len(1), .combine = 'rbind') %dopar% {
  print(paste0(round(i/ncol(asvs)*100, 2), ' % of ASVs and ', j, ' out of 3 functions'))
            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = all_data)))
                 }),
 system.time(my_model_output <- suppressMessages(suppressWarnings(do.call(rbind, mclapply(seq_len(10), fit_varComps, y = 1, all_data = all_data, mc.cores = num.cores))))),
 times = 10)
#fwrite(results, "../output/var_comp_out_brazil.csv")
 
##+ fit_var_comp_wet
#
#wet_data <- filter(all_data, Wetland == 'Wetland')
#
### Remove ASVs only present in 1 sample
#asvs <- as.matrix(wet_data[, grepl("^[asv]", colnames(wet_data))])
#vars <- wet_data[, !grepl("^[asv]", colnames(wet_data))]
#asvs <- asvs[, colSums(asvs) > 0]
#wet_data <- cbind(vars, asvs)
#
## Geographic similarity
##geo_sim <- calc_sim_matrix(wet_data[, c('X', 'Y')])
#
## environmental similarity'
#env_sim <- calc_sim_matrix(wet_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])
#
##x <- 6
##asv <- asvs[, x]
##asv_tax <- colnames(asvs)[x]
##com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
##
##model <- varComp(Low_final_k ~ asv, data = wet_data,
##                 varcov = list(com = com_sim, env = env_sim, geo = geo_sim) )
##summary(model)
#
#perm_lowk <- all_perm_lowk[all_data$Wetland == 'Wetland', ]
#
#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind') %:% 
#  foreach(j = seq_len(ncol(perm_lowk)), .combine = 'rbind') %dopar% {
#            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = wet_data)))
#                 })
#fwrite(results, "../output/var_comp_out_wet.csv")
#
##+ fit_var_comp_dry
#
#dry_data <- filter(all_data, Wetland == 'Upland')
#
### Remove ASVs only present in 1 sample
#asvs <- as.matrix(dry_data[, grepl("^[asv]", colnames(dry_data))])
#vars <- dry_data[, !grepl("^[asv]", colnames(dry_data))]
#asvs <- asvs[, colSums(asvs) > 0]
#dry_data <- cbind(vars, asvs)
#
## Geographic similarity
##geo_sim <- calc_sim_matrix(dry_data[, c('X', 'Y')])
#
## environmental similarity'
#env_sim <- calc_sim_matrix(dry_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])
#
##x <- 1
##asv <- asvs[, x]
##asv_tax <- colnames(asvs)[x]
##com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
##
##model <- varComp(Low_final_k ~ asv, data = dry_data,
##                 varcov = list(com = com_sim, env = env_sim, geo = geo_sim) )
##summary(model)
#
#perm_lowk <- all_perm_lowk[all_data$Wetland == 'Upland', ]
#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind') %:% 
#  foreach(j = seq_len(ncol(perm_lowk)), .combine = 'rbind') %dopar% {
#            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = dry_data)))
#                 })
#fwrite(results, "../output/var_comp_out_dry.csv")
