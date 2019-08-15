library(tidyverse)
library(broom)
library(varComp)
library(vegan)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(phyloseq)
#library(doSNOW)
#cl <- makeCluster(28, outfile = "")
#registerDoParallel(cl)
#registerDoSNOW(cl)
#import_data
registerDoParallel(cores = 4)
#taxon_table <- fread('../output/taxon_table.csv', data.table = FALSE, na.strings = "")
all_data <- as.data.frame(fread('../output/gab_all_vst_troph.csv'))
physeq <- readRDS('../output/physeq.rds')

physeq <- physeq %>% 
            subset_samples(Experiment == 'MO') %>% 
            tax_glom("Family")
# Variance Component Analysis

# sim_mat

calc_sim_matrix <- function(x) {
1/(1+as.matrix(dist(as.matrix(scale(x)))))
}

asvs <- as.matrix(otu_table(physeq))
# Total dataset n asvs
dim(asvs)
asvs <- asvs[, colSums(asvs) > 0]
# n asvs present in at least one sample
dim(asvs)


#### Cluster ASVs at the genus level
## n Unique ASVS
#length(unique(taxon_table$asv))
## n ASVs with no Genus assignment
#sum(is.na(taxon_table$Genus))
#dim(asvs)
#asvs <- cbind(Sample = all_data$Sample, all_data[, grepl("^[asv]", colnames(all_data))])
#asvs <- 
#  asvs %>% 
#  gather(asv, abundance, 2:ncol(asvs)) %>% 
#  left_join(taxon_table) %>% 
#  select(Sample, Genus, abundance) %>% 
#  group_by(Sample, Genus) %>% 
#  summarize(abundance = sum(abundance)) %>% 
#  spread(Genus, abundance)
#dim(asvs)
#asvs[1:5, 1:5]
#asvs <- as.matrix(asvs[, -1])
#dim(asvs)
#asvs <- asvs[, colSums(asvs) > 0]
#
###

#PA <- asvs
#PA[PA<2] = 0
#PA[PA>1] = 1
#PA[1:10, 1:10]
#sum(colSums(PA) == 44)
#sum(colSums(PA) == 0)
#sum(apply(PA,2,sum) == 44)
#sum(apply(PA,2,sum) == 0)
#PA <- PA[, !colSums(PA) %in% c(0, 44)]
#colnames(PA)
#asvs <- PA

#dim(asvs)
#sum(asvs > 0)
#sum(asvs > 3)
#abund <- asvs
#abund_ranks <- apply(asvs, 1, rank)
#t(abund[1:10, 1:10])
#abund_ranks[1:10, 1:10]
#abund_ranks <-abund_ranks - 500
#abund_ranks[1:10, 1:10]
#abund_ranks[abund_ranks < 1] <- 1
#abund_ranks[1:5,1:8]
#summary(abund_ranks)
#abund_ranks
#
#dim(asvs)
# Recombine
all_data <- as(sample_data(physeq), 'data.frame')
dim(all_data)



# Community similarity
#asv <- 3
#taxon <- asvs[, asv]
#colnames(asvs)[asv]
#com_sim <- 1 - as.matrix(vegdist(asvs[, -c(asv)]))

# Geographic similarity
geo_sim <- calc_sim_matrix(all_data[, c('X', 'Y')])

# environmental similarity'
env_sim <- calc_sim_matrix(all_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])

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

# Fit models

calc_varcomp_per <- function(model) {
  as_tibble(as.list(c(p.value = anova(model)$`Pr(>F)`[2], F.value = anova(model)$`F value`[2],
c(model$varComps, err = model$sigma2)/sum(c(model$varComps, model$sigma2)))))
}

fit_varComps <- function(x, y, all_data, func = "Low_final_k") {
  #print(paste0(round(x/ncol(asvs)*100, 2), ' % of ASVs and ', round(y/ncol(perm_lowk)*100, 2), ' % of permutations'))
  print(paste0(round(x/ncol(asvs)*100, 2), ' % of ASVs'))
asv <- asvs[, x]
asv_tax <- colnames(asvs)[x]
func <- all_data[, func]
#func <- perm_lowk[, y]
n_perm <- 1
#n_perm <- colnames(perm_lowk)[y]
com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
parts <- NULL

geo <- all_data$geocode
com <- cholRoot(com_sim)
env <- cholRoot(env_sim)

model <- varComp(func ~ asv ,
                 data = all_data,
                 ~ geo + com + env)
parts <- bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'gce'), calc_varcomp_per(model))

model <- varComp(func ~ asv, data = all_data,
                 ~ com + env)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'ce'), calc_varcomp_per(model)))

model <- varComp(func ~ asv, data = all_data,
                 ~ geo + com)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'gc'), calc_varcomp_per(model)))

model <- varComp(func ~ asv, data = all_data,
                 ~ geo + env)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'ge'), calc_varcomp_per(model)))

model <- varComp(func ~ asv, data = all_data,
                 ~ com)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'c'), calc_varcomp_per(model)))

model <- varComp(func ~ asv, data = all_data,
                 ~ env)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'e'), calc_varcomp_per(model)))

model <- varComp(func ~ asv, data = all_data,
                 ~ geo)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'g'), calc_varcomp_per(model)))

model <- varComp(func ~ asv, data = all_data)
parts <- bind_rows(parts, bind_cols(tibble(asv = asv_tax, perm = n_perm, comps = 'n'), calc_varcomp_per(model)))
parts
}
#perm_lowk <- all_perm_lowk

#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind', 
#                               .packages = c('vegan', 'varComp', 'tidyverse')) %:% 
#  foreach(j = seq_len(ncol(perm_lowk)), .combine = 'rbind') %dopar% {
#            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = all_data)))
#                 })
system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind', 
                               .packages = c('vegan', 'varComp', 'tidyverse')) %dopar% { 
            suppressMessages(suppressWarnings(fit_varComps(x = i, all_data = all_data)))
                 })
fwrite(results, "../output/var_comp_fam.csv")
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
##perm_lowk <- all_perm_lowk[all_data$Wetland == 'Wetland', ]
#
##system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind') %:% 
##  foreach(j = seq_len(ncol(perm_lowk)), .combine = 'rbind') %dopar% {
##            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = wet_data)))
##                 })
#
#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind') %dopar% {
#            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = wet_data, func = 'Vmax')))
#                 })
#fwrite(results, "../output/var_comp_out_wet.csv")

#+ fit_var_comp_dry
#
#dry_data <- filter(all_data, Wetland == 'Upland')
#
### Remove ASVs only present in 1 sample
#asvs <- as.matrix(dry_data[, grepl("^[asv]", colnames(dry_data))])
#vars <- dry_data[, !grepl("^[asv]", colnames(dry_data))]
#asvs <- asvs[, colSums(asvs) > 0]
##asvs <- decostand(asvs, method = 'pa')
##asvs <- asvs[,-4126]
#dry_data <- cbind(vars, asvs)
#
## Geographic similarity
##geo_sim <- calc_sim_matrix(dry_data[, c('X', 'Y')])
#
## environmental similarity'
#env_sim <- calc_sim_matrix(dry_data[, c('Bulk_dens', 'Mois_cont', 'N_percent', 'C_percent')])
#
#x <- 1
#asv <- asvs[, x]
#asv_tax <- colnames(asvs)[x]
#com_sim <- 1 - as.matrix(vegdist(asvs[, -x]))
#
#model <- varComp(Low_final_k ~ asv, data = dry_data,
#                 varcov = list(com = com_sim, env = env_sim, geo = geo_sim) )
#summary(model)

#perm_lowk <- all_perm_lowk[all_data$Wetland == 'Upland', ]

#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind') %:% 
#  foreach(j = seq_len(ncol(perm_lowk)), .combine = 'rbind') %dopar% {
#            suppressMessages(suppressWarnings(fit_varComps(x = i, y = j, all_data = dry_data)))
#                 })
#system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind') %dopar% {  
#            suppressMessages(suppressWarnings(fit_varComps(x = i, all_data = dry_data)))
#                 })
#
#fwrite(results, "../output/var_comp_out_dry.csv")
