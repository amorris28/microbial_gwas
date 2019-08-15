library(tidyverse)
library(broom)
library(varComp)
library(vegan)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(phyloseq)

registerDoParallel(cores = 4)

#import_data
physeq <- readRDS('../output/physeq_vst.rds')

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


# Fit models

calc_varcomp_per <- function(model) {
  as_tibble(as.list(c(p.value = anova(model)$`Pr(>F)`[2], F.value = anova(model)$`F value`[2],
c(model$varComps, err = model$sigma2)/sum(c(model$varComps, model$sigma2)))))
}

fit_varComps <- function(x, y, all_data, func = "Low_final_k") {
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
system.time(results <- foreach(i = seq_len(ncol(asvs)), .combine = 'rbind', 
                               .packages = c('vegan', 'varComp', 'tidyverse')) %dopar% { 
            suppressMessages(suppressWarnings(fit_varComps(x = i, all_data = all_data)))
                 })
fwrite(results, "../output/var_comp_abund.csv")
