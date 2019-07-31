# Initialize packages ###########
library(tidyverse)
library(phyloseq)
library(corncob)
library(breakaway)
library(DivNet)
library(vegan)
library(ape)

# Settings ############
theme_set(theme_classic())

# Import data #################3
physeq <- readRDS('../output/physeq.rds')

# Calculate positive methane oxidation for easier model interpretation
sample_data(physeq)$pos_lowk <- -sample_data(physeq)$Low_final_k

# Subset methane oxidation experiment
physeq <- physeq %>% 
            subset_samples(Experiment == 'MO')
sample_data(physeq)$X_km <- sample_data(physeq)$X / 1000
sample_data(physeq)$Y_km <- sample_data(physeq)$Y / 1000
# Center and scale environmental data
sample_data(physeq)[, c('X', 'Y', 'Mois_cont', 'N_percent', 'C_percent', 'Bulk_dens')] <-
  scale(sample_data(physeq)[, c('X', 'Y', 'Mois_cont', 'N_percent', 'C_percent', 'Bulk_dens')])

# Cluster taxa at the family level
phyl_fam <- physeq %>% 
            tax_glom("Family")

ncol(otu_table(phyl_fam))
rowSums(otu_table(phyl_fam))

# Diversity ####################
# Compute alpha and beta diversity metrics
divnet_fam <-  divnet(phyl_fam, ncores = 4)

# PCoA of Bray-Curtis
pcoa <- pcoa(divnet_fam$`bray-curtis`)
shan <- summary(divnet_fam$shannon)$estimate

# Screeplot
plot(seq_len(length(pcoa$values[, 1])), pcoa$values[, 1])

sample_data(phyl_fam) %>% 
  ggplot(mapping = aes(x = pmoa_copy_num, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = expression(italic('pmoA') ~ 'copy number'), y = expression('Methane oxidation rate (-k)')) + 
  annotate(geom = 'text', x = min(sample_data(physeq)$pmoa_copy_num, na.rm = TRUE), 
            y = max(sample_data(physeq)$pos_lowk), label = 'A', size = 6)

#ggsave(file = '../figures/lowk_pmoa.pdf', width = 4, height = 4)

sample_data(phyl_fam) %>% 
  ggplot(mapping = aes(x = shan, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Shannon Diversity', y = expression('Methane oxidation rate (-k)')) + 
  annotate(geom = 'text', x = min(shan),
           y = max(sample_data(physeq)$pos_lowk), label = 'B', size = 6)

#ggsave(file = '../figures/lowk_pc1.pdf', width = 4, height = 4)

fit <- lm(Low_final_k ~ pmoa_copy_num, data = as(sample_data(phyl_fam), 'data.frame'))
summary(fit)

fit <- lm(Low_final_k ~ shan, data = as(sample_data(phyl_fam), 'data.frame'))
summary(fit)

# Pull out important axes an plot pairwise with axes proportional to var expl
sample_data(phyl_fam)$pc1 <- pcoa$vectors[, 1]
sample_data(phyl_fam)$pc2 <- pcoa$vectors[, 2]
sample_data(phyl_fam)$pc3 <- pcoa$vectors[, 3]

summary(loess(Y ~ X, as(sample_data(phyl_fam), 'data.frame'))) 
Z <- predict(loess(Y ~ X, as(sample_data(phyl_fam), 'data.frame'))) 
loess
ggplot(sample_data(phyl_fam), aes(X, Y)) +
  geom_point() +
  stat_smooth()

ggplot(sample_data(phyl_fam), aes(x = X, y = Y, color = Z, shape = geocode)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[2, 2]/pcoa$values[1, 2]) +
  scale_color_gradient(name = "Distance") +
  scale_shape(name = 'Site')

ggplot(sample_data(phyl_fam), aes(x = pc1, y = pc2, color = Z, shape = geocode)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[2, 2]/pcoa$values[1, 2]) +
  scale_color_gradient(name = "Distance") +
  scale_shape(name = 'Site')
  #+
 # annotate("text", x = 0.03 + sample_data(phyl_fam)$pc1, y = sample_data(phyl_fam)$pc2, label = sample_names(phyl_fam))

ggplot(sample_data(phyl_fam), aes(x = pc2, y = pc3, color = Z, shape = geocode)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[2, 2]) +
  scale_color_gradient(name = "Distance") +
  scale_shape(name = 'Site')

ggplot(sample_data(phyl_fam), aes(x = pc1, y = pc3, color = Z, shape = geocode)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[1, 2]) +
  scale_color_gradient(name = "Distance") +
  scale_shape(name = 'Site')

# Environmental and spatial distance ################
# Environment
env_dis <- dist(as.matrix(sample_data(phyl_fam)[, c('Mois_cont', 'Bulk_dens', 'N_percent', 'C_percent')]))
env_sim <- 1/(1+dist(as.matrix(sample_data(phyl_fam)[, c('Mois_cont', 'Bulk_dens', 'N_percent', 'C_percent')])))
# Spatial
geo_dis <- dist(as.matrix(sample_data(phyl_fam)[, c('X', 'Y')]))
geo_sim <- calc_sim_matrix(sample_data(phyl_fam)[, c('X', 'Y')])

# Test multi-collinearity of community, environment, and geography ###############
com_dis <- divnet_fam$`bray-curtis`
com_dis <- vegdist(otu_table(physeq))
mantel(com_dis, env_dis)
mantel(com_dis, geo_dis)
mantel(env_dis, geo_dis)
mantel.partial(com_div, geo_dis, env_dis)

# Test differential abundance of taxa associated with changes in methane oxidation rate
# Rerun the test with different sets of covariates to get at different aspects of
# community structure
fit <- differentialTest(formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X * Y + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X * Y + pc1 + pc2 + pc3,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent + X * Y + pc1 + pc2 + pc3,
                        phi.formula_null = ~ Bulk_dens + Mois_cont + N_percent + C_percent + X * Y + pc1 + pc2 + pc3,
                        data = phyl_fam,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit)
length(fit$significant_taxa)

fit1 <- differentialTest(formula = ~ Low_final_k + Bulk_dens + Mois_cont + N_percent + C_percent + X * Y,
                        phi.formula = ~ Low_final_k + Bulk_dens + Mois_cont + N_percent + C_percent + X * Y,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent + X * Y,
                        phi.formula_null = ~ Bulk_dens + Mois_cont + N_percent + C_percent + X * Y,
                        data = phyl_fam,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit1)

length(fit1$significant_taxa)

fit2 <- differentialTest(formula = ~ Low_final_k + X * Y + pc1 + pc2,
                        phi.formula = ~ Low_final_k  + X * Y + pc1 + pc2,
                        formula_null = ~  X * Y + pc1 + pc2,
                        phi.formula_null = ~ X * Y + pc1 + pc2,
                        data = phyl_fam,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit2)


length(fit2$significant_taxa)

fit3 <- differentialTest(formula = ~ Low_final_k + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        phi.formula = ~ Low_final_k + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent + pc1 + pc2 + pc3,
                        phi.formula_null = ~ Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        data = phyl_fam,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit3)

length(fit3$significant_taxa)


fit4 <- differentialTest(formula = ~ Low_final_k,
                        phi.formula = ~ Low_final_k,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = phyl_fam,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit4)

length(fit4$significant_taxa)

