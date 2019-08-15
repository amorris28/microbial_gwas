# Initialize packages ###########
library(tidyverse)
library(phyloseq)
library(corncob)
#library(breakaway)
library(DivNet)
library(vegan)
library(ape)
library(ade4)
library(broom)
library(knitr)
library(kableExtra)
library(DESeq2)

# Settings ############
theme_set(theme_classic())

# Import data #################3
physeq <- readRDS('../output/physeq_raw.rds')

# Calculate positive methane oxidation for easier model interpretation
sample_data(physeq)$pos_lowk <- -sample_data(physeq)$Low_final_k

# Subset methane oxidation experiment
physeq <- physeq %>% subset_samples(Experiment == 'MO')
sample_data(physeq)$X_km <- sample_data(physeq)$X / 1000
sample_data(physeq)$Y_km <- sample_data(physeq)$Y / 1000

# Center and scale environmental data
#sample_data(physeq)[, c('X', 'Y', 'Mois_cont', 'N_percent', 'C_percent', 'Bulk_dens')] <-
#  scale(sample_data(physeq)[, c('X', 'Y', 'Mois_cont', 'N_percent', 'C_percent', 'Bulk_dens')])
##
# Cluster taxa at the family level
# phyl_gen <- physeq %>% tax_glom("Genus")
phyl <- physeq %>% tax_glom("Phylum")

# Diversity ####################
# Compute alpha and beta diversity metrics
#bray_fam <-  divnet(phyl, ncores = 3)
divnet_phyl <-  divnet(phyl, ncores = 3)
divnet_phyl$shannon

#saveRDS(bray_fam, '../output/divnet_fam.rds')

divnet <- readRDS('../output/divnet_fam.rds')

# PCoA of Bray-Curtis
pcoa <- pcoa(divnet$`bray-curtis`)
shan <- summary(divnet$shannon)$estimate

# Screeplot
#ggplot(mapping = aes(seq_len(length(pcoa$values[, 1])), pcoa$values[, 1])) +
#       geom_col() +
#       labs(y = "Eigenvalues", x = "Principal Coordinate")

sample_data(phyl)$pc1 <- pcoa$vectors[, 1]
sample_data(phyl)$pc2 <- pcoa$vectors[, 2]
sample_data(phyl)$pc3 <- pcoa$vectors[, 3]

#
sample_data(physeq) %>% 
  ggplot(aes(x = pmoa_copy_num, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = expression(italic('pmoA') ~ 'copy number'), y = expression('Methane oxidation rate (-k)')) + 
  annotate(geom = 'text', x = min(sample_data(physeq)$pmoa_copy_num, na.rm = TRUE), 
            y = max(sample_data(physeq)$pos_lowk), label = 'A', size = 6)

ggsave(file = '../figures/lowk_pmoa.pdf', width = 4, height = 4)

sample_data(physeq) %>% 
  ggplot(mapping = aes(x = shan, y = pos_lowk)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'Shannon Diversity', y = expression('Methane oxidation rate (-k)')) + 
  annotate(geom = 'text', x = min(shan),
           y = max(sample_data(physeq)$pos_lowk), label = 'B', size = 6)

ggsave(file = '../figures/lowk_shan.pdf', width = 4, height = 4)

pmoa_fit <- lm(Low_final_k ~ pmoa_copy_num, data = as(sample_data(physeq), 'data.frame'))
summary(pmoa_fit)
glance(pmoa_fit)

shan_fit <- lm(Low_final_k ~ shan, data = as(sample_data(physeq), 'data.frame'))
summary(shan_fit)
glance(shan_fit)

rbind(tidy(pmoa_fit)[2, ], tidy(shan_fit)[2, ])

# Pull out important axes an plot pairwise with axes proportional to var expl

Z <- predict(loess(Y_km ~ X_km, as(sample_data(physeq), 'data.frame'))) 

ggplot(sample_data(phyl), aes(X_km, Y_km)) +
  geom_point() +
  stat_smooth() +
  coord_fixed()

ggplot(sample_data(phyl), aes(x = X_km, y = Y_km, color = Mois_cont)) +
  geom_point(size = 3, alpha = 0.1) +
  scale_color_gradient(name = "Moisture") +
  scale_shape(name = 'Site') +
  coord_fixed()

ggplot(sample_data(phyl), aes(x = pc1, y = pc2, color = Mois_cont)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[2, 2]/pcoa$values[1, 2]) +
  scale_color_gradient(name = "Moisture")
ggplot(sample_data(phyl), aes(Bulk_dens, pos_lowk)) + 
  geom_point() +
  stat_smooth(method = 'lm')

#+
# annotate("text", x = 0.03 + sample_data(phyl_fam)$pc1, y = sample_data(phyl_fam)$pc2, label = sample_names(phyl_fam))

ggplot(sample_data(phyl), aes(x = pc2, y = pc3, color = Mois_cont)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[2, 2]) +
  scale_color_gradient(name = "Moisture") +
  scale_shape(name = 'Site')

ggplot(sample_data(phyl), aes(x = pc1, y = pc3, color = Mois_cont)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[3, 2]/pcoa$values[1, 2]) +
  scale_color_gradient(name = "Moisture") +
  scale_shape(name = 'Site')

ggplot(sample_data(phyl), aes(x = pc3, y = pc4, color = Z, shape = geocode)) +
  geom_point(size = 3) +
  coord_fixed(pcoa$values[4, 2]/pcoa$values[3, 2]) +
  scale_color_gradient(name = "Distance") +
  scale_shape(name = 'Site')

# Environmental and spatial distance ################
# Environment
env_dis <- dist(as.matrix(sample_data(physeq)[, c('Mois_cont', 'Bulk_dens', 'N_percent', 'C_percent')]))
env_sim <- 1/(1+dist(as.matrix(sample_data(physeq)[, c('Mois_cont', 'Bulk_dens', 'N_percent', 'C_percent')])))
# Spatial
geo_dis <- dist(as.matrix(sample_data(physeq)[, c('X', 'Y')]))
geo_sim <- 1/(1+dist(as.matrix(sample_data(physeq)[, c('X', 'Y')])))
# Community
com_dis <- divnet$`bray-curtis`

# Test multi-collinearity of community, environment, and geography ###############
#com_dis <- vegdist(otu_table(physeq))
mantel(com_dis, env_dis)
mantel(com_dis, geo_dis)
mantel(env_dis, geo_dis)

# Test differential abundance of taxa associated with changes in methane oxidation rate
# Rerun the test with different sets of covariates to get at different aspects of
# community structure
fit_gce <- differentialTest(formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X + Y + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X + Y + pc1 + pc2 + pc3,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent + X + Y + pc1 + pc2 + pc3,
                        phi.formula_null = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X + Y + pc1 + pc2 + pc3,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)

plot(fit_gce)
ggsave('../figures/dif_abund_taxa.pdf', width = 8, height = 4)
otu_to_taxonomy(fit_gce$significant_taxa, phyl)
fit_gce$significant_taxa %in% gce_
gce_%in%fit_gce$significant_taxa
otu_to_taxonomy(fit_n$significant_taxa, phyl)
fit_gce$significant_models[1]

fit <- bbdml(formula = OTU8926 ~ pos_lowk + X + Y + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + X + Y + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        data = phyl)

fit_ge <- differentialTest(formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X + Y,
                        phi.formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X + Y,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent + X + Y,
                        phi.formula_null = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + X + Y,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit_ge)

length(fit_ge$significant_taxa)
fit$significant_taxa %in% fit1$significant_taxa
fit1$significant_taxa %in% fit$significant_taxa

fit_gc <- differentialTest(formula = ~ pos_lowk + X + Y + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk  + X + Y + pc1 + pc2 + pc3,
                        formula_null = ~  X + Y + pc1 + pc2 + pc3,
                        phi.formula_null = ~ pos_lowk + X + Y + pc1 + pc2 + pc3,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit_gc)



length(fit_gc$significant_taxa)

fit_ce <- differentialTest(formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent + pc1 + pc2 + pc3,
                        phi.formula_null = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent + pc1 + pc2 + pc3,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit_ce)

length(fit_ce$significant_taxa)


fit_e <- differentialTest(formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent,
                        phi.formula = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent,
                        formula_null = ~ Bulk_dens + Mois_cont +N_percent + C_percent,
                        phi.formula_null = ~ pos_lowk + Bulk_dens + Mois_cont + N_percent + C_percent,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit_e)

length(fit_e$significant_taxa)

fit_c <- differentialTest(formula = ~ pos_lowk + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + pc1 + pc2 + pc3,
                        formula_null = ~ pc1 + pc2 + pc3,
                        phi.formula_null = ~ pos_lowk + pc1 + pc2 + pc3,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit_c)

length(fit_c$significant_taxa)


system.time(fit_g <- differentialTest(formula = ~ pos_lowk + X + Y,
                        phi.formula = ~ pos_lowk + X + Y,
                        formula_null = ~ X + Y,
                        phi.formula_null = ~ pos_lowk + X + Y,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05))
plot(fit_g)

length(fit_g$significant_taxa)

system.time(fit_n <- differentialTest(formula = ~ pos_lowk,
                        phi.formula = ~ pos_lowk,
                        formula_null = ~ 1,
                        phi.formula_null = ~ pos_lowk,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05))
plot(fit_n)
length(fit_n$significant_taxa)


# Example #####################################

data(soil_phylo)
soil <- soil_phylo %>%
phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
phyloseq::tax_glom("Phylum")
da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ DayAmdmt,
                                test = "Wald", boot = FALSE,
                                data = soil,
                                fdr_cutoff = 0.05)
plot(da_analysis)
da_analysis$significant_models

# None ################
diagdds = phyloseq_to_deseq2(phyl, ~  pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)

n_ <- rownames(sigtab)
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1/8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



# G ################
diagdds = phyloseq_to_deseq2(phyl, ~ X + Y + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)

g_ <- rownames(sigtab)
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1/8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



# C ################
diagdds = phyloseq_to_deseq2(phyl, ~ pc1 + pc2 + pc3 + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)

c_ <- rownames(sigtab)
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1/8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



# E ################
diagdds = phyloseq_to_deseq2(phyl, ~ Mois_cont + Bulk_dens + N_percent + C_percent + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)
e_ <- rownames(sigtab)

sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


# GE ################
diagdds = phyloseq_to_deseq2(phyl, ~ X+Y+Mois_cont + Bulk_dens + N_percent + C_percent + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)
ge_ <- rownames(sigtab)

sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


# GC ################
diagdds = phyloseq_to_deseq2(phyl, ~ X+Y+pc1 + pc2 + pc3 + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)
gc_ <- rownames(sigtab)

sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


# CE ################
diagdds = phyloseq_to_deseq2(phyl, ~pc1 + pc2 + pc3 + N_percent + C_percent + Mois_cont + Bulk_dens + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)
ce_ <- rownames(sigtab)

sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


# GCE ###################
diagdds = phyloseq_to_deseq2(phyl, ~ Bulk_dens + Mois_cont + N_percent + C_percent + X + Y + pc1 + pc2 + pc3 + pos_lowk)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyl)[rownames(sigtab), ], "matrix"))
head(sigtab)
nrow(sigtab)
gce_ <- rownames(sigtab)

theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, xmin = log2FoldChange - lfcSE, 
                      xmax = log2FoldChange + lfcSE, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  geom_errorbarh(height = 1/8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
#n_ <- fit_n$significant_taxa
#g_ <- fit_g$significant_taxa
#c_ <- fit_c$significant_taxa
#e_ <- fit_e$significant_taxa
#gc_ <- fit_gc$significant_taxa
#ge_ <- fit_ge$significant_taxa
#ce_ <- fit_ce$significant_taxa
#gce_ <- fit_gce$significant_taxa

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



fit <- differentialTest(formula = ~ Wetland,
                        phi.formula = ~ Wetland,
                        formula_null = ~ 1 ,
                        phi.formula_null = ~ Wetland,
                        data = phyl,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit)
sample_data(phyl)$Wetland
sample_data(physeq)$Wetland


samples <- sample_data(phyl)[, c('pos_lowk', 'pc1', 'pc2', 'pc3')]
otus <- otu_table(phyl)[, colSums(otu_table(phyl)) > 0]
test <- phyloseq(samples, otus)
dim(otu_table(test))

saveRDS(test, '../output/test_physeq.rds')

fit <- differentialTest(formula = ~ pos_lowk + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + pc1 + pc2 + pc3,
                        formula_null = ~ pc1 + pc2 + pc3,
                        phi.formula_null = ~ pos_lowk + pc1 + pc2 + pc3,
                        data = test,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(fit)

length(fit2$significant_taxa)

test2 <- readRDS('../../../../Downloads/test_physeq.rds')

fit2 <- differentialTest(formula = ~ pos_lowk + pc1 + pc2 + pc3,
                        phi.formula = ~ pos_lowk + pc1 + pc2 + pc3,
                        formula_null = ~ pc1 + pc2 + pc3,
                        phi.formula_null = ~ pos_lowk + pc1 + pc2 + pc3,
                        data = test2,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
test2
test
sample_data(test) == 
sample_data(test2)
otu_table(test)
otu_table(test2)
