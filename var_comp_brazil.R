library(tidyverse)
library(varComp)
library(vegan)
#library(broom)
#library(data.table)
library(amorris)


## Import and organize data

all_data <- read.table('output/dim_cleaned_data.tsv', header = TRUE)
#all_data <- fread('output/dim_cleaned_data.tsv')
#all_data[1:nrow(all_data), 1:40]
#colnames(all_data[, 1:40])
#sum(is.na(all_data))
all_data <- all_data[complete.cases(all_data$Lat), ]
all_data <- all_data[complete.cases(all_data$Long), ]
asv_table <- select(all_data, starts_with("asv"))
asv_table[asv_table>0] <- 1 # Presence/Absence
asv_mat <- as.matrix(asv_table)
all_data <- all_data[rowSums(asv_mat) > 1, ]
all_data <- cbind(all_data[, 1:39], all_data[, colnames(asv_mat[, colSums(asv_mat) > 1])])
#all_data[1:nrow(all_data), 1:ncol(all_data)]

### Plot ASV abundances
#
#i = 2
#j = ncol(asv_table)
#  asv_table[, c(1, i:j)] %>% 
#  gather(ASV, ABU, 2:ncol(asv_table[, c(1, i:j)])) %>% 
#ggplot(aes(y = ABU, x = ASV)) +
#  geom_jitter()  


### Estimate and plot variation in ASV abundance
#
#ColVar <- function(x, ...) {
#  colSums((x - colMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
#}
#asv_mat_var <- ColVar(asv_mat)
#plot(density(log10(asv_mat_var)))
#asv_mat_high_var<-asv_mat[,log10(asv_mat_var) > 0]
#ncol(asv_mat_high_var)
#rowSums(asv_mat)
#colSums(asv_mat)
#asv_mat_ns <- asv_mat[rowSums(asv_mat) > 1, ]
#asv_mat_ns <- asv_mat[, colSums(asv_mat) > 1]



# ASV matrix
taxa <- select(select(all_data, starts_with('asv')), -1)
# Pull out a single predictor taxon
taxon <- taxa[, 5437]
# Calculate similarity from all taxa except that one
com_dissim <- vegdist(taxa[, -5437])
com_dissim <- as.matrix(com_dissim)
com_sim <- 1-com_dissim

# Methane flux response variable
CH4 <- all_data$CH4
hist(CH4)
# Calculate geographic similarity from lat/long
xy <- all_data[, c('Lat', 'Long')]
geo_dis <- as.matrix(dist(as.matrix(scale(xy))))
geo_sim <- 1/(1+geo_dis)

# Calculate environmental similarity from environmental variables
env <- select(all_data, pH:Total_moisture)
env_dis <- as.matrix(dist(as.matrix(scale(env))))
env_sim <- 1/(1+env_dis)


model <- varComp(CH4 ~ taxon, varcov = list(com = com_sim, geo = geo_sim, env = env_sim) )
summary(model)


for (i in 11:20) {
  taxon <- asv_mat_ns[, i]
  com_sim <- 1- as.matrix(vegdist(asv_mat_ns[, -i]))
  model <- varComp(lowk ~ taxon, varcov = list(geo = geo_sim, env = env_sim))
  print(summary(model))
        }


# Environment
colnames(troph_attr_table)
summary(lm(Low_final_k ~ WFPS + N_percent + C_percent, data = troph_attr_table))

# Community
colSums(asv_mat)

PCA <- prcomp(asv_mat, center = T, scale. = T)
screeplot(PCA)
summary(PCA)
PC1 <- scores(PCA)[, 1]
summary(lm(lowk ~ PC1))
plot(lowk ~ PC1)
PC2 <- scores(PCA)[, 2]
summary(lm(lowk ~ PC2))
plot(lowk ~ PC2)
PC3 <- scores(PCA)[, 3]
summary(lm(lowk ~ PC3))
plot(lowk ~ PC3)
PC4 <- scores(PCA)[, 4]
summary(lm(lowk ~ PC4))
plot(lowk ~ PC4)
PC5 <- scores(PCA)[, 5]
summary(lm(lowk ~ PC5))
plot(lowk ~ PC5)
PC6 <- scores(PCA)[, 6]
summary(lm(lowk ~ PC6))
plot(lowk ~ PC6)



tail(PCA$rotation[, 4][order(abs(PCA$rotation[, 4]))])
tail(PCA$rotation[, 5][order(abs(PCA$rotation[, 5]))])
# Space
plot(lowk ~ xy[, 1])
plot(lowk ~ xy[, 2])
xyz <- cbind(xy, Z = xy[, 1] * xy[ ,2])
plot(lowk ~ z)
summary(lm(lowk ~ xy[, 1] * xy[, 2]))
xyz <- matrix(xy[, 1], xy[, 2], xy[, 1] * xy[, 2], nrow = 44, ncol = 3 )
RDA <- rda(asv_mat ~ xyz, scale = T)
plot(RDA)
summary(RDA)
RDA


# Abundance

summary(lm(Low_final_k ~ pmoa_copy_num, troph_attr_table))
summary(lm(Low_final_k ~ pmoa_transcript_num, troph_attr_table))
