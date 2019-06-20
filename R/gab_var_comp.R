library(tidyverse)
library(varComp)
library(vegan)
library(broom)
library(morris)

all_data <- read.csv('../output/gab_all_rare_troph.csv', stringsAsFactors = FALSE)

## Remove ASVs only present in 1 sample
asvs <- as.matrix(all_data[, grep("^asv", colnames(all_data))])
asvs[asvs>0] <- 1
all_data <- all_data[rowSums(asvs) > 1, ]
asv_mat <- as.matrix(all_data[, grep("^asv", colnames(all_data))])
dim(asvs); dim(asv_mat)

# Community similarity
asv <- 3
taxon <- asv_mat[, asv]
com_sim <- 1 - as.matrix(vegdist(asv_mat[, -asv]))

env_variables <- all_data[, 1:8]
# Function
lowk <- all_data$Low_final_k
hist(lowk)
hist(taxon)
# Geographic similarity
xy <- all_data[, c('X', 'Y')]
geo_dis <- as.matrix(dist(as.matrix(scale(xy))))
geo_sim <- 1/(1+geo_dis)

# environmental dissimilarity
wetland <- all_data[, "Wetland"]
env <- all_data[, c('WFPS', 'N_percent', 'C_percent')]
env_dis <- as.matrix(dist(as.matrix(scale(env))))
env_sim <- 1/(1+env_dis)

mantel(com_sim, geo_sim)
qplot(x = xy$X, y = xy$Y)

# Community and geography are correlated and sites cluster around three locations
# Also, geography uses up more of the similarity matrix space (see histograms)
# Therefore, I coded geography as three factors (one for each cluster of sites)
# Also, wetland/upland covaries with those three sites, but moisture content is
# included in the environmental similarity matrix and should capture that variation
all_data[all_data$Y > 30000, 'geocode'] <- 'A'
all_data[all_data$Y < 30000 & all_data$Y > 10000, 'geocode'] <- 'B'
all_data[all_data$Y < 10000, 'geocode'] <- 'C'
all_data$geocode <- factor(all_data$geocode)
hist(com_sim)
hist(env_sim)
hist(geo_sim)
qplot(x = xy$X, y = xy$Y) + geom_jitter(aes(color = all_data$Wetland))

model <- varComp(lowk ~ taxon + all_data$geocode, varcov = list(com = com_sim, env = env_sim) )
summary(model)

fit_varcomps <- function(x) {
  print(x)
  taxon <- asv_mat[, x]
  com_sim <- 1 - as.matrix(vegdist(asv_mat[, -x]))
  model <- varComp(Low_final_k ~ taxon + geocode, data = all_data, 
                   varcov = list(com = com_sim, env = env_sim))
  return(c(model$varComps, model$sigma2)/sum(c(model$varComps, model$sigma2)))
}


#system.time(suppressWarnings(mycomps <- lapply(seq(ncol(asv_mat[, 1:100])), fit_varcomps)))

