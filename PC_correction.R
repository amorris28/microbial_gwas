#### Calculate community and environment similarity scores using PCA
#### with the vegan::rda function. Then adjust the community matrix
#### and ecosystem function rates by correcting for the regression
#### coefficients calculated from genotype ~ community similarity,
#### phenotype ~ community similarity, genotype ~ environmental
#### similarity, and phenotype ~ environmental similarity
#### See: Price et al. 2006 in Nature Genetics

# Whether to correct by each factor
com_bool <- T
env_bool <- T
spa_bool <- T
# Set the number of axes to correct by
n_com_axes <- 6 
n_env_axes <- 3 
n_spa_axes <- 1
# Whether to center/scale the data
scale_com <- T
scale_env <- T
scale_spa <- T
# File names
input_sep <- ","
input_file <- "output/troph_total.csv"
output_csv <- "output/adj_data.csv"
com_output_csv <- 'output/com_adj_data.csv'
env_output_csv <- 'output/env_adj_data.csv'
spa_output_csv <- 'output/spa_adj_data.csv'
# Function column names
fun_cols <- c('Low_final_k', 'Vmax')
# Environment column names
env_cols <- c('WFPS', 'N_percent', 'C_percent')
# Spatial column names
spa_cols <- c('X', 'Y')
# For community data, this assumes the ASV columns all start with 'asv'
com_cols <- 'asv'
# Extra grouping variables to preserve
group_cols <- c('Sample', 'Wetland', 'Site')

#library(tidyverse)
library(vegan)

all_data <- read.table(input_file, header = TRUE, sep = input_sep)

all_data <- all_data[complete.cases(all_data[, spa_cols[1]]), ]
all_data <- all_data[complete.cases(all_data[, spa_cols[2]]), ]

non_asv_cols <- ncol(all_data[, !grepl( "asv" , names(all_data))])
asv_mat <- as.matrix(all_data[, grepl( "asv" , names(all_data))])
asv_mat[asv_mat > 0] <- 1
asv_mat <- asv_mat[, colSums(asv_mat) > 1]
all_data <- cbind(all_data[, 1:non_asv_cols], 
                  all_data[, colnames(asv_mat)])
n_samples <- nrow(all_data)

# Extract ASV table
community <- as.matrix(all_data[, grepl( "asv" , names(all_data))])
com_mat <- t(community)

# Initialize function vector
fun_mat <- as.matrix(all_data[, fun_cols])

fun_com <- cbind(fun_mat, community)
class(fun_com)

#### Functions

pc_adjust_mat <- function(raw_mat, aj, n_axes) {
  # com_mat is the community matrix
  # aj is the matrix of site scores from PCA
  # n_axes is the number of PC axes to correct by
  gamma <- vector(length = nrow(raw_mat))
	for (p in 1:n_axes) {
    a <- aj[, p] # select axis n for correction
    for (i in 1:nrow(raw_mat)) { # adjust each taxon abundance
      tax <- raw_mat[i, ] # select taxon abundances for j samples
      gamma[i] <- sum(a * tax)/sum(a^2) # calculate gamma (regression coef)
      raw_mat[i, ] <- tax - gamma[i] * a 
    }
  }
  return(raw_mat)
}
pc_adjust_mat2 <- function(raw_mat, aj, n_axes) {
  # com_mat is the community matrix
  # aj is the matrix of site scores from PCA
  # n_axes is the number of PC axes to correct by
  gamma <- vector(length = nrow(raw_mat))
	for (p in 1:n_axes) {
    a <- aj[, p] # select axis n for correction
    for (i in 1:ncol(raw_mat)) { # adjust each taxon abundance
      tax <- raw_mat[, i] # select taxon abundances for j samples
      gamma[i] <- sum(a * tax)/sum(a^2) # calculate gamma (regression coef)
      raw_mat[, i] <- tax - gamma[i] * a 
    }
  }
  return(raw_mat)
}
pc_adjust_vec <- function(raw_vec, aj, n_axes) {
  for (p in 1:n_axes) {
    a <- aj[, p]
    gamma <- sum(a * raw_vec)/sum(a^2)
    raw_vec <- raw_vec - gamma * a
  }
  return(raw_vec)
}

### Community similarity correction
if (com_bool) {
#center_com <- data.frame(scale(community, scale = FALSE))
# Initial parameterization
com_pca <- rda(community, scale = scale_com) # Calculate covariance matrix
#plot(com_pca, choices = c(1, 2))
#screeplot(com_pca)
com_aj <- summary(com_pca)$sites # pull out site scores, which represent
# the community similarity matrix (shared ancestry in GWAS language)
# for community similarity predicting genotype
}
#### Correct for environmental similarity
if (env_bool) {
	
env <- all_data[, env_cols]
#center_env <- data.frame(scale(env, scale = FALSE))
env_pca <- rda(env, scale = scale_env)
#plot(env_pca, choices = c(1, 2))
#screeplot(env_pca)
env_aj <- summary(env_pca)$sites # pull out site scores, which represent
# the community similarity matrix (shared ancestry in GWAS language)
#####
}
if (spa_bool) {
spa <- all_data[, spa_cols]
# spa_aj <- as.matrix(spa_aj)
##center_env <- data.frame(scale(env, scale = FALSE))
spa_pca <- rda(spa, scale = scale_spa)
#screeplot(spa_pca)
spa_aj <- summary(env_pca)$sites # pull out site scores, which represent
## the community similarity matrix (shared ancestry in GWAS language)
}
if (com_bool) {
# Community correction
com_adj <- pc_adjust_mat(com_mat, com_aj, n_com_axes) 
# Function correction
func_adj <- apply(func, 2, pc_adjust_vec, com_aj, n_com_axes) 
func_adj2 <- pc_adjust_mat2(func, com_aj, n_com_axes) 
com_adj2 <- pc_adjust_mat2(community, com_aj, n_com_axes) 
all_adj <- pc_adjust_mat2(fun_com, com_aj, n_com_axes) 
cbind(func_adj, func_adj2)
sum(func_adj != func_adj2)
# Recombine data
func_adj_a <- as.data.frame(func_adj)
com_adj_a <- as.data.frame(t(com_adj))
new_data <- cbind(func_adj_a, com_adj_a)


write.table(new_data, com_output_csv, sep = ",")
}

if (env_bool) {
com_adj <- pc_adjust_mat(com_mat, env_aj, n_env_axes) 

func_adj <- apply(func, 2, pc_adjust_vec, env_aj, n_env_axes) 

# Recombine data
func_adj_a <- as.data.frame(func_adj)
com_adj_a <- as.data.frame(t(com_adj))
new_data <- cbind(func_adj_a, com_adj_a)


write.table(new_data, env_output_csv, sep = ",")
}
if (spa_bool) {
com_adj <- pc_adjust_mat(com_mat, spa_aj, n_spa_axes) 

func_adj <- apply(func, 2, pc_adjust_vec, spa_aj, n_spa_axes) 
# Recombine data
func_adj_a <- as.data.frame(func_adj)
com_adj_a <- as.data.frame(t(com_adj))
new_data <- cbind(func_adj_a, com_adj_a)


write.table(new_data, spa_output_csv, sep = ",")
}
if (com_bool & env_bool & spa_bool) {
## Community  
# for community similarity predicting genotype
com_adj <- pc_adjust_mat(com_mat, com_aj, n_com_axes) 

# Same correction as the com matrix now for vector of process rates
func_adj <- apply(func, 2, pc_adjust_vec, com_aj, n_com_axes) 
## Environment
com_adj <- pc_adjust_mat(com_adj, env_aj, n_env_axes) 

func_adj <- apply(func_adj, 2, pc_adjust_vec, env_aj, n_env_axes) 
## Spatial
com_adj <- pc_adjust_mat(com_adj, spa_aj, n_spa_axes) 

func_adj <- apply(func_adj, 2, pc_adjust_vec, spa_aj, n_spa_axes) 
# Recombine data
func_adj_a <- as.data.frame(func_adj)
com_adj_a <- as.data.frame(t(com_adj))
new_data <- cbind(func_adj_a, all_data[, group_cols], com_adj_a)


write.table(new_data, output_csv, sep = ",")
}

pc_correction <- function(com_mat, func_mat) {
  com_adj <- pc_adjust_mat(com_mat, com_aj, n_com_axes) 
  func_adj <- apply(func, 2, pc_adjust_vec, com_aj, n_com_axes) 

}
