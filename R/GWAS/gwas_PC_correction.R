#### Calculate community and environment similarity scores using PCA
#### with the vegan::rda function. Then adjust the community matrix
#### and ecosystem function rates by correcting for the regression
#### coefficients calculated from genotype ~ community similarity,
#### phenotype ~ community similarity, genotype ~ environmental
#### similarity, and phenotype ~ environmental similarity
#### See: Price et al. 2006 in Nature Genetics

library(vegan)
source('functions.R')
# Whether to correct by each factor
com_bool <- T
env_bool <- T
spa_bool <- T
# Set the number of axes to correct by
n_com_axes <- 5 
n_env_axes <- 2 
n_spa_axes <- 1
# Whether to center/scale the data
scale_com <- T
scale_env <- T
scale_spa <- T
# File names
input_sep <- ","
input_file <- "../output/gab_troph_total.csv"
output_file <- "../output/gab_adj"
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


all_data <- read.table(input_file, header = TRUE, sep = input_sep)

# Remove observations with missing location information
all_data <- all_data[complete.cases(all_data[, spa_cols[1]]), ]
all_data <- all_data[complete.cases(all_data[, spa_cols[2]]), ]

non_asv_cols <- ncol(all_data[, !grepl(com_cols, names(all_data))])
asv_mat <- as.matrix(all_data[, grepl(com_cols, names(all_data))])
asv_mat[asv_mat > 0] <- 1
asv_mat <- asv_mat[, colSums(asv_mat) > 1]
all_data <- cbind(all_data[, 1:non_asv_cols], 
                  all_data[, colnames(asv_mat)])
n_samples <- nrow(all_data)

# Extract ASV table
com_mat <- as.matrix(all_data[, grepl(com_cols, names(all_data))])

fun_mat <- as.matrix(all_data[, fun_cols])
fun_com <- cbind(fun_mat, com_mat)


### Community similarity correction
if (com_bool) {
#center_com <- data.frame(scale(community, scale = FALSE))
# Initial parameterization
com_pca <- rda(com_mat, scale = scale_com) # Calculate covariance matrix
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
all_adj <- pc_adjust_mat(fun_com, com_aj, n_com_axes) 
write.table(all_adj, paste0(output_file, '_c.tsv'))
}

if (env_bool) {
all_adj <- pc_adjust_mat(fun_com, env_aj, n_env_axes) 
write.table(all_adj, paste0(output_file, '_e.tsv'))
}

if (spa_bool) {
all_adj <- pc_adjust_mat(fun_com, spa_aj, n_spa_axes) 
write.table(all_adj, paste0(output_file, '_s.tsv'))
}

if (com_bool & env_bool) {
all_adj <- pc_adjust_mat(fun_com, com_aj, n_com_axes) 
all_adj <- pc_adjust_mat(all_adj, env_aj, n_env_axes) 
all_adj <- cbind(all_data[, group_cols], all_adj)
write.table(all_adj, paste0(output_file, '_ce.tsv'))
}
if (com_bool & spa_bool) {
all_adj <- pc_adjust_mat(fun_com, com_aj, n_com_axes) 
all_adj <- pc_adjust_mat(all_adj, spa_aj, n_spa_axes) 
all_adj <- cbind(all_data[, group_cols], all_adj)
write.table(all_adj, paste0(output_file, '_cs.tsv'))
}
if (env_bool & spa_bool) {
all_adj <- pc_adjust_mat(fun_com, env_aj, n_env_axes) 
all_adj <- pc_adjust_mat(all_adj, spa_aj, n_spa_axes) 
all_adj <- cbind(all_data[, group_cols], all_adj)
write.table(all_adj, paste0(output_file, '_es.tsv'))
}
if (com_bool & env_bool & spa_bool) {
all_adj <- pc_adjust_mat(fun_com, com_aj, n_com_axes) 
all_adj <- pc_adjust_mat(all_adj, env_aj, n_env_axes) 
all_adj <- pc_adjust_mat(all_adj, spa_aj, n_spa_axes) 
all_adj <- cbind(all_data[, group_cols], all_adj)
write.table(all_adj, paste0(output_file, '_ces.tsv'))
}

write.table(all_data, paste0(output_file, '_raw.tsv'))
