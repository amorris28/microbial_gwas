
varCompGE <- function(model, drop = FALSE) {
  if (length(model$varComps) <= 1) {
h2G(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = drop))
  } else {
h2GE(c(model$varComps, err = model$sigma2), vcov(model, what = 'varComp', drop = drop))
  }
}

permute_no_replace_for <- function(perm_vec, nperm = 999) {
samples <- matrix(ncol = nperm, nrow = length(perm_vec))
samples <- cbind(samples, truth = perm_vec)
for (i in seq_len(nperm)) {
  sample_i <- sample(perm_vec)
  while (any(apply(samples, 2, function(x) {all(sample_i == x)}), na.rm = T)) {
    sample_i <- sample(perm_vec)
  }
  samples[, i] <- sample_i
  colnames(samples)[i] <- paste0('perm.', i)
}
return(samples)
}

#taxon_table$taxon <- apply(taxon_table, 1, function(r) { r[ 
#                      if( any(is.na(r))){ 
#                          (which(is.na(r))[1]-1) 
#                       }else{
#                          (length(r)-0)}
#                                      ] }
#       )
##### Functions

lowk_model <- function(df) {
  lm(Low_final_k ~ abund, data = df)
}
vmax_model <- function(df) {
  lm(Vmax ~ abund, data = df)
}
# Create nested data set by asv
  # Calculate linear model for each asv between abundance and function
fit_models <- function(data, model) {
  data %>% 
    gather(asv, abund, starts_with('asv_')) %>% 
    group_by(asv) %>% 
    nest() %>% 
    mutate(model = map(data, model)) %>% 
    return()
}
#### Functions
pc_adjust_mat <- function(raw_mat, aj, n_axes) {
  # raw_mat is the matrix of function and community
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

man_plot <- function(by_asv, taxon_table) {

by_asv %>% 
  mutate(tidy = map(model, broom::tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'abund') %>%
  left_join(taxon_table, by = 'asv') %>% 
  ggplot(aes(x = reorder(asv, Phylum, median), y = -log10(p.value))) +
  geom_point(aes(color=Phylum), size = 1) +
  geom_hline(yintercept = -log10(0.05/nrow(by_asv))) +
  scale_color_discrete(guide = F) +
  labs(x = 'Amplicon Sequence Variant (ASV)', 
       y = expression(Significance ~ (-log[10](p.value)))) +
    #  geom_text_repel(aes(label=ifelse(-log10(p.value)>-log10(alpha),as.character(genus),'')), size = 1) +
    # geom_text(aes(label=ifelse(-log10(p.value)>-log10(alpha),as.character(genus),'')), position = position_dodge(width = 10)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0, 10)  %>% 
  return()
}
# Identify ASVs with -log(p value) > 
id_asvs <- function(model_output, alpha = 0.05/nrow(model_output)) {
lowk_asvs <- model_output %>% 
  mutate(tidy = map(model, broom::tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'abund') %>% 
  filter(-log10(p.value) > -log10(alpha)) %>% 
  pull(asv)
}

# Pull out taxonomy for significant taxa
id_taxonomy <- function(sig_asvs, taxon_table) {
  taxon_table %>% 
    filter(asv %in% sig_asvs) 
}
