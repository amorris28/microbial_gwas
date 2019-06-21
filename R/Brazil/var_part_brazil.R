library(tidyverse)
library(vegan)

all_data <- read_tsv('output/cleaned_data.tsv')
dim(all_data)
all_data <- all_data[complete.cases(all_data$Lat), ]
dim(all_data)
all_data <- all_data[complete.cases(all_data$Long), ]
dim(all_data)
colnames(all_data[, 1:40])

spa <- select(all_data, Lat, Long)
geo_dist <- as.matrix(dist(scale(as.matrix(spa))))

env <- select(all_data, pH:Total_moisture)
env_dissim <- as.matrix(dist(scale(as.matrix(env))))

com <- select(all_data, starts_with('asv'))
com_dissim <- as.matrix(vegdist(as.matrix(com)))

fun <- select(all_data, CH4)
fit <- varpart(fun$CH4, ~ geo_dist, ~ env_dissim)
fit
plot(fit)
