# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

cran_packages <- c('devtools', 'tidyverse', 'sp', 'rgdal', 'broom', 'vegan',
                   'foreach', 'doParallel', 'data.table',  
                   'kableExtra', 'gap', 'ape')
github_packages <- c('amorris28/morris', 'kbroman/broman', "adw96/breakaway")
bioc_packages <- c('DESeq2', 'phyloseq', 'qvalue')

install.packages(cran_packages)
devtools::install_version("varComp", version = "0.2-0", repos = "http://cran.us.r-project.org")
BiocManager::install(bioc_packages)
devtools::install_github(github_packages)
