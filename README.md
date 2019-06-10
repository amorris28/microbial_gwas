# Royal Society Paper
- Processes data from Gabon and Brazil
- Performs variance partitioning
- Performs GWAS-esque CWAS analysis
- Compiles the manuscript
- Also compiles a presentation

## Terms
Dim, Dimensions, Brazil - Files associated with the Dimensions project in Brazil

Gab, Gabon, gab - Files associated with the Gabon project

adj, adjusted - Data that has had the PC correction performed

## Scripts
`brazil_cleaning.R` - Process Brazil data and combine into one data file
including the ASV table, environmental metadata and lat long coordinates

`gabon_cleaning.R` - Process Gabon data and combine into one data file including
ASV table, environmental metadata, and lat long coordiantes

`geodist_processing.R` - Takes Gabon lat long coordinates and converts them to
meters in UTM

`rarefy_asv_table.R` - Rarefies Gabon ASV table

`lm_model.R` - Takes adjusted function and community data and performs linear
models for all taxa

`id_taxa.R` - Takes Gabon GWAS model output and returns tables of identified
taxa

`var_comp_brazil.R` - 
