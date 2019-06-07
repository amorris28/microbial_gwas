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

`combine_data.R` - Process Gabon data and combine into one data file including
ASV table, environmental metadata, and lat long coordiantes

`geodist_processing.R` - Takes Gabon lat long coordinates and converts them to
meters in UTM
