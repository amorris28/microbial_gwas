# Royal Society Paper
This is the project directory for the Royal Society paper. 

How to cite this manuscript:

Morris AH, Meyer KM, Bohannan BJM (2019) Linking microbial communities
to ecosystem functions: what we can learn from genotype-phenotype mapping
in organisms 

Data are from:

Meyer KM, Hopple AM, Klein A, Morris AH, Bridgham
SD, Bohannan BJM (2019) Community structure – ecosystem function
relationships in the Congo Basin methane cycle depend on the physiological
scale of function. 
[bioRxiv](https://www.biorxiv.org/) doi:[10.1101/639989](https://doi.org/10.1101/639989)


## Data

`raw_data/` contains raw data as received from colleagues. This is never modified.

`output/` is where modified data files are saved. 

`talapas-output/` is where model output from the talapas HPC cluster are saved.

`Figures/` contains figures output from scripts

`R/` Contains scripts that run plus helper functions

## Scripts

### R
`R/` directory contains `R` scripts that run.

`0_setup.R` has code to install all required packages.

`1_clean_sample_data.R` takes in the raw sample date file and munges it for analysis.

`2_clean_asvs.R` takes in the raw ASV table and munges it for analysis including

`3_geodist_processing.R` takes in the raw GPS data and converts it to meters using
UTM.

`4_var_stab_transform.R` performs the `varianceStabilizingTransformation` from `DESeq2`.

`5_combine_data.R` pulls in the asvs, sample data, and GPS data and combines them
into a single data structure.

`6_var_comp.R` is the script that fits variance component models from the raw
data. This is typically submitted to the cluster and the output from the cluster is
`scp`ed into `talapas-output/`.

`7_analysis.R` is the primary file that presents the results and creates the
figures and tables.

`functions.R` contains holds the user-made functions specific to this project.

