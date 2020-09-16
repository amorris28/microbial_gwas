# Linking microbial communities to ecosystem functions: what we can learn from genotype-phenotype mapping in organisms
This is the project directory for the Royal Society paper cited below. All intermediate
data files (in `output/`) and all figures (in `figures/`) have been
created so the final script `7_analysis.R` can be run immediately to recreate
figures and results. Otherwise, run scripts in order of numbering to recreate
the output.

## How to cite this manuscript:

Morris Andrew, Meyer Kyle, and Bohannan Brendan 2020. Linking microbial
communities to ecosystem functions: what we can learn from genotype–phenotype
mapping in organisms. [Phil. Trans. R. Soc.
B](https://royalsocietypublishing.org/journal/rstb) 375:20190244.
[https://doi.org/10.1098/rstb.2019.0244](https://doi.org/10.1098/rstb.2019.0244)

## Data are from:

Meyer, KM, Hopple, AM, Klein, AM, Morris, AH, Bridgham, SD, Bohannan, BJM.
Community structure – Ecosystem function relationships in the Congo Basin
methane cycle depend on the physiological scale of function. [Mol Ecol.](https://onlinelibrary.wiley.com/journal/1365294x) 2020; 29:
1806– 1819.
[https://doi.org/10.1111/mec.15442](https://doi.org/10.1111/mec.15442)


## Data

`raw_data/` contains raw data as received from colleagues. This is never modified.

`output/` is where modified data files are saved. 

`talapas-output/` is where model output from the talapas HPC cluster are saved.

`figures/` contains figures output from scripts

`R/` Contains scripts that run plus helper functions

## Scripts

### R/

`0_setup.R` has code to install all required packages. May need to modify to
address your local dependencies.

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

