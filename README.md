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

`data/` contains raw data as received from colleagues. This is never modified,
but is under version control.

`output/` is where modified data files are saved. These are not version
controlled and she be recreated from source.

`talapas-output/` is where model output from the talapas HPC cluster.

## Scripts

### R
`R/` directory contains `R` scripts that run.

`clean_asvs.R` takes in the raw ASV table and munges it for analysis including
performing the `varianceStabilizingTransformation` from `DESeq2`.

`clean_metadata.R` takes in the raw sample date file and munges it for analysis.

`geodist_processing.R` takes in the raw GPS data and converts it to meters using
UTM.

`combine_data.R` pulls in the asvs, sample data, and GPS data and combines them
into a single data structure.

`var_comp.R` is the script that fits variance component models from the raw
data. This is typically submitted to talapas and the output from talapas is
`scp`ed into `talapas-output/`.

`analysis.R` is the primary file that presents the results and creates the
figures and tables. This script can be rendered into `html` using `knitr::render`. 

`functions.R` contains holds the user-made functions specific to this project.

### Manuscript
`manuscript/` directory contains `.tex` files to compile the manuscript

`ecology.bst` is the Ecology journal bibliography style file.

`latexmkrc` removes margins for easier viewing of the PDF on computers. This can
be renamed or removed to return to normal `latex` formatting.

`library.bib` is the bibliography file.

`main.tex` contains the manuscript text.

### Tables and Figures
These directories are not version controlled. They should always we able to be
recreated from source.
`tables/` is where tables are saved from `R` and imported to `latex`
`figures/` is where figures are saved from `R` and imported to `latex`
