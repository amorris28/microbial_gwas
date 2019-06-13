R_OPTS=--no-save --no-restore --no-init-file --no-site-file

manuscript/main.pdf: manuscript/main.tex
	cd $(<D); latexmk -pdf -M -MP -MF $*.deps $(<F)

output/lowk_fits%.rds: R/gab_lm_model.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

R/gab_lm_model.R: output/gab_adj%.tsv

output/gab_adj%.tsv: R/PC_correction.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

output/gab_troph_attr_table.csv output/gab_gen_attr_table.csv: R/gab_cleaning.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

output/gab_asv_table.csv: R/gab_clean_asvs.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

output/gab_rare_asv_table.csv: R/rarefy_asv_table.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

clean: 
	\rm -f *.aux *.bbl *.blg *.log *.bak *.Rout */*.Rout */*.aux */*.log

cleanall:
	\rm  -f *.aux *.bbl *.blg *.log *.bak *.Rout */*.Rout */*.aux */*.log *.pdf Figures/*.pdf

clean-latex:
	cd manuscript; latexmk -C

-include *.deps

