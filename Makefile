R_OPTS=--no-save --no-restore --no-init-file --no-site-file

manuscript/main.pdf: manuscript/main.tex
	cd $(@D); latexmk -pdf -M -MP -MF $*.deps $(<F)

LM_FIT := output/lowk_gits_c.rds

PC_ADJUST := output/gab_adj_com.tsv output/gab_adj_env.tsv output/gab_adj_spa.tsv
PC_ADJUST += output/gab_adj_com_env.tsv output/gab_adj_com_spa.tsv
PC_ADJUST += output/gab_adj_env_spa.tsv output/gab_adj_com_env_spa.tsv

$(PC_ADJUST): R/PC_correction.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

output/gab_troph_attr_table.csv output/gab_gen_attr_table.csv: R/gab_cleaning.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

output/gab_asv_table.csv: R/gab_clean_asvs.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

output/gab_rare_asv_table.csv: R/rarefy_asv_table.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

clean: 
	\rm -f *.aux *.bbl *.blg *.log *.bak *.Rout */*.Rout */*.aux */*.log

clean-latex:
	cd manuscript; latexmk -C

-include *.deps

