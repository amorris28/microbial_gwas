
manuscript/main.pdf: manuscript/main.tex
	cd $(@D); latexmk -pdf -M -MP -MF $*.deps $(<F)

clean-latex:
	cd manuscript; latexmk -C

-include *.deps
