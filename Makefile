
main.pdf: main.tex
	latexmk -pdf -M -MP -MF $*.deps $<

clean-latex:
	latexmk -C

-include *.deps
