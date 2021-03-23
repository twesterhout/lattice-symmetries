all: pdf jpg

pdf: 01_basis_construction.pdf 02_operator_application.pdf
jpg: 01_basis_construction.jpg 02_operator_application.jpg

%.pdf: %.tex cn71/*.dat
	latexrun $<
	rm -r latex.out/

%.jpg: %.pdf
	pdftoppm -jpeg -r 600 $< $*
	mv $*-1.jpg $*.jpg

clean: *.pdf *.jpg latex.out/
	rm -rf latex.out
	rm -f *.pdf
	rm -f *.jpg

.PHONY: clean
