
all: manual_5.2.pdf 

manual_5.2.pdf:manual_5.2.tex
	latex  manual_5.2
	bibtex manual_5.2
	latex  manual_5.2
	latex  manual_5.2
	dvipdf  manual_5.2.dvi
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg