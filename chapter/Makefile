# chapter/Makefile
DOC = Chapter
FIGURES := $(wildcard figures/*.pdf)
.PHONY: all clean
all: $(DOC).pdf
$(DOC).pdf: $(DOC).tex sections/*.tex References.bib $(FIGURES)
	pdflatex -interaction=nonstopmode $(DOC).tex
	-bibtex $(DOC)
	pdflatex -interaction=nonstopmode $(DOC).tex
	pdflatex -interaction=nonstopmode $(DOC).tex
clean:
	rm -f $(DOC).aux $(DOC).bbl $(DOC).blg $(DOC).log $(DOC).out $(DOC).pdf
