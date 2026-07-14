# chapter/Makefile
DOC = Chapter
.PHONY: all clean
all: $(DOC).pdf
$(DOC).pdf: $(DOC).tex sections/*.tex References.bib
	pdflatex -interaction=nonstopmode $(DOC).tex
	-bibtex $(DOC)
	pdflatex -interaction=nonstopmode $(DOC).tex
	pdflatex -interaction=nonstopmode $(DOC).tex
clean:
	rm -f $(DOC).aux $(DOC).bbl $(DOC).blg $(DOC).log $(DOC).out $(DOC).pdf
