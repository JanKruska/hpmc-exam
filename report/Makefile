#Settings
LATEXMK ?= latexmk
LATEXMKFLAGS = -use-make -xelatex -file-line-error -bibtex

#Set targets
SOURCES := $(wildcard *.tex)
FILES := $(SOURCES:.tex=.pdf)

all:$(FILES)

%.pdf: %.tex
	$(LATEXMK) $(LATEXMKFLAGS) $<

clean:
	$(LATEXMK) $(LATEXMKFLAGS) -C