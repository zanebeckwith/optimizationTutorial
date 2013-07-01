MAIN=optimizationTutorial
MAINPDF:=$(MAIN:=.pdf)

ifeq ($(findstring physics.unc.edu,$(shell hostname -f)),physics.unc.edu)
CC=pdflatex
CLEAN=rm -f *.log *.aux *.out
CLEANALL=rm -f $(MAINPDF)
else
CC=latexmk -pdf
CLEAN=$(CC) -c
CLEANALL=$(CC) -C
endif

.PHONY:	all run clean cleanall

all:	$(MAINPDF)

%.pdf:	%.tex
	$(CC) $<

run:	$(MAINPDF)
	okular 2>/dev/null $(MAINPDF) &

clean:
	$(CLEAN)

cleanall:	clean
	$(CLEANALL)
