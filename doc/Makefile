MAIN=optimizationTutorial
MAINPDF:=$(MAIN:=.pdf)

CC=latexmk

.PHONY:	all run clean cleanall

all:	$(MAINPDF)

%.pdf:	%.tex
	$(CC) -pdf $<

run:	$(MAINPDF)
	okular 2>/dev/null $(MAINPDF) &

clean:
	$(CC) -c

cleanall:	clean
	$(CC) -C
