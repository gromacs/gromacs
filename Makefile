#
# Makefile for man pages
#
RM	=	/bin/rm -f
SHELL	=	/bin/csh

CHMOD	=	chmod 664
TEX	=	latex
BIB	=	bibtex
IDX	=	makeindex
DVIPS	=	dvips

all:		gromacs.ps

gromacs.aux:	
		$(TEX) gromacs

gromacs.bbl:	gromacs.aux
		$(BIB) gromacs
		$(TEX) gromacs

gromacs.idx:	gromacs.aux
		$(IDX) gromacs

gromacs.dvi:	gromacs.aux	gromacs.bbl	gromacs.idx
		$(TEX) gromacs

gromacs.ps:	gromacs.dvi
		dvips -M -o $@ $^

clean:
		$(RM) *.log *.lof *.lot *.bbl *.blg *.toc *.dvi *.aux *~ #*# *.idx *.ilg *.ind



