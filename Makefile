#
# Makefile for man pages
#
RM	=	/bin/rm -f
SHELL	=	/bin/csh

CHMOD	=	chmod 664
TEX	=	latex
BIB	=	bibtex
IDX	=	makeindex -s flush.ist
DVIPS	=	dvips

AUXFILES = 	algorithms.aux	analyse.aux	averages.aux	\
		defunits.aux	forcefield.aux	implement.aux	\
		install.aux	intro.aux	lr-corr.aux	\
		special.aux	tables.aux	topology.aux

all:		gromacs.ps

dvi:		gromacs.dvi

#
# make a booklet, i.e. 4 pages onto one double-sided page.
# To get the booklet, rearrange the pages according to page numbering 
# and fold in the middle
#
booklet:	gromacs.ps
		pstops "4:-3L@.7(21cm,0)+0L@.7(21cm,14.85cm),1L@.7(21cm,0)+-2L@.7(21cm,14.85cm)" $^ booklet.ps

gromacs.aux:	gromacs.tex $(AUXFILES)
		$(TEX) gromacs

gromacs.bbl:	*.tex
		$(BIB) gromacs
		$(TEX) gromacs

gromacs.ind:	*.tex
		$(IDX) gromacs
		$(TEX) gromacs

gromacs.dvi:	gromacs.aux	gromacs.bbl	gromacs.ind
		$(TEX) gromacs

gromacs.ps:	gromacs.dvi
		dvips -M -o $@ $^

algorithms.aux:	algorithms.tex
		$(TEX) gromacs

analyse.aux:	analyse.tex
		$(TEX) gromacs

averages.aux:	averages.tex
		$(TEX) gromacs

defunits.aux:	defunits.tex
		$(TEX) gromacs

forcefield.aux:	forcefield.tex
		$(TEX) gromacs

implement.aux:	implement.tex
		$(TEX) gromacs

install.aux:	install.tex
		$(TEX) gromacs

intro.aux:	intro.tex
		$(TEX) gromacs

lr-corr.aux:	lr-corr.tex
		$(TEX) gromacs

special.aux:	special.tex
		$(TEX) gromacs

tables.aux:	tables.tex
		$(TEX) gromacs

topology.aux:	topology.tex
		$(TEX) gromacs

clean:
		$(RM) *.log *.lof *.lot *.bbl *.blg *.toc *.dvi *.aux *~ #*# *.idx *.ilg *.ind



