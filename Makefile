#!gmake
#
# Makefile for man pages
#
RM	=	/bin/rm -f
RMDIR	=	/bin/rm -rf
TOUCH	=	/usr/bin/touch
SHELL	=	/bin/csh

CHMOD	=	chmod 664
TEX	=	latex
BIB	=	bibtex
IDX	=	makeindex -s hfill.ist
DVIPS	=	dvips

LOCAL	=	$(GMXHOME)/src/local
HTML	=	$(GMXHOME)/html

FILES = intro		defunits	algorithms	par-md		\
	forcefield	topology	special		programs	\
	analyse		install		implement	tables		\
	lr-corr		averages	progman

AUXFILES = $(foreach FILE,$(FILES), $(FILE).aux)
TEXFILES = $(foreach FILE,$(FILES), $(FILE).tex)

all:		gromacs.ps

full:		man all

dvi:		gromacs.dvi

#
# make a booklet, i.e. 4 pages onto one double-sided page.
# To get the booklet, rearrange the pages according to page numbering 
# and fold in the middle
#
booklet:	gromacs.ps
		pstops "4:-3L@.7(21cm,0)+0L@.7(21cm,14.85cm),1L@.7(21cm,0)+-2L@.7(21cm,14.85cm)" $^ booklet.ps

gromacs.tex:	$(TEXFILES)

gromacs.aux:	gromacs.tex $(AUXFILES)
		$(TEX) gromacs

bib+idx:	gromacs.tex
		$(TEX) gromacs
		$(BIB) gromacs
		$(IDX) gromacs

gromacs.dvi:	bib+idx		gromacs.aux

gromacs.ps:	gromacs.dvi
		dvips -M -o $@ $^

%.aux:		%.tex

prog:		mdp_opt.tex proglist.tex

man:		./mkman

progman.tex:	
		$(TOUCH) progman.tex

mdp_opt.tex:	./mkmdp $(HTML)/progman.html
		./mkmdp $(GMXHOME)

proglist.tex:	$(LOCAL)/mkonline $(LOCAL)/programs.txt
		cd $(LOCAL) ; ./mkonline $(GMXHOME)

man:		
		mkman

clean:
		$(RM) *.log *.lof *.lot *.bbl *.blg *.toc *.dvi *.aux *.ps *~ #*# *.idx *.ilg *.ind 
		$(RM) progman.tex
		$(RMDIR) progman
