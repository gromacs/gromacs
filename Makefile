#!gmake
#
# Makefile for man pages
#
RM	=	/bin/rm -f
RMDIR	=	/bin/rm -rf
TOUCH	=	touch
SHELL	=	/bin/csh -f

CHMOD	=	chmod 664
TEX	=	latex
BIB	=	bibtex
IDX	=	makeindex -s hfill.ist
DVIPS	=	dvips

LOCAL	=	$(GMXHOME)/src/local
HTML	=	$(GMXHOME)/html
COPYRGT	=	$(LOCAL)/copyrgt

TEXFS = algorithms	analyse		averages			\
	defunits	files		forcefield	ieee		\
	implement	install		intro		lr-corr		\
	lr_elstat	macros		mdp_opt		par-md		\
	proglist	progman		programs	special		\
	sqrt		tables		topolfig	topology	\
	virial

AUXFS = algorithms	analyse		averages	defunits	\
	forcefield	implement	install		intro		\
	lr-corr		progman		programs	special		\
	topology	

AUXFILES = $(foreach FILE,$(AUXFS), $(FILE).aux)
TEXFILES = $(foreach FILE,$(TEXFS), $(FILE).tex)

all:		ps

ps:		gromacs.ps
pdf:		gromacs.ps
		ps2pdf gromacs.ps gromacs.pdf

full:		man all

dvi:		gromacs.dvi

#
# make a booklet, i.e. 4 pages onto one double-sided page.
# To get the booklet, rearrange the pages according to page numbering 
# and fold in the middle
#
booklet.ps:	gromacs.ps
		psbook $^ | psnup -2 > ! $@

gromacs.tex:	$(TEXFILES)

gromacs.aux:	gromacs.tex $(AUXFILES)
		$(TEX) gromacs

bib+idx:	gromacs.tex
		$(TEX) gromacs
		$(BIB) gromacs
		$(IDX) gromacs
		./subindex gromacs.ind > gromacs.sind
		mv gromacs.sind gromacs.ind

gromacs.dvi:	bib+idx		gromacs.aux

gromacs.ps:	gromacs.dvi
		dvips -M -o $@ $^

letter.ps:	gromacs.dvi
		dvips -M -t Letter -O 0cm,-0.9cm -o $@ $^

%.aux:		%.tex

prog:		mdp_opt.tex proglist.tex

man:		./mkman

files.tex:	
		$(LOCAL)/prfn; $(RM) files.html; ./mkfiles

progman.tex:	
		$(TOUCH) progman.tex

mdp_opt.tex:	./mkmdp
		./mkmdp $(GMXHOME)

proglist.tex:	$(LOCAL)/mkonline $(LOCAL)/programs.txt
		cd $(LOCAL) ; ./mkonline $(GMXHOME)

man:		
		mkman

copyrgt:
		$(COPYRGT) *.tex

clean:
		$(RM) *.log *.lof *.lot *.bbl *.blg *.toc *.dvi *.aux *.ps *~ #*# *.idx *.ilg *.ind 
		$(RM) progman.tex
		$(RMDIR) progman
