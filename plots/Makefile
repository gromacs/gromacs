# Makefile for figures
MKEPS	= xmgrace -hardcopy -hdevice EPS -printfile 
MKPDF	= xmgrace -hardcopy -hdevice PDF -printfile 

all:	pdf

eps:	f_angle.eps	f_bham.eps	f_bond.eps	f_dih.eps\
	f_dr.eps	f_imps.eps	f_lj.eps	f_morse.eps\
	f_pr.eps	f_rbs.eps	maxwell.eps	msdwater.eps\
	rama.eps	rdfO_O.eps	softcore.eps	vcrf.eps

pdf:	f_angle.pdf	f_bham.pdf	f_bond.pdf	f_dih.pdf\
	f_dr.pdf	f_imps.pdf	f_lj.pdf	f_morse.pdf\
	f_pr.pdf	f_rbs.pdf	maxwell.pdf	msdwater.pdf\
	rama.pdf	rdfO_O.pdf	softcore.pdf	vcrf.pdf

%.eps:	%.xvg
	$(MKEPS) $@ $<
	chmod +x ./myvi
	./myvi $@
%.pdf:	%.xvg
	$(MKPDF) $@ $<
