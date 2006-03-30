# Makefile for figures
MKEPS	= xmgrace -hardcopy -hdevice EPS -printfile 
MKPNG	= xmgrace -hardcopy -hdevice PNG -printfile 

all:	eps png

eps:	f_angle.eps	f_bham.eps	f_bond.eps	f_dih.eps\
	f_dr.eps	f_imps.eps	f_lj.eps	f_morse.eps\
	f_pr.eps	f_rbs.eps	maxwell.eps	msdwater.eps\
	rama.eps	rdfO_O.eps	softcore.eps	vcrf.eps

png:	f_angle.png	f_bham.png	f_bond.png	f_dih.png\
	f_dr.png	f_imps.png	f_lj.png	f_morse.png\
	f_pr.png	f_rbs.png	maxwell.png	msdwater.png\
	rama.png	rdfO_O.png	softcore.png	vcrf.png

%.eps:	%.xvg
	$(MKEPS) $@ $<
	chmod +x ./myvi
	./myvi $@
%.png:	%.xvg
	$(MKPNG) $@ $<
