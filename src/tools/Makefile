#
#	@(#) Makefile 1.44 11/3/97
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990, 1991, 1992, Groningen University
#
# 	Makefile for tools
#
include		$(GMXHOME)/src/makef/Makefile.def
#
LIB		=
POBJ		= g_hbond.o dah.o hbond.o list.o parse.o
PROGRAMS	= g_rdf		g_msd		g_dist		\
		  g_hbond	g_lifetime	merge_freq	\
		  genbox 	genconf 	genion		\
		  make_ndx	genpr 		\
		  trjconv       mk_angndx	eneconv		\
		  editconf	protonate	gendr		\
		  g_rms		g_rama		g_disre		\
		  g_gyrate	g_mindist	g_sgangle	\
		  g_angle	g_dih		\
		  g_bond	g_dipoles			\
		  g_mdmat	g_energy	g_sas		\
		  do_dssp	sas2mat         g_rmsdist	\
		  g_density     g_potential     g_order         \
                  g_saltbr	g_rmsf		g_coord         \
                  g_h2order     g_com		g_rdens		\
		  g_velacc	g_rotacf	g_enemat	\
		  g_covar       g_anaeig	\
		  g_nmeig	g_chi		g_dielectric	\
		  g_helix	wheel	        		\
                  g_confrms     xpm2ps 

ACOBJ	=	autocorr.o	expfit.o	polynomials.o	\
		levenmar.o
#
include		$(STDTARGET)

genbox:		genbox.o	addconf.o
		$(LD) -o  $@ $^ $(SYSLIBS)

genconf:	genconf.o
		$(LD) -o  $@ $^ $(SYSLIBS)

gendr:		gendr.o
		$(LD) -o  $@ $^ $(SYSLIBS)

protonate:	protonate.o
		$(LD) -o  $@ $^ $(SYSLIBS)

eneconv:	eneconv.o
		$(LD) -o  $@ $^ $(SYSLIBS)

genion:		genion.o
		$(LD) -o  $@ $^ $(SYSLIBS)

gen_xh:		gen_xh.o
		$(LD) -o  $@ $^ $(SYSLIBS)

centerx:	centerx.o
		$(LD) -o  $@ $^ $(SYSLIBS)

patchmagic:	patchmagic.o 
		$(LD) -o $@ $^ $(SYSLIBS)

genpr:		genpr.o $(OBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

econv:		econv.o
		$(LD) -o $@ $^ $(SYSLIBS)

stat2asc:	stat2asc.o $(OBJ)
		$(LD) -o $@ $^  $(SYSLIBS)

make_ndx:	make_ndx.o
		$(LD) -o  $@ $^ $(SYSLIBS)

mk_angndx:	mk_angndx.o
		$(LD) -o  $@ $^ $(SYSLIBS)

g_relax:	g_relax.o	$(ACOBJ)
		$(LD) -o  $@ $^ $(SYSLIBS)

truncstat:	truncstat.o
		$(LD) -o $@ $^ $(SYSLIBS)

average:	average.o
		$(LD) -o $@ $^ $(SYSLIBS)

trjconv:	trjconv.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_rmsdist:	g_rmsdist.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_diffuse:	g_diffuse.o
		$(LD) -o $@ $^ $(SYSLIBS)

editconf:	editconf.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_energy:	g_energy.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_mdmat:	g_mdmat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_disre:	g_disre.o
		$(LD) -o $@ $^ $(SYSLIBS) $(FFTLIBS) -lm

g_gyrate:	g_gyrate.o  
		$(LD) -o $@ $^ $(SYSLIBS)

g_chi:		g_chi.o		anadih.o	pp2shift.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_velacc:	g_velacc.o	$(ACOBJ)
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_angle:	g_angle.o  	anadih.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_dipoles:	g_dipoles.o  	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_rotacf:	g_rotacf.o 	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

testacf:	testacf.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_dielectric:	g_dielectric.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_bond:		g_bond.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_com:		g_com.o  
		$(LD) -o $@ $^ $(SYSLIBS)

g_rms:	 	g_rms.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_rmsf:	 	g_rmsf.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_mindist:	g_mindist.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_dih:		g_dih.o	
		$(LD) -o $@ $^ $(SYSLIBS)

do_dssp:	do_dssp.o
		$(LD) -o $@ $^ $(SYSLIBS)

xpm2ps:		xpm2ps.o
		$(LD) -o $@ $^ $(SYSLIBS)

sas2mat:	sas2mat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_rama:		g_rama.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_sgangle:	g_sgangle.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_order:	g_order.o
	 	$(LD) -o $@ $^ $(SYSLIBS)

g_h2order:      g_h2order.o
	        $(LD) -o $@ $^ $(SYSLIBS)

g_potential:	g_potential.o
	 	$(LD) -o $@ $^ $(SYSLIBS)

g_density:	g_density.o
	 	$(LD) -o $@ $^ $(SYSLIBS)

g_coord:	g_coord.o
	        $(LD) -o $@ $^ $(SYSLIBS)

g_saltbr:	g_saltbr.o
		$(LD) -o  $@ $^ $(SYSLIBS)

g_diff:         g_diff.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_rdens:        g_rdens.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_nmeig.o:	g_nmeig.c

ifeq ($(USEF77),yes)
g_nmeig:        g_nmeig.o
	        $(LD) -o  $@ $^ $(SYSLIBS)
else
ql77.o:		ql77.c
g_nmeig:        g_nmeig.o ql77.o
	        $(LD) -o  $@ $^ $(SYSLIBS)
endif
g_nmeig.o:	g_nmeig.c

ifeq ($(USEF77),yes)
g_covar:          g_covar.o
	        $(LD) -o  $@ $^ $(SYSLIBS)
else
ql77.o:		ql77.c
g_covar:          g_covar.o ql77.o
	        $(LD) -o  $@ $^ $(SYSLIBS)
endif

g_anaeig:       g_anaeig.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_dist:         g_dist.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

ipgro:          ipgro.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

my_rdf:         my_rdf.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_confrms:      g_confrms.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

#
# Helix stuff
#
g_helix:	g_helix.o	hxprops.o	fitahx.o
		$(LD) -o $@ $^ $(SYSLIBS)

proptim:	proptim.o	hxprops.o	readev.o	\
		mcprop.o	orise.o		pinput.o	\
		recomb.o
		$(LD) -o $@ $^ $(SYSLIBS)

wheel:		wheel.o
		$(LD) -o $@ $^ $(SYSLIBS)

mkice:		mkice.o
		$(LD) -o $@ $^ $(SYSLIBS)

fcount:		fcount.o
		$(LD) -o $@ $^ $(SYSLIBS)
waterins:	waterins.o
		$(LD) -o $@ $^ $(SYSLIBS)

anaeig:		anaeig.o	readev.o
	        $(LD) -o $@ $^ $(SYSLIBS)

conveig:	conveig.o
		$(F77) -o $@ $^ -lm

g_rdf.o:	g_rdf.cc	do_rdf.h
do_rdf.o:	do_rdf.cc	do_rdf.h
rdfang.o:	rdfang.cc	do_rdf.h
rdfn.o:		rdfn.cc		do_rdf.h
g_rdf:		g_rdf.o 	do_rdf.o	rdfn.o	rdfang.o
		$(CCLD) -o $@ $^ $(SYSLIBS)

corrutil.o:	corrutil.h	corrutil.cc
g_msd.o:	corrutil.h	g_msd.cc

g_msd:		corrutil.o	g_msd.o
		$(CCLD) -o $@ $^ $(SYSLIBS)

g_hbond:	g_hbond.o hbond.o list.o dah.o parse.o 
		$(CCLD) -g -o $@ $^ $(SYSLIBS)


g_lifetime:	g_lifetime.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_quadrup:	g_quadrup.o
		$(LD) -o $@ $^ $(SYSLIBS)

merge_freq:	merge_freq.o
		$(LD) -o $@ $^ $(SYSLIBS)
g_enemat:	g_enemat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_hbond.o: g_hbond.cc
g_hbond.o: g_hbond.h
g_hbond.o: parse.h
g_hbond.o: hbond.h
g_hbond.o: list.h
g_hbond.o: dah.h
g_hbond.o: water.h

g_sas:		g_sas.o	nsc.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_multipoles:	g_multipoles.o
		$(LD) -o $@ $^ $(SYSLIBS)

html_luck:	html_luck.o
		$(LD) -Bstatic -o $@ $^ $(SYSLIBS)
include	Make.dep
