#!gmake
#
#	@(#) Makefile 1.44 11/3/97
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990, 1991, 1992, Groningen University
#
# 	Makefile for tools

ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/makef/Makefile.def

LIB		=
ifneq ($(PARALLEL_ONLY_MDRUN),yes)
PROGRAMS	= average	do_dssp		editconf	\
		eneconv		g_analyze	g_anaeig	\
		g_angle		g_bond		g_chi		\
		g_cluster	\
		g_confrms	g_covar		g_density	\
		g_dielectric	g_dih		g_dipoles	\
		g_disre		g_dist		g_enemat	\
		g_energy	g_gyrate	g_h2order	\
		g_hbond		g_helix		g_mdmat		\
		g_mindist	g_msd		\
		g_nmeig		g_nmens		\
		g_order		g_potential	g_rama		\
		g_rdens		g_rdf		g_rms		\
		g_rmsdist	g_rmsf		g_rotacf	\
		g_saltbr	g_sas		g_sgangle	\
		g_traj		\
		g_velacc	genbox		genconf		\
		gendr		genion		genpr		\
		make_ndx	mk_angndx	protonate	\
		trjconv		          trjorder	\
		wheel		xpm2ps
endif

ACOBJ	=	autocorr.o	expfit.o	polynomials.o	\
		levenmar.o

include		$(STDTARGET)

movebin:
	mv $(EXECUTABLES) $(BINDIR)

disco:		disco.o dlist.o dtools.o correct.o
		$(LD) -o  $@ $^ $(SYSLIBS)

cdist:		cdist.o bondlist.o smooth.o residues.o virtual.o
		$(LD) -o  $@ $^ $(SYSLIBS)

genbox:		genbox.o	addconf.o	
		$(LD) $(RUNLDOPT) -o $@ $^ $(RUNSYSLIBS)

genconf:	genconf.o
		$(LD) -o  $@ $^ $(SYSLIBS)

gendr:		gendr.o
		$(LD) -o  $@ $^ $(SYSLIBS)

protonate:	protonate.o
		$(LD) -o  $@ $^ $(SYSLIBS)

eneconv:	eneconv.o
		$(LD) -o  $@ $^ $(SYSLIBS)

genion:		genion.o calcpot.o edittop.o
		$(LD) $(RUNLDOPT) -o $@ $^ $(RUNSYSLIBS)

genpr:		genpr.o $(OBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

make_ndx:	make_ndx.o
		$(LD) -o  $@ $^ $(SYSLIBS)

mk_angndx:	mk_angndx.o
		$(LD) -o  $@ $^ $(SYSLIBS)

g_relax:	g_relax.o	$(ACOBJ)
		$(LD) -o  $@ $^ $(SYSLIBS)

trjcat:		trjcat.o
		$(LD) -o $@ $^ $(SYSLIBS)

trjconv:	trjconv.o
		$(LD) -o $@ $^ $(SYSLIBS)

trjorder:	trjorder.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_rmsdist:	g_rmsdist.o
		$(LD) -o $@ $^ $(SYSLIBS)

editconf:	editconf.o	addconf.o
		$(LD) $(RUNLDOPT) -o $@ $^ $(RUNSYSLIBS)

g_energy:	g_energy.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_mdmat:	g_mdmat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_disre:	g_disre.o
		$(LD) -o $@ $^ $(RUNSYSLIBS)

g_gyrate:	g_gyrate.o  
		$(LD) -o $@ $^ $(SYSLIBS)

g_chi:		g_chi.o		anadih.o pp2shift.o	dlist.o $(ACOBJ)
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

g_bond:		g_bond.o 	lsq.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_traj:		g_traj.o  
		$(LD) -o $@ $^ $(SYSLIBS)

g_rms:	 	g_rms.o cmat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_rmsf:	 	g_rmsf.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_diffuse:	g_diffuse.o
		$(LD) -o $@ $^ $(SYSLIBS)

average:	average.o
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

g_saltbr:	g_saltbr.o
		$(LD) -o  $@ $^ $(SYSLIBS)

g_rdens:        g_rdens.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_nmeig:        g_nmeig.o eigio.o ql77.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_nmens:	g_nmens.o eigio.o
		$(LD) -o  $@ $^ $(SYSLIBS)

g_covar:        g_covar.o eigio.o ql77.o
		$(LD) -o  $@ $^ $(SYSLIBS)

g_anaeig:       g_anaeig.o eigio.o
	        $(LD) -o  $@ $^ $(SYSLIBS)

g_dist:         g_dist.o
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

g_rdf:		g_rdf.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_msd:		g_msd.o	lsq.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_hbond:	g_hbond.o	$(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

g_enemat:	g_enemat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_sas:		g_sas.o	nsc.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_cluster:	g_cluster.o ql77.o cmat.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_multipoles:	g_multipoles.o
		$(LD) -o $@ $^ $(SYSLIBS)

g_analyze:	g_analyze.o $(ACOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

include	Make.dep
endif
