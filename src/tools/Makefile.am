## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#
include $(srcdir)/../Makefile.inc

# This is necessary for VPATH to work, 
# but it can probably be done in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include
LDFLAGS = @LDFLAGS@ -L${top_builddir}/src/gmxlib -L${top_builddir}/src/mdlib 
LDADD = ../gmxlib/libgmxXXX_SUFFIX_XXX.la ../mdlib/libmdXXX_SUFFIX_XXX.la

bin_PROGRAMS = \
	averageXXX_SUFFIX_XXX		do_dsspXXX_SUFFIX_XXX		\
	editconfXXX_SUFFIX_XXX		eneconvXXX_SUFFIX_XXX		\
	g_analyzeXXX_SUFFIX_XXX		g_anaeigXXX_SUFFIX_XXX		\
	g_angleXXX_SUFFIX_XXX		g_bondXXX_SUFFIX_XXX		\
	g_bundleXXX_SUFFIX_XXX		g_chiXXX_SUFFIX_XXX		\
	g_clusterXXX_SUFFIX_XXX		g_confrmsXXX_SUFFIX_XXX		\
	g_covarXXX_SUFFIX_XXX		g_densityXXX_SUFFIX_XXX		\
	g_dielectricXXX_SUFFIX_XXX	g_dihXXX_SUFFIX_XXX		\
	g_dipolesXXX_SUFFIX_XXX		g_disreXXX_SUFFIX_XXX		\
	g_distXXX_SUFFIX_XXX		g_dyndomXXX_SUFFIX_XXX		\
	g_enematXXX_SUFFIX_XXX		g_energyXXX_SUFFIX_XXX		\
	g_lieXXX_SUFFIX_XXX						\
	g_gyrateXXX_SUFFIX_XXX		g_h2orderXXX_SUFFIX_XXX		\
	g_hbondXXX_SUFFIX_XXX		g_helixXXX_SUFFIX_XXX		\
	g_mdmatXXX_SUFFIX_XXX		g_mindistXXX_SUFFIX_XXX		\
	g_msdXXX_SUFFIX_XXX		g_morphXXX_SUFFIX_XXX		\
	g_nmeigXXX_SUFFIX_XXX		g_nmensXXX_SUFFIX_XXX		\
	g_orderXXX_SUFFIX_XXX		g_potentialXXX_SUFFIX_XXX	\
	g_ramaXXX_SUFFIX_XXX						\
	g_rdfXXX_SUFFIX_XXX		g_rmsXXX_SUFFIX_XXX		\
	g_rmsdistXXX_SUFFIX_XXX		g_rmsfXXX_SUFFIX_XXX		\
	g_rotacfXXX_SUFFIX_XXX		g_saltbrXXX_SUFFIX_XXX		\
	g_sasXXX_SUFFIX_XXX		g_sgangleXXX_SUFFIX_XXX		\
	g_sorientXXX_SUFFIX_XXX 	g_tcafXXX_SUFFIX_XXX		\
	g_trajXXX_SUFFIX_XXX		g_velaccXXX_SUFFIX_XXX		\
	genboxXXX_SUFFIX_XXX		genconfXXX_SUFFIX_XXX		\
	genionXXX_SUFFIX_XXX		genprXXX_SUFFIX_XXX		\
	make_ndxXXX_SUFFIX_XXX		mk_angndxXXX_SUFFIX_XXX		\
	trjcatXXX_SUFFIX_XXX		trjconvXXX_SUFFIX_XXX        	\
	trjorderXXX_SUFFIX_XXX		wheelXXX_SUFFIX_XXX		\
	xpm2psXXX_SUFFIX_XXX


AUTOCORR_SRC         = 	autocorr.c  	expfit.c  	polynomials.c  	\
			levenmar.c

averageXXX_SUFFIX_XXX_SOURCES      =    average.c

do_dsspXXX_SUFFIX_XXX_SOURCES      =    do_dssp.c

editconfXXX_SUFFIX_XXX_SOURCES     =    editconf.c

eneconvXXX_SUFFIX_XXX_SOURCES      =    eneconv.c

g_analyzeXXX_SUFFIX_XXX_SOURCES    = 	g_analyze.c 	$(AUTOCORR_SRC)

g_anaeigXXX_SUFFIX_XXX_SOURCES     = 	g_anaeig.c	eigio.c

g_angleXXX_SUFFIX_XXX_SOURCES      = 	g_angle.c  	anadih.c 	$(AUTOCORR_SRC)

g_bondXXX_SUFFIX_XXX_SOURCES       = 	g_bond.c	lsq.c

g_bundleXXX_SUFFIX_XXX_SOURCES     = 	g_bundle.c

g_chiXXX_SUFFIX_XXX_SOURCES        = 	g_chi.c  	anadih.c  	pp2shift.c  	\
			pp2shift.h	dlist.c  	$(AUTOCORR_SRC)

g_clusterXXX_SUFFIX_XXX_SOURCES    = 	g_cluster.c	ql77.c		ql77.h		\
			cmat.c

g_confrmsXXX_SUFFIX_XXX_SOURCES    =    g_confrms.c

g_covarXXX_SUFFIX_XXX_SOURCES      = 	g_covar.c	eigio.c		ql77.c		\
			ql77.h

g_densityXXX_SUFFIX_XXX_SOURCES    =    g_density.c

g_dielectricXXX_SUFFIX_XXX_SOURCES = 	g_dielectric.c	$(AUTOCORR_SRC)

g_dihXXX_SUFFIX_XXX_SOURCES        =    g_dih.c

g_dipolesXXX_SUFFIX_XXX_SOURCES    = 	g_dipoles.c  	lsq.c  	  	$(AUTOCORR_SRC)

g_disreXXX_SUFFIX_XXX_SOURCES      =    g_disre.c

g_distXXX_SUFFIX_XXX_SOURCES       =    g_dist.c

g_dyndomXXX_SUFFIX_XXX_SOURCES     =    g_dyndom.c

g_enematXXX_SUFFIX_XXX_SOURCES     =    g_enemat.c

g_energyXXX_SUFFIX_XXX_SOURCES     = 	g_energy.c  	$(AUTOCORR_SRC)

g_lieXXX_SUFFIX_XXX_SOURCES        = 	g_lie.c

g_gyrateXXX_SUFFIX_XXX_SOURCES     =    g_gyrate.c

g_h2orderXXX_SUFFIX_XXX_SOURCES    =    g_h2order.c

g_hbondXXX_SUFFIX_XXX_SOURCES      = 	g_hbond.c	$(AUTOCORR_SRC)

g_helixXXX_SUFFIX_XXX_SOURCES      = 	g_helix.c	hxprops.c	hxprops.h	\
			fitahx.c	fitahx.h

g_mdmatXXX_SUFFIX_XXX_SOURCES      =    g_mdmat.c

g_mindistXXX_SUFFIX_XXX_SOURCES    =    g_mindist.c

g_msdXXX_SUFFIX_XXX_SOURCES        = 	g_msd.c		lsq.c

g_morphXXX_SUFFIX_XXX_SOURCES      = 	g_morph.c

g_nmeigXXX_SUFFIX_XXX_SOURCES      = 	g_nmeig.c	eigio.c		eigio.h		\
			ql77.c		ql77.h

g_nmensXXX_SUFFIX_XXX_SOURCES      = 	g_nmens.c	eigio.c	

g_orderXXX_SUFFIX_XXX_SOURCES      = 	g_order.c

g_potentialXXX_SUFFIX_XXX_SOURCES  = 	g_potential.c

g_ramaXXX_SUFFIX_XXX_SOURCES       = 	g_rama.c

g_rdfXXX_SUFFIX_XXX_SOURCES        = 	g_rdf.c		$(AUTOCORR_SRC)

g_rmsXXX_SUFFIX_XXX_SOURCES        = 	g_rms.c		cmat.c		cmat.h

g_rmsdistXXX_SUFFIX_XXX_SOURCES    = 	g_rmsdist.c

g_rmsfXXX_SUFFIX_XXX_SOURCES       = 	g_rmsf.c	ql77.c		ql77.h

g_rotacfXXX_SUFFIX_XXX_SOURCES     = 	g_rotacf.c	$(AUTOCORR_SRC)

g_saltbrXXX_SUFFIX_XXX_SOURCES     = 	g_saltbr.c	

g_sasXXX_SUFFIX_XXX_SOURCES        = 	g_sas.c		nsc.c		nsc.h

g_sgangleXXX_SUFFIX_XXX_SOURCES    = 	g_sgangle.c	

g_sorientXXX_SUFFIX_XXX_SOURCES    = 	g_sorient.c	

g_tcafXXX_SUFFIX_XXX_SOURCES       = 	g_tcaf.c	$(AUTOCORR_SRC)

g_trajXXX_SUFFIX_XXX_SOURCES       = 	g_traj.c

g_velaccXXX_SUFFIX_XXX_SOURCES     = 	g_velacc.c	$(AUTOCORR_SRC)

genboxXXX_SUFFIX_XXX_SOURCES       = 	genbox.c  	addconf.c	addconf.h

genconfXXX_SUFFIX_XXX_SOURCES      = 	genconf.c  

genionXXX_SUFFIX_XXX_SOURCES       = 	genion.c  	calcpot.c  	calcpot.h	edittop.c

genprXXX_SUFFIX_XXX_SOURCES        = 	genpr.c  

make_ndxXXX_SUFFIX_XXX_SOURCES     = 	make_ndx.c  

mk_angndxXXX_SUFFIX_XXX_SOURCES    = 	mk_angndx.c  

trjcatXXX_SUFFIX_XXX_SOURCES       = 	trjcat.c

trjconvXXX_SUFFIX_XXX_SOURCES      = 	trjconv.c

trjorderXXX_SUFFIX_XXX_SOURCES     = 	trjorder.c

wheelXXX_SUFFIX_XXX_SOURCES        = 	wheel.c

xpm2psXXX_SUFFIX_XXX_SOURCES       = 	xpm2ps.c


# clean things explicitly, since the target names might have changed
CLEANFILES   = 	${bin_PROGRAMS}	${EXTRA_PROGRAMS}	\
		*_mpi *_d *~ \\\#*


