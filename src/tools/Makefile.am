## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libgmxana@LIBSUFFIX@.la

pkgconfigdir = ${libdir}/pkgconfig
pkgconfig_DATA = libgmxana@LIBSUFFIX@.pc

libgmxana@LIBSUFFIX@_la_LIBADD         = ../gmxlib/libgmx@LIBSUFFIX@.la	../mdlib/libmd@LIBSUFFIX@.la	
libgmxana@LIBSUFFIX@_la_DEPENDENCIES   = ../gmxlib/libgmx@LIBSUFFIX@.la	../mdlib/libmd@LIBSUFFIX@.la	
libgmxana@LIBSUFFIX@_la_LDFLAGS        = -no-undefined -version-info @SHARED_VERSION_INFO@ @GSL_LIBS@


libgmxana@LIBSUFFIX@_la_SOURCES = \
	geminate.c	geminate.h					\
	autocorr.c 	expfit.c 	polynomials.c 	levenmar.c	\
	anadih.c 	pp2shift.c 	pp2shift.h 	dlist.c		\
	eigio.c		cmat.c 		cmat.h		\
	eigensolver.c   eigensolver.h	nsc.c 		nsc.h		\
	hxprops.c 	hxprops.h 	fitahx.c 	fitahx.h	\
	gmx_analyze.c	gmx_anaeig.c	gmx_bar.c	\
	gmx_angle.c	gmx_bond.c	\
	gmx_bundle.c	gmx_chi.c	gmx_cluster.c	gmx_confrms.c	\
	gmx_covar.c	gmx_current.c	gmx_dos.c	\
	gmx_density.c	gmx_densmap.c   gmx_dih.c	\
	gmx_dielectric.c	gmx_kinetics.c gmx_spatial.c	\
	gmx_dipoles.c	gmx_disre.c	gmx_dist.c	gmx_dyndom.c	\
	gmx_enemat.c	gmx_energy.c	gmx_lie.c	gmx_filter.c	\
	gmx_gyrate.c	gmx_h2order.c	gmx_hbond.c	gmx_helix.c	\
	gmx_mindist.c	gmx_msd.c	gmx_morph.c	gmx_nmeig.c	\
	gmx_nmens.c	gmx_order.c	gmx_principal.c \
	gmx_polystat.c	gmx_potential.c	gmx_rama.c	\
	gmx_rdf.c	gmx_rms.c	gmx_rmsdist.c	gmx_rmsf.c	\
	gmx_rotacf.c	gmx_rotmat.c	gmx_saltbr.c	gmx_sas.c	\
	gmx_select.c	gmx_pme_error.c	gmx_options.c\
	gmx_sgangle.c	gmx_sorient.c 	gmx_spol.c	gmx_tcaf.c	\
	gmx_traj.c	gmx_velacc.c	gmx_helixorient.c \
	gmx_clustsize.c	gmx_mdmat.c	gmx_wham.c	eigio.h		\
	correl.c	correl.h	gmx_sham.c	gmx_nmtraj.c	\
	gmx_trjconv.c	gmx_trjcat.c	gmx_trjorder.c	gmx_xpm2ps.c	\
	gmx_editconf.c	gmx_genbox.c	gmx_genion.c	gmx_genconf.c	\
	gmx_genpr.c	gmx_eneconv.c	gmx_vanhove.c	gmx_wheel.c	\
	addconf.c 	addconf.h	gmx_tune_pme.c  gmx_membed.c     \
	calcpot.c 	calcpot.h 	edittop.c	\
	dens_filter.c	dens_filter.h	\
	powerspect.c	powerspect.h	interf.h	\
	gmx_densorder.c	gmx_hydorder.c	\
	binsearch.h	binsearch.c

bin_PROGRAMS = \
	do_dssp		editconf	eneconv		\
	genbox		genconf		genrestr	g_nmtraj	\
	make_ndx	mk_angndx	trjcat		trjconv        	\
	trjorder	g_wheel		xpm2ps		genion		\
	g_anadock	make_edi	\
	g_analyze   	g_anaeig    	g_bar		\
	g_angle     	g_bond      	\
	g_bundle    	g_chi       	g_cluster   	g_confrms   	\
	g_covar     	g_current	g_dos \
	g_density   	g_densmap       g_dih       	\
	g_dielectric 	g_helixorient 	g_principal \
	g_dipoles   	g_disre     	g_dist      	g_dyndom    	\
	g_enemat    	g_energy    	g_lie       	g_filter    	\
	g_gyrate    	g_h2order   	g_hbond     	g_helix     	\
	g_mindist   	g_msd       	g_morph     	g_nmeig     	\
	g_nmens     	g_order     	g_kinetics	\
	g_polystat	g_potential 	g_rama      	\
	g_rdf       	g_rms       	g_rmsdist   	g_rmsf      	\
	g_rotacf    	g_rotmat	g_saltbr    	g_sas       	\
	g_select	g_sgangle   	\
	g_sham 		g_sorient   	g_spol		\
	g_spatial	g_pme_error		g_options	\
	g_tcaf      	g_traj      	g_tune_pme   \
	g_vanhove	g_velacc    	g_membed      \
	g_clustsize 	g_mdmat     	g_wham		g_kinetics	\
	g_sigeps	g_densorder	g_hydorder

EXTRA_DIST		=	libgmxana.pc.cmakein

LDADD = $(lib_LTLIBRARIES) ../mdlib/libmd@LIBSUFFIX@.la \
	../gmxlib/libgmx@LIBSUFFIX@.la


CLEANFILES   = *.la *~ \\\#*


