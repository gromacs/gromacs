## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libgmxana@LIBSUFFIX@.la

libgmxana@LIBSUFFIX@_la_LIBADD         = 
libgmxana@LIBSUFFIX@_la_DEPENDENCIES   = 
libgmxana@LIBSUFFIX@_la_LDFLAGS        = -version-info @SHARED_VERSION_INFO@


libgmxana@LIBSUFFIX@_la_SOURCES = \
	autocorr.c 	expfit.c 	polynomials.c 	levenmar.c	\
	anadih.c 	pp2shift.c 	pp2shift.h 	dlist.c		\
	eigio.c		lsq.c		cmat.c 		cmat.h		\
	eigensolver.c   eigensolver.h	nsc.c 		nsc.h		\
	hxprops.c 	hxprops.h 	fitahx.c 	fitahx.h	\
	gmx_analyze.c	gmx_anaeig.c	gmx_angle.c	gmx_bond.c	\
	gmx_bundle.c	gmx_chi.c	gmx_cluster.c	gmx_confrms.c	\
	gmx_covar.c	gmx_density.c	gmx_densmap.c   gmx_dih.c	\
	gmx_dielectric.c	gmx_kinetics.c \
	gmx_dipoles.c	gmx_disre.c	gmx_dist.c	gmx_dyndom.c	\
	gmx_enemat.c	gmx_energy.c	gmx_lie.c	gmx_filter.c	\
	gmx_gyrate.c	gmx_h2order.c	gmx_hbond.c	gmx_helix.c	\
	gmx_mindist.c	gmx_msd.c	gmx_morph.c	gmx_nmeig.c	\
	gmx_nmens.c	gmx_order.c	gmx_potential.c	gmx_rama.c	\
	gmx_rdf.c	gmx_rms.c	gmx_rmsdist.c	gmx_rmsf.c	\
	gmx_rotacf.c	gmx_saltbr.c	gmx_sas.c	gmx_sdf.c	\
	gmx_sgangle.c	gmx_sorient.c 	gmx_spatial.c	gmx_tcaf.c	\
	gmx_traj.c	gmx_velacc.c	\
	gmx_clustsize.c	gmx_mdmat.c	gmx_wham.c	eigio.h		\
	correl.c	correl.h	gmx_sham.c	gmx_nmtraj.c	\
	gmx_trjconv.c	gmx_trjcat.c	gmx_trjorder.c	gmx_xpm2ps.c	\
	gmx_editconf.c	gmx_genbox.c	gmx_genion.c	gmx_genconf.c	\
	gmx_genpr.c	gmx_eneconv.c	gmx_wheel.c	\
	addconf.c 	addconf.h	\
	calcpot.c 	calcpot.h 	edittop.c

bin_PROGRAMS = \
	average		do_dssp		editconf	eneconv		\
	genbox		genconf		genpr		g_nmtraj	\
	make_ndx	mk_angndx	trjcat		trjconv        	\
	trjorder	wheel		xpm2ps		genion		\
	anadock		make_edi	\
	g_analyze   	g_anaeig    	g_angle     	g_bond      	\
	g_bundle    	g_chi       	g_cluster   	g_confrms   	\
	g_covar     	g_density   	g_densmap       g_dih       	\
	g_dielectric 	\
	g_dipoles   	g_disre     	g_dist      	g_dyndom    	\
	g_enemat    	g_energy    	g_lie       	g_filter    	\
	g_gyrate    	g_h2order   	g_hbond     	g_helix     	\
	g_mindist   	g_msd       	g_morph     	g_nmeig     	\
	g_nmens     	g_order     	g_potential 	g_rama      	\
	g_rdf       	g_rms       	g_rmsdist   	g_rmsf      	\
	g_rotacf    	g_saltbr    	g_sas       	g_sgangle   	\
	g_sham 		g_sorient   	g_sdf		g_spatial	\
	g_tcaf      	g_traj      	g_velacc    	\
	g_clustsize 	g_mdmat     	g_wham		g_kinetics


LDADD = $(lib_LTLIBRARIES) ../mdlib/libmd@LIBSUFFIX@.la \
	../gmxlib/libgmx@LIBSUFFIX@.la


# link the mpi library to non-mpi names if the latter are not present
install-exec-hook:
	libname="libgmxana@LIBSUFFIX@"; \
	nompi="`echo $$libname | sed -e 's,_mpi,,'`"; \
	libdir="$(libdir)"; \
	if echo $$libname | grep mpi >/dev/null ; then \
	  (cd $$libdir && test -e $$libname.a -a ! -e $$nompi.a && $(LN_S) $$libname.a $$nompi.a ; exit 0); \
	  (cd $$libdir && test -e $$libname.so -a ! -e $$nompi.so && $(LN_S) $$libname.so $$nompi.so ; exit 0); \
	fi; 

CLEANFILES   = *.la *~ \\\#*


