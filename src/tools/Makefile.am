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
	ql77.c 		ql77.h		nsc.c 		nsc.h		\
	hxprops.c 	hxprops.h 	fitahx.c 	fitahx.h	\
	g_analyze.c	g_anaeig.c	g_angle.c	g_bond.c	\
	g_bundle.c	g_chi.c		g_cluster.c	g_confrms.c	\
	g_covar.c	g_density.c	g_dielectric.c	g_dih.c		\
	g_dipoles.c	g_disre.c	g_dist.c	g_dyndom.c	\
	g_enemat.c	g_energy.c	g_lie.c		g_filter.c	\
	g_gyrate.c	g_h2order.c	g_hbond.c	g_helix.c	\
	g_mindist.c	g_msd.c		g_morph.c	g_nmeig.c	\
	g_nmens.c	g_order.c	g_potential.c	g_rama.c	\
	g_rdf.c		g_rms.c		g_rmsdist.c	g_rmsf.c	\
	g_rotacf.c	g_saltbr.c	g_sas.c		g_sgangle.c	\
	g_sorient.c 	g_tcaf.c	g_traj.c	g_velacc.c	\
	g_clustsize.c	g_mdmat.c	g_wham.c	eigio.h


bin_PROGRAMS = \
	average		do_dssp		editconf	eneconv		\
	gmx-ana		genbox		genconf		genpr		\
	make_ndx	mk_angndx	trjcat		trjconv        	\
	trjorder	wheel		xpm2ps		genion		\
	anadock		cdist		disco           

dist_bin_SCRIPTS = \
	g_energy	g_wham		g_velacc	g_traj		\
	g_tcaf		g_sorient	g_sgangle	g_sas		\
	g_saltbr	g_rotacf	g_rmsf		g_rmsdist	\
	g_rms		g_rdf		g_rama		g_potential	\
	g_order		g_nmens		g_nmeig		g_msd		\
	g_morph		g_mindist	g_mdmat		g_lie		\
	g_helix		g_hbond		g_h2order	g_gyrate	\
	g_filter	g_enemat	g_dyndom	g_dist		\
	g_disre		g_dipoles	g_dih		g_dielectric	\
	g_density	g_covar		g_confrms	g_clustsize	\
	g_cluster	g_chi		g_bundle	g_bond		\
	g_angle		g_analyze	g_anaeig	

genbox_SOURCES       = 	genbox.c addconf.c addconf.h
genion_SOURCES       = 	genion.c calcpot.c calcpot.h edittop.c
cdist_SOURCES        =  cdist.c smooth.c bondlist.c residues.c 
disco_SOURCES        =  disco.c correct.c discopar.c dtools.c 
anadock_SOURCES	     =  anadock.c  

LDADD = $(lib_LTLIBRARIES) ../mdlib/libmd@LIBSUFFIX@.la \
	../gmxlib/libgmx@LIBSUFFIX@.la

# Create scripts. This is never done in a normal build - the script files
# should be in CVS, so we can do it in the source directory.
scripts:	
	./gmx-ana.pl $(dist_bin_SCRIPTS)

# link the mpi library to non-mpi names if the latter are not present
install-exec-hook:
	libname="libgmxana@LIBSUFFIX@"; \
	nompi="`echo $$libname | sed -e 's,_mpi,,'`"; \
	libdir="$(libdir)"; \
	if echo $$libname | grep mpi >/dev/null ; then \
	  (cd $$libdir && test -e $$libname.a -a ! -e $$nompi.a && $(LN_S) $$libname.a $$nompi.a ; exit 0); \
	  (cd $$libdir && test -e $$libname.so -a ! -e $$nompi.so && $(LN_S) $$libname.so $$nompi.so ; exit 0); \
	fi; 
	/bin/cp $(xxx) $(bindir)

CLEANFILES   = *.la *~ \\\#*


