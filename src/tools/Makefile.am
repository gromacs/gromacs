## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la


bin_PROGRAMS = \
	average		do_dssp		editconf	eneconv		\
	g_analyze	g_anaeig	g_angle		g_bond		\
	g_bundle	g_chi		g_cluster	g_confrms	\
	g_covar		g_density	g_dielectric	g_dih		\
	g_dipoles	g_disre		g_dist		g_dyndom	\
	g_enemat	g_energy	g_lie		g_gyrate	\
	g_h2order	g_hbond		g_helix		g_mdmat		\
	g_mindist	g_msd		g_morph		g_nmeig		\
	g_nmens		g_order		g_potential	g_rama		\
	g_rdf		g_rms		g_rmsdist	g_rmsf		\
	g_rotacf	g_saltbr	g_sas		g_sgangle	\
	g_sorient 	g_tcaf		g_traj		g_velacc	\
	g_clustsize	genbox		genconf		genpr		\
	make_ndx	mk_angndx	trjcat		trjconv        	\
	trjorder	wheel		xpm2ps		genion


AUTOCORR_SRC         = 	autocorr.c expfit.c polynomials.c levenmar.c

# Note: you don't have to list sources for "prog" if it is the single file prog.c

g_analyze_SOURCES    = 	g_analyze.c $(AUTOCORR_SRC)
g_anaeig_SOURCES     = 	g_anaeig.c eigio.c
g_angle_SOURCES      = 	g_angle.c anadih.c $(AUTOCORR_SRC)
g_bond_SOURCES       = 	g_bond.c lsq.c
g_chi_SOURCES        = 	g_chi.c anadih.c pp2shift.c pp2shift.h dlist.c $(AUTOCORR_SRC)
g_cluster_SOURCES    = 	g_cluster.c cmat.c ql77.c ql77.h	
g_clustsize_SOURCES  = 	g_clustsize.c cmat.c 
g_covar_SOURCES      = 	g_covar.c eigio.c ql77.c ql77.h
g_dielectric_SOURCES = 	g_dielectric.c $(AUTOCORR_SRC)
g_dipoles_SOURCES    = 	g_dipoles.c lsq.c $(AUTOCORR_SRC)
g_energy_SOURCES     = 	g_energy.c $(AUTOCORR_SRC)
g_hbond_SOURCES      = 	g_hbond.c $(AUTOCORR_SRC)
g_helix_SOURCES      = 	g_helix.c hxprops.c hxprops.h fitahx.c fitahx.h
g_msd_SOURCES        = 	g_msd.c lsq.c
g_nmeig_SOURCES      = 	g_nmeig.c eigio.c eigio.h ql77.c ql77.h
g_nmens_SOURCES      = 	g_nmens.c eigio.c	
g_rdf_SOURCES        = 	g_rdf.c $(AUTOCORR_SRC)
g_rms_SOURCES        = 	g_rms.c cmat.c cmat.h
g_rmsf_SOURCES       = 	g_rmsf.c ql77.c ql77.h
g_rotacf_SOURCES     = 	g_rotacf.c $(AUTOCORR_SRC)
g_sas_SOURCES        = 	g_sas.c nsc.c nsc.h
g_tcaf_SOURCES       = 	g_tcaf.c $(AUTOCORR_SRC)
g_velacc_SOURCES     = 	g_velacc.c $(AUTOCORR_SRC)
genbox_SOURCES       = 	genbox.c addconf.c addconf.h
genion_SOURCES       = 	genion.c calcpot.c calcpot.h edittop.c

CLEANFILES   = 	*~ \\\#*


