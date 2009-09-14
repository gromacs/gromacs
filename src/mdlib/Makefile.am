## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

libmd@LIBSUFFIX@_la_LIBADD         = ../gmxlib/libgmx@LIBSUFFIX@.la	
libmd@LIBSUFFIX@_la_DEPENDENCIES   = ../gmxlib/libgmx@LIBSUFFIX@.la	
libmd@LIBSUFFIX@_la_LDFLAGS        = -version-info @SHARED_VERSION_INFO@

lib_LTLIBRARIES = libmd@LIBSUFFIX@.la

libmd@LIBSUFFIX@_la_SOURCES = \
	adress.c	adress.h	\
	calcmu.c	calcvir.c	constr.c	\
	coupling.c	\
	domdec.c	domdec_box.c	domdec_con.c	\
	domdec_network.c domdec_setup.c	domdec_top.c	\
	ebin.c		\
	edsam.c		ewald.c         fftgrid.c	\
	force.c  	ghat.c		init.c		\
	mdatom.c	mdebin.c	minimize.c	\
	mvxvf.c		ns.c		nsgrid.c	\
	perf_est.c	genborn.c			\
	genborn_sse2_single.c				\
	genborn_sse2_single.h				\
	genborn_sse2_double.c				\
	genborn_sse2_double.h				\
	genborn_allvsall.c                              \
	genborn_allvsall.h                              \
	genborn_allvsall_sse2_single.c			\
	genborn_allvsall_sse2_single.h			\
	gmx_qhop_parm.c	gmx_qhop_parm.h			\
	gmx_qhop_xml.c	gmx_qhop_xml.h			\
	gmx_qhop_db.c	gmx_qhop_db.h			\
	pme.c		pme_pp.c	pppm.c		\
	partdec.c	pull.c 		pullutil.c	\
	rf_util.c	shakef.c	sim_util.c	\
	shellfc.c	stat.c 		\
	tables.c	tgroup.c	\
	update.c	vcm.c		vsite.c	\
	wall.c		wnblist.c	\
	csettle.c       clincs.c	\
	qmmm.c		gmx_fft.c 	gmx_parallel_3dfft.c	\
	gmx_wallcycle.c \
	qm_gaussian.c	qm_mopac.c	qm_gamess.c		\
	gmx_fft_fftw2.c	gmx_fft_fftw3.c	gmx_fft_fftpack.c	\
	gmx_fft_mkl.c

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la 

EXTRA_PROGRAMS = gmx_qhop_db_test

gmx_qhop_db_test_LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la ../kernel/libgmxpreprocess@LIBSUFFIX@.la

# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#*


















