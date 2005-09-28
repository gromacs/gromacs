## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

libmd@LIBSUFFIX@_la_LIBADD         = @MDLIB_LIBOBJS@	
libmd@LIBSUFFIX@_la_DEPENDENCIES   = @MDLIB_LIBOBJS@	
libmd@LIBSUFFIX@_la_LDFLAGS        = -version-info @SHARED_VERSION_INFO@

lib_LTLIBRARIES = libmd@LIBSUFFIX@.la

libmd@LIBSUFFIX@_la_SOURCES = \
	calcmu.c	calcvir.c	constr.c	\
	coupling.c	ebin.c		\
	edsam.c		ewald.c         \
	force.c  	ghat.c		init.c		\
	mdatom.c	mdebin.c	minimize.c	\
	ns.c		nsb.c		nsgrid.c	\
	pme.c		pppm.c		fftgrid.c	\
	pull.c 		pullinit.c	\
	pullio.c	pullutil.c	\
	rf_util.c	shakef.c	sim_util.c	\
	splittop.c	tables.c	tgroup.c	\
	update.c	vcm.c		vsite.c	\
	wnblist.c	splittop.h	wnblist.h	\
	pull_internal.h	csettle.c       clincs.c	\
	qmmm.c		gmx_fft.c	gmx_parallel_3dfft.c


EXTRA_libmd@LIBSUFFIX@_la_SOURCES = \
	qm_gaussian.c	qm_mopac.c	qm_gamess.c		\
	gmx_fft_fftw2.c	gmx_fft_fftw3.c	gmx_fft_fftpack.c	\
	gmx_fft_mkl.c	



LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la 

EXTRA_PROGRAMS = compnl
compnl_SOURCES = compnl.c
compnl_LDADD   = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la


# link the mpi library to non-mpi names if the latter are not present
install-exec-hook:
#	libname="libmd@LIBSUFFIX@"; \
#	nompi="`echo $$libname | sed -e 's,_mpi,,'`"; \
#	libdir="$(libdir)"; \
#	if echo $$libname | grep mpi >/dev/null ; then \
#	  (cd $$libdir && test -e $$libname.a -a ! -e $$nompi.a && $(LN_S) $$libname.a $$nompi.a ; exit 0); \
#	  (cd $$libdir && test -e $$libname.so -a ! -e $$nompi.so && $(LN_S) $$libname.so $$nompi.so ; exit 0); \
#	fi; 


# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#*


















