
## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

if !GMX_EXTERNAL_BLAS
  BLAS_DIR       = gmx_blas
  BLAS_LIBOBJS   = gmx_blas/libblas.la       
endif

if !GMX_EXTERNAL_LAPACK
  LAPACK_DIR     = gmx_lapack
  LAPACK_LIBOBJS = gmx_lapack/liblapack.la
endif


SUBDIRS = nonbonded $(BLAS_DIR) $(LAPACK_DIR)

AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libgmx@LIBSUFFIX@.la

#
# Use a utility library for all the nonbonded kernels and wrapper routines.
#
libgmx@LIBSUFFIX@_la_LIBADD =       @GMXLIB_LIBOBJS@                  \
			            nonbonded/libnonbonded.la         \
	                            $(BLAS_LIBOBJS) $(LAPACK_LIBOBJS)

libgmx@LIBSUFFIX@_la_DEPENDENCIES = @GMXLIB_LIBOBJS@                  \
			            nonbonded/libnonbonded.la         \
				    $(BLAS_LIBOBJS) $(LAPACK_LIBOBJS)

#	
#
libgmx@LIBSUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

libgmx@LIBSUFFIX@_la_SOURCES = \
	3dview.c	atomprop.c	block_tx.c	bondfree.c	\
	calcgrid.c	calch.c		\
	confio.c	copyrite.c	disre.c		do_fit.c	\
	enxio.c		ewald_util.c	gmx_fatal.c	ffscanf.c	\
	filenm.c	futil.c		gbutil.c	\
	gmxfio.c	ifunc.c		index.c		cinvsqrtdata.c	\
	invblock.c	macros.c	orires.c	sparsematrix.c  \
	main.c		maths.c		matio.c		mtxio.c		\
	mshift.c	mvdata.c	mvxvf.c		names.c		\
	network.c	nrama.c		nrjac.c		nrnb.c		\
	pargs.c		pbc.c		pdbio.c		princ.c		\
	rando.c		random.c	gmx_random.c	rbin.c		\
	readinp.c	replace.c	rmpbc.c		shift_util.c	\
	sortwater.c	smalloc.c	stat.c		statutil.c	\
	strdb.c		string2.c	symtab.c	tpxio.c		\
	trnio.c		trxio.c		txtdump.c	typedefs.c	\
	viewit.c	wgms.c		wman.c		writeps.c	\
	xtcio.c		xvgr.c   	replace.h	\
	libxdrf.c	xmlio.h		gmx_arpack.c			\
	dihres.c	xmlio.c		gmx_random_gausstable.h		\
	tcontrol.c

EXTRA_libgmx@LIBSUFFIX@_la_SOURCES = \
	mgmx.c		widget.c	widget.h	\
	gmx_system_xdr.c


# link the mpi library to non-mpi names if the latter are not present after install
#install-exec-hook:
#	libname="libgmx@LIBSUFFIX@"; \
#	nompi="`echo $$libname | sed -e 's,_mpi,,'`"; \
#	libdir="$(libdir)"; \
#	if echo $$libname | grep mpi >/dev/null ; then \
#	  (cd $$libdir && test -e $$libname.a -a ! -e $$nompi.a && $(LN_S) $$libname.a $$nompi.a ; exit 0); \
#	  (cd $$libdir && test -e $$libname.so -a ! -e $$nompi.so && $(LN_S) $$libname.so $$nompi.so ; exit 0); \
#	fi; 

# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#* innerc.c innerf.f mkinl 







