
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

if THREAD_PARALLEL
    THREAD_MPI_LIBOBJS = thread_mpi/libthread_mpi.la 
    THREAD_MPI_DIR     = thread_mpi
endif


SUBDIRS = nonbonded selection statistics trajana $(THREAD_MPI_DIR) \
	  $(BLAS_DIR) $(LAPACK_DIR)

AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libgmx@LIBSUFFIX@.la

#
# Use a utility library for all the nonbonded kernels and wrapper routines.
#
libgmx@LIBSUFFIX@_la_LIBADD =       nonbonded/libnonbonded.la         \
				    selection/libselection.la \
				    statistics/libstatistics.la \
				    trajana/libtrajana.la \
				    $(THREAD_MPI_LIBOBJS) \
	                            $(BLAS_LIBOBJS) $(LAPACK_LIBOBJS)

libgmx@LIBSUFFIX@_la_DEPENDENCIES = nonbonded/libnonbonded.la         \
				    selection/libselection.la \
				    statistics/libstatistics.la \
				    trajana/libtrajana.la \
				    $(THREAD_MPI_LIBOBJS) \
				    $(BLAS_LIBOBJS) $(LAPACK_LIBOBJS)

#	
#
libgmx@LIBSUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ @DLOPEN_LIBS@ $(PTHREAD_LIBS)

libgmx@LIBSUFFIX@_la_SOURCES = \
	3dview.c	atomprop.c	bondfree.c	\
	calcgrid.c	calch.c		chargegroup.c	checkpoint.c	\
	confio.c	copyrite.c	disre.c		do_fit.c	\
	enxio.c		ewald_util.c	ffscanf.c	\
	filenm.c	futil.c		gbutil.c	gmx_fatal.c	\
	gmx_sort.c	gmxcpp.c \
	gmxfio.c	ifunc.c		index.c		inputrec.c	\
	cinvsqrtdata.c	\
	invblock.c	macros.c	orires.c	sparsematrix.c  \
	main.c		maths.c		matio.c		mshift.c	\
	mtop_util.c	mtxio.c		mvdata.c	names.c		\
	network.c	nrama.c		nrjac.c		nrnb.c		\
	pargs.c		pbc.c		pdbio.c		princ.c		\
	rando.c		random.c	gmx_random.c	rbin.c		\
	readinp.c	replace.c	rmpbc.c		shift_util.c	\
	sortwater.c	smalloc.c	statutil.c	sfactor.c	\
	strdb.c		string2.c	symtab.c	\
	topsort.c	tpxio.c		\
	trnio.c		trxio.c		txtdump.c	typedefs.c	\
	viewit.c	warninp.c	\
	wgms.c		wman.c		writeps.c	\
	xdrd.c		xtcio.c		xvgr.c   	replace.h	\
	libxdrf.c	gmx_arpack.c	gmx_matrix.c		\
	dihres.c	gmx_random_gausstable.h	gmxfio_int.h\
	tcontrol.c	splitter.c	gmx_cyclecounter.c		\
	gmx_system_xdr.c md5.c vmdio.c vmddlopen.c	sighandler.c 	\
	oenv.c		gmxfio_rw.c	gmxfio_asc.c	gmxfio_bin.c	\
	gmxfio_xdr.c

pkgconfigdir = ${libdir}/pkgconfig
pkgconfig_DATA = libgmx@LIBSUFFIX@.pc

# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#* innerc.c innerf.f mkinl

if USE_VERSION_H
#version.c contains generated git version information
libgmx@LIBSUFFIX@_la_SOURCES += version.c
CLEANFILES += version.c
#The empty target FORCE forces make to run the command every time. But
#version.c is only changed if the version actually has changed, and hence
#rebuilds are only triggered when they are needed.
version.c: FORCE
	$(top_srcdir)/src/gmxlib/genversion.sh @VERSION@ $(top_srcdir)
FORCE:
endif

EXTRA_DIST = version.h libgmx.pc.cmakein
