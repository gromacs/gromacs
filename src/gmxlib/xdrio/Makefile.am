#
## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libxdrio@LIBSUFFIX@.la

#
# Use a utility library for all the nonbonded kernels and wrapper routines.
#
libxdrio@LIBSUFFIX@_la_LIBADD =       @XDRIOLIB_LIBOBJS@                  \
			            nonbonded/libnonbonded.la         \
	                            $(BLAS_LIBOBJS) $(LAPACK_LIBOBJS)

libxdrio@LIBSUFFIX@_la_DEPENDENCIES = @XDRIOLIB_LIBOBJS@                  \
			            nonbonded/libnonbonded.la         \
				    $(BLAS_LIBOBJS) $(LAPACK_LIBOBJS)

#	
#
libxdrio@LIBSUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

libxdrio@LIBSUFFIX@_la_SOURCES = \
	xdrfile.o

EXTRA_PROGRAMS = test_xdrio_c test_xdrio_f
test_xdrio_c_SOURCES = test_xdrio_c.c
test_xdrio_f_SOURCES = test_xdrio_f.f

# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#*







