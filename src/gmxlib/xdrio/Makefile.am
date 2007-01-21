#
## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libxdrio@LIBSUFFIX@.la

#
libmd@LIBSUFFIX@_la_LIBADD         = @XDRIO_LIBOBJS@	
libmd@LIBSUFFIX@_la_DEPENDENCIES   = @XDRIO_LIBOBJS@	
libxdrio@LIBSUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

libxdrio@LIBSUFFIX@_la_SOURCES = \
	xdrfile.c

EXTRA_PROGRAMS = test_xdrio_c test_xdrio_f
test_xdrio_c_SOURCES = test_xdrio_c.c
test_xdrio_f_SOURCES = test_xdrio_f.f
test_xdrio_c_LDADD   = libxdrio@LIBSUFFIX@.la
test_xdrio_f_LDADD   = libxdrio@LIBSUFFIX@.la

# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#*







