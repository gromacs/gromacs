## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

templatedir = ${pkgdatadir}/template

# we don't build this, it's meant for the user to play with.
# Since a Makefile is used for the installation, we call the
# template makefile Template.mak, and rename it to Makefile when installing.
# We add the GMX_DOUBLE define to CPPFLAGS, so users can do without config.h
XLDFLAGS = -L@libdir@ @LDFLAGS@
XLIBS = -lmd@LIBSUFFIX@ -lgmx@LIBSUFFIX@ @FFT_LIBS@ @XML_LIBS@ @GSL_LIBS@ @LIBS@ 
XCC = @CC@
if GMX_DOUBLE
MYCFLAGS = @CFLAGS@ @INCLUDES@ -I@includedir@ -DGMX_DOUBLE
else
MYCFLAGS = @CFLAGS@ @INCLUDES@ -I@includedir@
endif

if GMX_DOUBLE
template_DATA = template.c README Makefile.@host@_double Makefile.pkg
else
template_DATA = template.c README Makefile.@host@ Makefile.pkg
endif
EXTRA_DIST = template.c README Template.mak Makefile.pkg CMakeLists.txt.template

BUILT_SOURCES = gromacs

gromacs:
	$(LN_S) $(top_srcdir)/include gromacs

Makefile.@host@ Makefile.@host@_double: Template.mak Makefile 
	cat $(srcdir)/Template.mak | \
	sed 's,@LDFLAGS\@,$(XLDFLAGS),' | \
	sed 's,@LIBS\@,$(XLIBS),' | \
	sed 's,@CFLAGS\@,$(MYCFLAGS),' | \
	sed 's,@CC\@,$(XCC),' | \
	sed 's,@host\@,@host@,' >$@

CLEANFILES = Makefile.@host@ Makefile.@host@_double *~ \\\#*

template_SOURCES = template.c
noinst_PROGRAMS = template
LDADD = ../../src/mdlib/libmd@LIBSUFFIX@.la ../../src/gmxlib/libgmx@LIBSUFFIX@.la
#template.$(OBJEXT): gromacs

clean-local:
	rm -rf gromacs
