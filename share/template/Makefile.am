## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

templatedir = ${pkgdatadir}/template

# we don't build this, it's meant for the user to play with.
# Since a Makefile is used for the installation, we call the
# template makefile Template.mak, and rename it to Makefile when installing.
# We add the GMX_DOUBLE define to CPPFLAGS, so users can do without config.h
LDFLAGS = -L@libdir@ @LDFLAGS@
LIBS = -lmd@LIBSUFFIX@ -lgmx@LIBSUFFIX@ @LIBS@ 
CC = @CC@
if GMX_DOUBLE
MYCFLAGS = @CFLAGS@ @INCLUDES@ -I@includedir@/gromacs -DGMX_DOUBLE
else
MYCFLAGS = @CFLAGS@ @INCLUDES@ -I@includedir@/gromacs
endif

if GMX_DOUBLE
template_DATA = template.c README Makefile.@host@_double
else
template_DATA = template.c README Makefile.@host@
endif
EXTRA_DIST = template.c README Template.mak

Makefile.@host@ Makefile.@host@_double: Template.mak Makefile 
	cat $(srcdir)/Template.mak | \
	sed 's,@LDFLAGS\@,$(LDFLAGS),' | \
	sed 's,@LIBS\@,$(LIBS),' | \
	sed 's,@CFLAGS\@,$(MYCFLAGS),' | \
	sed 's,@CC\@,$(CC),' | \
	sed 's,@host\@,@host@,' >$@

CLEANFILES = Makefile.@host@ Makefile.@host@_double *~ \\\#*




