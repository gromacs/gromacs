## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

templatedir = ${datadir}/template

# we don't build this, it's meant for the user to play with.
# Since a Makefile is used for the installation, we call the
# template makefile Template.mak, and rename it to Makefile when installing.
# We add the DOUBLE define to CPPFLAGS, so users can do without config.h
LDFLAGS = -L@libdir@ @LDFLAGS@
LIBS = -lmd@LIBSUFFIX@ -lgmx@LIBSUFFIX@ @LIBS@ 
CC = @CC@
if DOUBLE
CFLAGS = @CFLAGS@ -DDOUBLE @INCLUDES@ -I@includedir@ 
else
CFLAGS = @CFLAGS@ @INCLUDES@ -I@includedir@ 
endif

if DOUBLE
template_DATA = template.c README Makefile.@host@_double
else
template_DATA = template.c README Makefile.@host@
endif
EXTRA_DIST = template.c README Template.mak

Makefile.@host@ Makefile.@host@_double: Template.mak Makefile 
	cat $(srcdir)/Template.mak | \
	sed 's,@LDFLAGS\@,$(LDFLAGS),' | \
	sed 's,@LIBS\@,$(LIBS),' | \
	sed 's,@CFLAGS\@,$(CFLAGS),' | \
	sed 's,@CC\@,$(CC),' | \
	sed 's,@host\@,@host@,' >$@

install-data-hook:
if DOUBLE
	(cd $(DESTDIR)$(templatedir) && test ! -e Makefile && $(LN_S) Makefile.@host@_double Makefile ; exit 0)
else
	(cd $(DESTDIR)$(templatedir) && test ! -e Makefile && $(LN_S) Makefile.@host@ Makefile ; exit 0)
endif
CLEANFILES = Makefile.@host@ Makefile.@host@_double *~ \\\#*




