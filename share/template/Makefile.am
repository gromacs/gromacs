## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

templatedir = ${datadir}/template

# we don't build this, it's meant for the user to play with.
# Since a Makefile is used for the installation, we call the
# template makefile Template.mak, and rename it to Makefile when installing.
LDFLAGS = @LDFLAGS@
LIBS = @LIBS@
CPPFLAGS = @CPPFLAGS@
CFLAGS = @CFLAGS@
CC = @CC@

template_DATA = template.c README Makefile.@host@
EXTRA_DIST = template.c README Template.mak

#BUILT_SOURCES = Makefile.@host@

Makefile.@host@: Template.mak Makefile 
	cat $(srcdir)/Template.mak | sed 's,@libdir\@,$(libdir),' | \
	sed 's,@LDFLAGS\@,$(LDFLAGS),' | \
	sed 's,@LIBS\@,$(LIBS),' | \
	sed 's,@CPPFLAGS\@,$(CPPFLAGS),' | \
	sed 's,@CFLAGS\@,$(CFLAGS),' | \
	sed 's,@CC\@,$(CC),' | \
	sed 's,@includedir\@,$(includedir),' | \
	sed 's,@host\@,@host@,' >$@

install-data-hook:
	(cd $(DESTDIR)$(templatedir) && test ! -e Makefile && $(LN_S) Makefile.@host@ Makefile ; exit 0)

CLEANFILES = Makefile.@host@ *~ \\\#*




