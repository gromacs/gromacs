# Makefile.in generated automatically by automake 1.4h from Makefile.am.

# Copyright 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001
# Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.



#
# Don't edit - this file is generated automatically from Makefile.am
#

SHELL = /bin/sh

srcdir = .
top_srcdir = ..

prefix = /usr/local/gromacs
exec_prefix = ${prefix}/mips-sgi-irix6.5/r12000

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include
pkgdatadir = $(datadir)/gromacs
pkglibdir = $(libdir)/gromacs
pkgincludedir = $(includedir)/gromacs
top_builddir = ..

ACLOCAL = ${SHELL} /tmp/gmx/config/missing --run aclocal
AUTOCONF = ${SHELL} /tmp/gmx/config/missing --run autoconf
AUTOMAKE = ${SHELL} /tmp/gmx/config/missing --run automake
AUTOHEADER = ${SHELL} /tmp/gmx/config/missing --run autoheader

INSTALL = /usr/freeware/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL}
INSTALL_HEADER = $(INSTALL_DATA)
transform = s,x,x,
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
host_alias = 
host_triplet = mips-sgi-irix6.5
AMTAR = ${SHELL} /tmp/gmx/config/missing --run tar
AS = cc
ASFLAGS = -n32 -r12000 -mips4 -O3 -OPT:IEEE_arithmetic=3 -OPT:rsqrt=ON -SWP:loop_overhead -INLINE:=ON -LNO:opt=1 -LNO:ou_further=3 -OPT:Olimit=0:roundoff=3:alias=typed -woff 1174 -D__INLINE_INTRINSICS
AWK = mawk
CC = cc
CC_FOR_BUILD = cc
CPP = cc -E
CXX = 
CXXCPP = 
DEPDIR = .deps
DLLTOOL = @DLLTOOL@
ECHO = print -r
EXEEXT = 
F77 = f77
FLIBS =  -L/usr5/home/lindahl/fftw/lib -L/usr/lib32/mips4/r10000 -L/usr/lib32/mips4 -L/usr/lib32 -lftn -lm
GCJ = @GCJ@
GCJFLAGS = @GCJFLAGS@
GMXLIB_LIBOBJS =  innerf.lo f77_wrappers.lo libnet.lo mgmx.lo widget.lo
IDENT = /usr/sbin/ident
INCLUDES =    
INSTALL_STRIP_PROGRAM = ${SHELL} `CDPATH=: && cd config && pwd`/install-sh -c -s
LIBSUFFIX = 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LN_S = ln -s
MDLIB_LIBOBJS =  flincs.lo fsettle.lo fshake.lo
MPICC = 
OBJDUMP = @OBJDUMP@
OBJEXT = o
PACKAGE = gromacs
RANLIB = ranlib
RPM = rpm
SHARED_VERSION_INFO = 1:0:0
STRIP = strip
VERSION = 3.0
X_CFLAGS = 
X_EXTRA_LIBS = 
X_LIBS = 
X_PRE_LIBS =  -lSM -lICE
_am_include = include
_am_quote = 
install_sh = config/install-sh

EXTRA_SCRIPTS = mkcompl  mkhtml  mknroff  mktex 

EXTRA_DIST = gromacs.spec  gromacs-mpi.spec  programs.txt README.patches $(EXTRA_SCRIPTS)
subdir = admin
mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs
CONFIG_HEADER = $(top_builddir)/src/config.h
CONFIG_CLEAN_FILES =
DIST_SOURCES =
DIST_COMMON = Makefile.am Makefile.in
all: all-am

.SUFFIXES:

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

distclean-libtool:
	-rm -f libtool
$(srcdir)/Makefile.in:  Makefile.am  $(top_srcdir)/configure.ac $(ACLOCAL_M4)
	cd $(top_srcdir) && \
	  $(AUTOMAKE) --gnu  admin/Makefile
Makefile:  $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) && \
	  CONFIG_HEADERS= CONFIG_LINKS= \
	  CONFIG_FILES=$(subdir)/$@ $(SHELL) ./config.status
uninstall-info-am:
tags: TAGS
TAGS:


DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)

top_distdir = ..
distdir = $(top_distdir)/$(PACKAGE)-$(VERSION)

distdir: $(DISTFILES)
	@for file in $(DISTFILES); do \
	  if test -f $$file; then d=.; else d=$(srcdir); fi; \
	  dir=`echo "$$file" | sed -e 's,/[^/]*$$,,'`; \
	  if test "$$dir" != "$$file" && test "$$dir" != "."; then \
	    $(mkinstalldirs) "$(distdir)/$$dir"; \
	  fi; \
	  if test -d $$d/$$file; then \
	    cp -pR $$d/$$file $(distdir) \
	    || exit 1; \
	  else \
	    test -f $(distdir)/$$file \
	    || cp -p $$d/$$file $(distdir)/$$file \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
check: check-am
all-am: Makefile

installdirs:

install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES) stamp-h stamp-h[0-9]*

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-generic clean-libtool mostlyclean-am

distclean: distclean-am

distclean-am: clean-am distclean-generic distclean-libtool

dvi: dvi-am

dvi-am:

info: info-am

info-am:

install-data-am:

install-exec-am:

install-info: install-info-am

install-man:

installcheck-am:

maintainer-clean: maintainer-clean-am

maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-generic mostlyclean-libtool

uninstall-am: uninstall-info-am

.PHONY: all all-am check check-am clean clean-generic clean-libtool \
	distclean distclean-generic distclean-libtool distdir dvi \
	dvi-am info info-am install install-am install-data \
	install-data-am install-exec install-exec-am install-info \
	install-info-am install-man install-strip installcheck \
	installcheck-am installdirs maintainer-clean \
	maintainer-clean-generic mostlyclean mostlyclean-generic \
	mostlyclean-libtool uninstall uninstall-am uninstall-info-am

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
