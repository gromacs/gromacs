# Generated automatically from Makefile.in by configure.
# Makefile.in generated automatically by automake 1.4d from Makefile.am

# Copyright 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001
# Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /usr/local/gromacs
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin/i686-pc-linux-gnu
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib/i686-pc-linux-gnu
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

pkgdatadir = $(datadir)/gromacs
pkglibdir = $(libdir)/gromacs
pkgincludedir = $(includedir)/gromacs

top_builddir = ..

ACLOCAL = ${SHELL} /usr/local/gromacs/cvs/gmx/config/missing --run aclocal
AUTOCONF = ${SHELL} /usr/local/gromacs/cvs/gmx/config/missing --run autoconf
AUTOMAKE = ${SHELL} /usr/local/gromacs/cvs/gmx/config/missing --run automake
AUTOHEADER = ${SHELL} /usr/local/gromacs/cvs/gmx/config/missing --run autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
INSTALL_STRIP_FLAG =
INSTALL_HEADER = $(INSTALL_DATA)
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :


host_alias = i686-pc-linux-gnu
host_triplet = i686-pc-linux-gnu
AMDEP = 
AMTAR = ${SHELL} /usr/local/gromacs/cvs/gmx/config/missing --run tar
AWK = gawk
AXP_ASM_OBJ = 
BUILD_CC = cc
CC = cc
CPP = cc -E
DEPDIR = .deps
F77 = 
FLIBS = 
IDENT = /usr/bin/ident
INCLUDES =    -I/usr/X11R6/include 
INNER_C_OBJ = ${inner_c_obj}
INNER_F77_OBJ = 
MAKEINFO = ${SHELL} /usr/local/gromacs/cvs/gmx/config/missing --run makeinfo
MOTIF_OBJ = ${motif_obj}
MPICC = 
NASM = /usr/bin/nasm
NASMFLAGS = -f elf
PACKAGE = gromacs
PAR_OBJ = ${libnet_obj}
RANLIB = ranlib
SSE_OBJ = ${sse_obj}
SUFFIX = 
TDN_OBJ = ${tdn_obj}
VERSION = 3.0b
X86_ASM_OBJ = ${x86_asm_obj}
XDR_OBJ = ${xdr_obj}
_am_include = include
bindir = ${exec_prefix}/bin/i686-pc-linux-gnu
install_sh = ./config/install-sh
libdir = ${exec_prefix}/lib/i686-pc-linux-gnu

#
# Don't edit - this file is generated automatically from Makefile.am
#


#
# Don't edit - this file is generated automatically from Makefile.am
#

# The asm suffix is for intel syntax assembly, and
# the s suffix for at & t syntax.
# S suffix files will be preprocessed by cpp, nasm
# on the other hand can do this directly.


SUFFIXES = .asm .S .F

#################
# We need a second compile command producing executables 
# that can be run on the local host to make the innerloops.
# In most cases this is identical to the normal one, but not
# for cross-compilation. We also need it on parallel machines 
# where the MPI executables cannot be run outside a batch queue.

BUILD_COMPILE = $(BUILD_CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)

SUBDIRS = include gmxlib mdlib kernel tools ngmx
# don't do anything for the local directory
EXTRA_DIST = Makefile.inc

CLEANFILES = *~ \\\#*
subdir = src
mkinstalldirs = $(SHELL) $(top_srcdir)/./config/mkinstalldirs
CONFIG_HEADER = ../src/include/config.h
CONFIG_CLEAN_FILES = 
CFLAGS = -O6 -fomit-frame-pointer -finline-functions -funroll-all-loops -Wall -Wno-unused -malign-double
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
DIST_SOURCES = 
DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)

GZIP_ENV = --best
DIST_SUBDIRS =  $(SUBDIRS)
all: all-recursive
.SUFFIXES:
.SUFFIXES: .asm .S .F
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) $(srcdir)/Makefile.inc
	cd $(top_srcdir) && $(AUTOMAKE) --gnu src/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= CONFIG_LINKS= $(SHELL) ./config.status


# This directory's subdirectories are mostly independent; you can cd
# into them and run `make' without going through this Makefile.
# To change the values of `make' variables: instead of editing Makefiles,
# (1) if the variable is set in `config.status', edit `config.status'
#     (which will cause the Makefiles to be regenerated when you run `make');
# (2) otherwise, pass the desired values on the `make' command line.
all-recursive install-data-recursive install-exec-recursive \
installdirs-recursive install-recursive uninstall-recursive  \
check-recursive installcheck-recursive info-recursive dvi-recursive:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	dot_seen=no; \
	target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    dot_seen=yes; \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	   || case "$$amf" in *=*) exit 1;; *k*) fail=yes;; *) exit 1;; esac; \
	done; \
	if test "$$dot_seen" = "no"; then \
	  $(MAKE) $(AM_MAKEFLAGS) "$$target-am" || exit 1; \
	fi; test -z "$$fail"

mostlyclean-recursive clean-recursive distclean-recursive \
maintainer-clean-recursive:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	dot_seen=no; \
	case "$@" in \
	  distclean-* | maintainer-clean-*) list='$(DIST_SUBDIRS)' ;; \
	  *) list='$(SUBDIRS)' ;; \
	esac; \
	rev=''; for subdir in $$list; do \
	  if test "$$subdir" = "."; then :; else \
	    rev="$$subdir $$rev"; \
	  fi; \
	done; \
	rev="$$rev ."; \
	target=`echo $@ | sed s/-recursive//`; \
	for subdir in $$rev; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	   || case "$$amf" in *=*) exit 1;; *k*) fail=yes;; *) exit 1;; esac; \
	done && test -z "$$fail"
tags-recursive:
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  test "$$subdir" = . || (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) tags); \
	done

tags: TAGS

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	mkid -fID $$unique $(LISP)

TAGS: tags-recursive $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
   if test "$$subdir" = .; then :; else \
	    test -f $$subdir/TAGS && tags="$$tags -i $$here/$$subdir/TAGS"; \
   fi; \
	done; \
	list='$(SOURCES) $(HEADERS) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(ETAGS_ARGS)$$unique$(LISP)$$tags" \
	  || etags $(ETAGS_ARGS) $$tags  $$unique $(LISP)

GTAGS:
	here=`CDPATH=: && cd $(top_builddir) && pwd` \
	  && cd $(top_srcdir) \
	  && gtags -i $$here

distclean-tags:
	-rm -f TAGS ID

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)


distdir: $(DISTFILES)
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pR $$d/$$file $(distdir) \
	    || exit 1; \
	  else \
	    test -f $(distdir)/$$file \
	    || cp -p $$d/$$file $(distdir)/$$file \
	    || exit 1; \
	  fi; \
	done
	for subdir in $(SUBDIRS); do \
	  if test "$$subdir" = .; then :; else \
	    test -d $(distdir)/$$subdir \
	    || mkdir $(distdir)/$$subdir \
	    || exit 1; \
	    (cd $$subdir && \
	      $(MAKE) $(AM_MAKEFLAGS) top_distdir=../$(top_distdir) distdir=../$(distdir)/$$subdir distdir) \
	      || exit 1; \
	  fi; \
	done
info-am:
info: info-recursive
dvi-am:
dvi: dvi-recursive
check-am: all-am
check: check-recursive
installcheck-am:
installcheck: installcheck-recursive
install-exec-am:
install-exec: install-exec-recursive

install-data-am:
install-data: install-data-recursive

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-recursive
uninstall-am:
uninstall: uninstall-recursive
all-am: Makefile
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_STRIP_FLAG=-s install
installdirs: installdirs-recursive
installdirs-am:


mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
	-rm -f Makefile.in
clean: clean-recursive

clean-am: clean-generic mostlyclean-am

distclean: distclean-recursive

distclean-am: clean-am distclean-generic distclean-tags

maintainer-clean: maintainer-clean-recursive

maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-recursive

mostlyclean-am: mostlyclean-generic

.PHONY: all-recursive check check-am check-recursive clean clean-generic \
	clean-recursive distclean distclean-generic distclean-recursive \
	distclean-tags distdir dvi dvi-am dvi-recursive info info-am \
	info-recursive install install-am install-data install-data-am \
	install-data-recursive install-exec install-exec-am \
	install-exec-recursive install-recursive install-strip \
	installcheck installcheck-am installcheck-recursive installdirs \
	installdirs-am installdirs-recursive maintainer-clean \
	maintainer-clean-generic maintainer-clean-recursive mostlyclean \
	mostlyclean-generic mostlyclean-recursive tags tags-recursive \
	uninstall uninstall-am uninstall-recursive

.asm.o:
	$(NASM) $(NASMFLAGS) $< -o $@
.S.s:
	$(CPP) $< > $@
.F.f:	
	$(CPP) $< > $@

mdrun:		
	(cd ${top_builddir}/src/gmxlib && $(MAKE))
	(cd ${top_builddir}/src/mdlib && $(MAKE))
	(cd ${top_builddir}/src/kernel && $(MAKE) mdrun)

install-mdrun:
	cd ${top_builddir}/src/kernel 
	make install-mdrun

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
