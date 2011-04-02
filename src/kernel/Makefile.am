## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES  = libgmxpreprocess@LIBSUFFIX@.la 

pkgconfigdir = ${libdir}/pkgconfig
pkgconfig_DATA = libgmxpreprocess@LIBSUFFIX@.pc

EXTRA_DIST = libgmxpreprocess.pc.cmakein

LDADD = $(lib_LTLIBRARIES) ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la 

libgmxpreprocess@LIBSUFFIX@_la_LIBADD       = ../mdlib/libmd@LIBSUFFIX@.la
libgmxpreprocess@LIBSUFFIX@_la_DEPENDENCIES =

libgmxpreprocess@LIBSUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ $(PTHREAD_LIBS)

libgmxpreprocess@LIBSUFFIX@_la_SOURCES = \
add_par.c	add_par.h	\
compute_io.c	compute_io.h	\
convparm.c 	convparm.h	\
fflibutil.c	fflibutil.h	\
gen_ad.c 	\
gen_vsite.c	gen_vsite.h	\
genhydro.c	genhydro.h	\
gpp_atomtype.c		\
gpp_bond_atomtype.c	gpp_bond_atomtype.h	\
h_db.c		h_db.h		\
hackblock.c		\
hizzie.c	hizzie.h	\
pdb2top.c	\
pgutil.c	pgutil.h	\
readir.c	readir.h	\
readpull.c	\
resall.c	\
sorting.c 	sorting.h	\
specbond.c	specbond.h	\
ter_db.c	ter_db.h	\
tomorse.c	gpp_tomorse.h	\
topdirs.c	topdirs.h	\
topexcl.c	topexcl.h	\
topio.c 	topio.h		\
toppush.c 	toppush.h	\
topshake.c 	topshake.h	\
toputil.c	\
tpbcmp.c	tpbcmp.h	\
vsite_parm.c	vsite_parm.h	\
xlate.c		xlate.h	

bin_PROGRAMS = \
	grompp		mdrun		tpbconv	  	pdb2gmx		\
	g_protonate	g_luck		gmxdump				\
	gmxcheck	g_x2top

g_x2top_SOURCES = g_x2top.c nm2type.c g_x2top.h

mdrun_SOURCES = \
	gctio.c 	\
	ionize.c 	ionize.h 	xmdrun.h	\
	do_gct.c 	repl_ex.c	repl_ex.h	\
	xutils.c	runner.c	md.c		mdrun.c		\
	genalg.c	genalg.h	md_openmm.h	md_openmm.c

if GMX_FAHCORE
  noinst_LTLIBRARIES = libfahcore.la
  libfahcore_la_SOURCES = $(mdrun_SOURCES) 
  libfahcore_la_CFLAGS = $(AM_CFLAGS)
  libfahcore_la_LDFLAGS = -all-static $(AM_LDFLAGS)
  libfahcore_la_LIBADD = libgmxpreprocess@LIBSUFFIX@.la
endif

fahcore: libfahcore.la 

# shortcut to install only mdrun

install-mdrun: mdrun$(EXEEXT)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	if test -f mdrun$(EXEEXT); then \
	  f=`echo mdrun$(EXEEXT)|sed 's/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/'`; \
	  echo " $(INSTALL_PROGRAM_ENV) $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) mdrun$(EXEEXT) $(DESTDIR)$(bindir)/$$f"; \
	  $(INSTALL_PROGRAM_ENV) $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) mdrun$(EXEEXT) $(DESTDIR)$(bindir)/$$f; \
	else :; fi 

CLEANFILES = *~ \\\#*

