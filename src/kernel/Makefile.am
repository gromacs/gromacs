## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

include $(srcdir)/../Makefile.inc

# This is necessary for VPATH to work, 
# but it can probably be done in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include
LDFLAGS = @LDFLAGS@ -L${top_builddir}/src/gmxlib -L${top_builddir}/src/mdlib
LDADD = -lmdXXX_SUFFIX_XXX -lgmxXXX_SUFFIX_XXX


bin_PROGRAMS = \
	gromppXXX_SUFFIX_XXX		mdrunXXX_SUFFIX_XXX		\
	tpbconvXXX_SUFFIX_XXX	  	pdb2gmxXXX_SUFFIX_XXX		\
	protonateXXX_SUFFIX_XXX		nmrunXXX_SUFFIX_XXX	  	\
	luckXXX_SUFFIX_XXX		gmxdumpXXX_SUFFIX_XXX		\
	gmxcheckXXX_SUFFIX_XXX  	x2topXXX_SUFFIX_XXX				



gromppXXX_SUFFIX_XXX_SOURCES = \
	topio.c 	toppush.c 	topcat.c 	topshake.c 	\
	convparm.c 	tomorse.c	sorting.c 	splitter.c	\
	dum_parm.c	readir.c  	add_par.c	topexcl.c	\
	toputil.c	topdirs.c	grompp.c	add_par.h	\
	splitter.h	topio.h		toputil.h	convparm.h	\
	sorting.h	topdirs.h	toppush.h	dum_parm.h	\
	readir.h	topcat.h	topexcl.h	topshake.h

mdrunXXX_SUFFIX_XXX_SOURCES = \
	mdrun.c

tpbconvXXX_SUFFIX_XXX_SOURCES = \
	readir.c	toputil.c	topdirs.c	add_par.c	\
	topexcl.c	tpbconv.c	add_par.h	toputil.h	\
	topdirs.h	readir.h	topexcl.h

pdb2gmxXXX_SUFFIX_XXX_SOURCES = \
	hizzie.c	xlate.c		specbond.c	ter_db.c	\
	h_db.c		genhydro.c     	pdb2top.c 	gen_ad.c 	\
	gen_dum.c	pgutil.c	resall.c	hackblock.c	\
	add_par.c	topexcl.c 	toputil.c	topdirs.c	\
	pdb2gmx.c	gen_ad.h	add_par.h	h_db.h		\
	pdb2top.h	resall.h	toputil.h	xlate.h		\
	gen_dum.h	hackblock.h	pgutil.h	ter_db.h	\
	topdirs.h	genhydro.h	hizzie.h	specbond.h	\
	topexcl.h

protonateXXX_SUFFIX_XXX_SOURCES = \
	hackblock.c	ter_db.c	h_db.c		genhydro.c	\
	pgutil.c	resall.c	add_par.c	topexcl.c	\
	toputil.c	topdirs.c	protonate.c	add_par.h	\
	h_db.h		resall.h	toputil.h	hackblock.h	\
	pgutil.h	ter_db.h	topdirs.h	topexcl.h

nmrunXXX_SUFFIX_XXX_SOURCES = \
	nmrun.c

luckXXX_SUFFIX_XXX_SOURCES = \
	luck.c

gmxdumpXXX_SUFFIX_XXX_SOURCES = \
	gmxdump.c

gmxcheckXXX_SUFFIX_XXX_SOURCES = \
	gmxcheck.c	tpbcmp.c	tpbcmp.h

x2topXXX_SUFFIX_XXX_SOURCES = \
	toppush.c 	nm2type.c	pdb2top.c 	gen_ad.c 	\
	gen_dum.c	pgutil.c	resall.c	hackblock.c	\
	add_par.c	topexcl.c	toputil.c	topdirs.c	\
	x2top.c		add_par.h	gen_ad.h	pdb2top.h	\
	resall.h	toputil.h	gen_dum.h	hackblock.h	\
	pgutil.h	topdirs.h	toppush.h	topexcl.h	\
	x2top.h

install-mdrun: mdrunXXX_SUFFIX_XXX
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	if test -f mdrunXXX_SUFFIX_XXX; then \
		echo "$(INSTALL_PROGRAM) $(INSTALL_STRIP_FLAG) mdrunXXX_SUFFIX_XXX $(DESTDIR)$(bindir)/mdrunXXX_SUFFIX_XXX"; \
	     $(INSTALL_PROGRAM) $(INSTALL_STRIP_FLAG) mdrunXXX_SUFFIX_XXX $(DESTDIR)$(bindir)/mdrunXXX_SUFFIX_XXX; \
	fi

# clean things explicitly, since the target names might have changed

CLEANFILES   = 	${bin_PROGRAMS}	${EXTRA_PROGRAMS}	\
		*_mpi *_d *~ \\\#*


