## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la 


bin_PROGRAMS = \
	grompp		mdrun		tpbconv	  	pdb2gmx		\
	protonate	luck		gmxdump				\
	gmxcheck  	x2top		ffscan


grompp_SOURCES = \
	topio.c 	toppush.c 	topcat.c 	topshake.c 	\
	convparm.c 	tomorse.c	sorting.c 	splitter.c	\
	vsite_parm.c	readir.c  	add_par.c	topexcl.c	\
	toputil.c	topdirs.c	grompp.c	add_par.h	\
	splitter.h	topio.h		toputil.h	convparm.h	\
	sorting.h	topdirs.h	toppush.h	vsite_parm.h	\
	readir.h	topcat.h	topexcl.h	topshake.h	\
	compute_io.h	compute_io.c

tpbconv_SOURCES = \
	readir.c	toputil.c	topdirs.c	add_par.c	\
	topexcl.c	tpbconv.c	add_par.h	toputil.h	\
	topdirs.h	readir.h	topexcl.h

pdb2gmx_SOURCES = \
	hizzie.c	xlate.c		specbond.c	ter_db.c	\
	h_db.c		genhydro.c     	pdb2top.c 	gen_ad.c 	\
	gen_vsite.c	pgutil.c	resall.c	hackblock.c	\
	add_par.c	topexcl.c 	toputil.c	topdirs.c	\
	pdb2gmx.c	gen_ad.h	add_par.h	h_db.h		\
	pdb2top.h	resall.h	toputil.h	xlate.h		\
	gen_vsite.h	hackblock.h	pgutil.h	ter_db.h	\
	topdirs.h	genhydro.h	hizzie.h	specbond.h	\
	topexcl.h

#rtp2xml_SOURCES = \
#	hizzie.c	xlate.c		specbond.c	ter_db.c	\
#	h_db.c		genhydro.c     	pdb2top.c 	gen_ad.c 	\
#	gen_dum.c	pgutil.c	resall.c	hackblock.c	\
#	add_par.c	topexcl.c 	toputil.c	topdirs.c	\
#	rtp2xml.c	gen_ad.h	add_par.h	h_db.h		\
#	pdb2top.h	resall.h	toputil.h	xlate.h		\
#	gen_dum.h	hackblock.h	pgutil.h	ter_db.h	\
#	topdirs.h	genhydro.h	hizzie.h	specbond.h	\
#	topexcl.h

protonate_SOURCES = \
	hackblock.c	ter_db.c	h_db.c		genhydro.c	\
	pgutil.c	resall.c	add_par.c	topexcl.c	\
	toputil.c	topdirs.c	protonate.c	add_par.h	\
	h_db.h		resall.h	toputil.h	hackblock.h	\
	pgutil.h	ter_db.h	topdirs.h	topexcl.h

gmxcheck_SOURCES = \
	gmxcheck.c	tpbcmp.c	tpbcmp.h

x2top_SOURCES = \
	toppush.c 	nm2type.c	pdb2top.c 	gen_ad.c 	\
	gen_vsite.c	pgutil.c	resall.c	hackblock.c	\
	add_par.c	topexcl.c	toputil.c	topdirs.c	\
	x2top.c		add_par.h	gen_ad.h	pdb2top.h	\
	resall.h	toputil.h	gen_vsite.h	hackblock.h	\
	pgutil.h	topdirs.h	toppush.h	topexcl.h	\
	x2top.h

mdrun_SOURCES = \
	glaasje.c 	glaasje.h 	gctio.c 	init_sh.c 	\
	ionize.c 	ionize.h 	xmdrun.h	\
	do_gct.c 	relax_sh.c	repl_ex.c	repl_ex.h	\
	xutils.c	compute_io.h	compute_io.c \
	md.c		mdrun.c		genalg.c	genalg.h

ffscan_SOURCES = \
	ffscan.c	\
	glaasje.c 	glaasje.h 	gctio.c 	init_sh.c 	\
	ionize.c 	ionize.h 	md.c		genalg.c	\
	do_gct.c 	relax_sh.c	repl_ex.c	repl_ex.h	\
	xutils.c	compute_io.h	compute_io.c

# shortcut to install only mdrun

install-mdrun: mdrun$(EXEEXT)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	if test -f mdrun$(EXEEXT); then \
	  f=`echo mdrun$(EXEEXT)|sed 's/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/'`; \
	  echo " $(INSTALL_PROGRAM_ENV) $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) mdrun$(EXEEXT) $(DESTDIR)$(bindir)/$$f"; \
	  $(INSTALL_PROGRAM_ENV) $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) mdrun$(EXEEXT) $(DESTDIR)$(bindir)/$$f; \
	else :; fi 

CLEANFILES = *~ \\\#*


