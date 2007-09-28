## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la

EXTRA_DIST              = README  

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

EXTRA_PROGRAMS 	=  	sigeps 		copyrgt 	mkyaw 	\
			prfn 		options 	hrefify \
			mkice 		total 		hexamer \
			optwat 		addquote 	do_shift \
			g_anavel 	my_dssp 	anaf 	\
			ehole 		pmetest 	gen_table

# Note: you don't have to list sources for "prog" if it is the single file prog.c

my_dssp_SOURCES 	= my_dssp.c dsspcore.c

total_SOURCES		= total.f

ehole_SOURCES		= ehole.c ehdata.c ehanal.c

pmetest_SOURCES   	= pmetest.c

CLEANFILES   = 	*~ \\\#*










