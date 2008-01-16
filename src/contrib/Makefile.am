## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

LDADD = ../kernel/libgmxpreprocess@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la ../tools/libgmxana@LIBSUFFIX@.la


EXTRA_DIST              = README  

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

EXTRA_PROGRAMS 	=  	copyrgt 	mkyaw 		\
			prfn 		options 	hrefify 	\
			mkice 		total 		hexamer 	\
			optwat 		addquote 	do_shift 	\
			g_anavel  	anaf 		tune_dip	\
			ehole 		pmetest 	gen_table	\
			tune_pol

# Note: you don't have to list sources for "prog" if it is the single file prog.c

total_SOURCES		= total.f

ehole_SOURCES		= ehole.c ehdata.c ehanal.c

pmetest_SOURCES   	= pmetest.c

tune_pol_SOURCES	= tune_pol.c tune_pol.h pol_xml.c poldata.c \
			merge_xml.c

CLEANFILES   = 	*~ \\\#*










