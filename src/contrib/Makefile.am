## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\" -I/sw/include

LDADD = ../kernel/libgmxpreprocess@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la ../tools/libgmxana@LIBSUFFIX@.la


EXTRA_DIST              = README  

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

EXTRA_PROGRAMS 	=  	copyrgt 	\
			hrefify 	\
			compnl	\
			gen_table	ehole   g_sdf\
			do_multiprot

compnl_SOURCES = compnl.c
compnl_LDADD   = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la

g_sdf_SOURCES  = gmx_sdf.c g_sdf.c

ehole_SOURCES = ehdata.c ehdata.h ehanal.c ehole.c

# Note: you don't have to list sources for "prog" if it is the single file prog.c

CLEANFILES   = 	*~ \\\#*
