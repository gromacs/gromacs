# Convenience library for statistics routine - not installed

AM_CPPFLAGS= -I$(top_srcdir)/include 

noinst_LTLIBRARIES = libstatistics.la

libstatistics_la_SOURCES =	\
	gmx_statistics.c	histogram.c

LDADD =	../libgmx@LIBSUFFIX@.la

EXTRA_PROGRAMS = gmx_statistics_test

CLEANFILES     = *.la *~ \\\#* 
