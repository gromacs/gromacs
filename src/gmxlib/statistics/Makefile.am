# Convenience library for statistics routine - not installed

AM_CPPFLAGS= -I$(top_srcdir)/include 

noinst_LTLIBRARIES = libstatistics.la

libgmx_statistics_la_SOURCES =	\
	gmx_statistics.c

LDADD =	../libgmx@LIBSUFFIX@.la

EXTRA_PROGRAMS = gmx_statistics_test

CLEANFILES     = *.la *~ \\\#* 
