# Convenience library for parallel processing - not installed

AM_CPPFLAGS= -I$(top_srcdir)/include 

noinst_LTLIBRARIES = libparallel.la

if THREAD_PARALLEL
if THREAD_PTHREADS
libparallel_la_SOURCES =	gmx_thread_pthreads.c 	thread_mpi.c	
else
libparallel_la_SOURCES =	gmx_thread_no.c 	thread_mpi.c	
endif
else
libparallel_la_SOURCES =	gmx_thread_no.c 	
endif 


LDADD =	../libgmx@LIBSUFFIX@.la

#EXTRA_PROGRAMS = gmx_statistics_test

CLEANFILES     = *.la *~ \\\#* 

