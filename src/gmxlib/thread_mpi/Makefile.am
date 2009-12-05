# Cross-platform threading library/MPI implementation.
#
# This Makefile.am can be used both for linking the final executables against
# the libthread_mpi.a object created here, or for linking into libraries.
# 


AM_CPPFLAGS= -I$(top_srcdir)/include

noinst_LTLIBRARIES = libthread_mpi.la

# Note that not all .c files are compiled directly: some of them 
# are #included (some multiple times) from other source files.
libthread_mpi_la_SOURCES = threads.c tmpi_init.c errhandler.c type.c group.c\
			   comm.c topology.c send_recv.c collective.c list.c

CLEANFILES     = *.la *~ \\\#*

