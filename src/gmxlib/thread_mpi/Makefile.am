# Cross-platform threading library/MPI implementation.
#
# This Makefile.am can be used both for linking the final executables against
# the libthread_mpi.a object created here, or for linking into libraries.
# 


AM_CPPFLAGS= -I$(top_srcdir)/include

noinst_LTLIBRARIES = libthread_mpi.la

# Note that not all .c files are compiled directly: some of them 
# are #included (some multiple times) from other source files.
libthread_mpi_la_SOURCES = alltoall.h      lock.c          reduce_fast.c \
			   barrier.c       nothreads.c     scatter.h \
			   bcast.h         once.c          settings.h \
			   collective.c    p2p.c           threads.c \
			   comm.c          p2p_protocol.h  tmpi_init.c \
			   errhandler.c    p2p_send_recv.h tmpi_ops.h \
			   event.c         p2p_wait.h      topology.c \
			   gather.h        profile.c       type.c \
			   group.c         profile.h       winthreads.c \
			   hwinfo.c        pthreads.c      winthreads.h \
			   impl.h          pthreads.h \
			   list.c          reduce.h

CLEANFILES     = *.la *~ \\\#*

