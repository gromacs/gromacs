# Cross-platform threading library/MPI implementation.
#
# This Makefile.am can be used both for linking the final executables against
# the libthread_mpi.a object created here, or for linking into libraries.
# 


AM_CPPFLAGS= -I$(top_srcdir)/include

noinst_LTLIBRARIES = libthread_mpi.la

# again, we assume that we're using pthreads if we're using autotools.
libthread_mpi_la_SOURCES = alltoall.c      impl.h          pthreads.h \
			   barrier.c       list.c          reduce.c \
			   bcast.c         lock.c          reduce_fast.c \
			   collective.c    once.c          scatter.c \
			   collective.h    p2p.h           settings.h \
			   comm.c          p2p_protocol.c  tmpi_init.c \
			   errhandler.c    p2p_send_recv.c tmpi_ops.h \
			   event.c         p2p_wait.c      topology.c \
			   gather.c        profile.c       type.c \
			   group.c         profile.h       winthreads.c \
			   pthreads.c      winthreads.h	   numa_malloc.c


CLEANFILES     = *.la *~ \\\#*

