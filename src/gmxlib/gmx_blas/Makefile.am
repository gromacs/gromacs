# Convenience library for optional built-in BLAS routines - not installed.

AM_CPPFLAGS= -I$(top_srcdir)/include 


noinst_LTLIBRARIES = libblas.la

libblas_la_SOURCES =	\
	dasum.c  	dcopy.c  	dgemm.c  	dger.c		\
	drot.c		dswap.c  	dsyr2.c   	dtrmm.c		\
	dtrsm.c		daxpy.c  	ddot.c   	dgemv.c  	\
	dnrm2.c  	dscal.c  	dsymv.c  	dsyr2k.c  	\
	dtrmv.c		idamax.c					\
	sasum.c         scopy.c         sgemm.c         sger.c          \
	srot.c          sswap.c         ssyr2.c         strmm.c         \
	strsm.c         saxpy.c         sdot.c          sgemv.c         \
	snrm2.c         sscal.c         ssymv.c         ssyr2k.c        \
	strmv.c         isamax.c                                        

CLEANFILES     = *.la *~ \\\#* 
