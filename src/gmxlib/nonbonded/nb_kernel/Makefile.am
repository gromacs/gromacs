
AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

noinst_LTLIBRARIES = libnb_kernel.la

#
# We cannot use the standard "built sources" directive,
# since the make_nb_kernels program needs to be compiled with a non-mpi
# compiler in some cases. Supercomputer centers often prevent
# initializing MPI libraries on frontend nodes, and we need
# to run the small program to create all the sources.

FORTRAN_KERNEL_SOURCES = \
	nb_kernel010_f.f nb_kernel020_f.f nb_kernel030_f.f nb_kernel100_f.f \
	nb_kernel101_f.f nb_kernel102_f.f nb_kernel103_f.f nb_kernel104_f.f \
	nb_kernel110_f.f nb_kernel111_f.f nb_kernel112_f.f nb_kernel113_f.f \
	nb_kernel114_f.f nb_kernel120_f.f nb_kernel121_f.f nb_kernel122_f.f \
	nb_kernel123_f.f nb_kernel124_f.f nb_kernel130_f.f nb_kernel131_f.f \
	nb_kernel132_f.f nb_kernel133_f.f nb_kernel134_f.f nb_kernel200_f.f \
	nb_kernel201_f.f nb_kernel202_f.f nb_kernel203_f.f nb_kernel204_f.f \
	nb_kernel210_f.f nb_kernel211_f.f nb_kernel212_f.f nb_kernel213_f.f \
	nb_kernel214_f.f nb_kernel220_f.f nb_kernel221_f.f nb_kernel222_f.f \
	nb_kernel223_f.f nb_kernel224_f.f nb_kernel230_f.f nb_kernel231_f.f \
	nb_kernel232_f.f nb_kernel233_f.f nb_kernel234_f.f nb_kernel300_f.f \
	nb_kernel301_f.f nb_kernel302_f.f nb_kernel303_f.f nb_kernel304_f.f \
	nb_kernel310_f.f nb_kernel311_f.f nb_kernel312_f.f nb_kernel313_f.f \
	nb_kernel314_f.f nb_kernel320_f.f nb_kernel321_f.f nb_kernel322_f.f \
	nb_kernel323_f.f nb_kernel324_f.f nb_kernel330_f.f nb_kernel331_f.f \
	nb_kernel332_f.f nb_kernel333_f.f nb_kernel334_f.f nb_kernel400_f.f \
	nb_kernel410_f.f nb_kernel420_f.f nb_kernel430_f.f 

#

C_KERNEL_SOURCES = \
	nb_kernel010_c.c nb_kernel020_c.c nb_kernel030_c.c nb_kernel100_c.c \
	nb_kernel101_c.c nb_kernel102_c.c nb_kernel103_c.c nb_kernel104_c.c \
	nb_kernel110_c.c nb_kernel111_c.c nb_kernel112_c.c nb_kernel113_c.c \
	nb_kernel114_c.c nb_kernel120_c.c nb_kernel121_c.c nb_kernel122_c.c \
	nb_kernel123_c.c nb_kernel124_c.c nb_kernel130_c.c nb_kernel131_c.c \
	nb_kernel132_c.c nb_kernel133_c.c nb_kernel134_c.c nb_kernel200_c.c \
	nb_kernel201_c.c nb_kernel202_c.c nb_kernel203_c.c nb_kernel204_c.c \
	nb_kernel210_c.c nb_kernel211_c.c nb_kernel212_c.c nb_kernel213_c.c \
	nb_kernel214_c.c nb_kernel220_c.c nb_kernel221_c.c nb_kernel222_c.c \
	nb_kernel223_c.c nb_kernel224_c.c nb_kernel230_c.c nb_kernel231_c.c \
	nb_kernel232_c.c nb_kernel233_c.c nb_kernel234_c.c nb_kernel300_c.c \
	nb_kernel301_c.c nb_kernel302_c.c nb_kernel303_c.c nb_kernel304_c.c \
	nb_kernel310_c.c nb_kernel311_c.c nb_kernel312_c.c nb_kernel313_c.c \
	nb_kernel314_c.c nb_kernel320_c.c nb_kernel321_c.c nb_kernel322_c.c \
	nb_kernel323_c.c nb_kernel324_c.c nb_kernel330_c.c nb_kernel331_c.c \
	nb_kernel332_c.c nb_kernel333_c.c nb_kernel334_c.c nb_kernel400_c.c \
	nb_kernel410_c.c nb_kernel420_c.c nb_kernel430_c.c

if GMX_FORTRAN
  KERNEL_SOURCES = $(FORTRAN_KERNEL_SOURCES) nb_kernel_f77sync.c
else
  KERNEL_SOURCES = $(C_KERNEL_SOURCES)
endif




libnb_kernel_la_SOURCES = \
	nb_kernel010.h   nb_kernel020.h   nb_kernel030.h   nb_kernel100.h   \
	nb_kernel101.h   nb_kernel102.h   nb_kernel103.h   nb_kernel104.h   \
	nb_kernel110.h   nb_kernel111.h   nb_kernel112.h   nb_kernel113.h   \
	nb_kernel114.h   nb_kernel120.h   nb_kernel121.h   nb_kernel122.h   \
	nb_kernel123.h   nb_kernel124.h   nb_kernel130.h   nb_kernel131.h   \
	nb_kernel132.h   nb_kernel133.h   nb_kernel134.h   nb_kernel200.h   \
	nb_kernel201.h   nb_kernel202.h   nb_kernel203.h   nb_kernel204.h   \
	nb_kernel210.h   nb_kernel211.h   nb_kernel212.h   nb_kernel213.h   \
	nb_kernel214.h   nb_kernel220.h   nb_kernel221.h   nb_kernel222.h   \
	nb_kernel223.h   nb_kernel224.h   nb_kernel230.h   nb_kernel231.h   \
	nb_kernel232.h   nb_kernel233.h   nb_kernel234.h   nb_kernel300.h   \
	nb_kernel301.h   nb_kernel302.h   nb_kernel303.h   nb_kernel304.h   \
	nb_kernel310.h   nb_kernel311.h   nb_kernel312.h   nb_kernel313.h   \
	nb_kernel314.h   nb_kernel320.h   nb_kernel321.h   nb_kernel322.h   \
	nb_kernel323.h   nb_kernel324.h   nb_kernel330.h   nb_kernel331.h   \
	nb_kernel332.h   nb_kernel333.h   nb_kernel334.h   nb_kernel400.h   \
	nb_kernel410.h   nb_kernel420.h   nb_kernel430.h   		    \
	nb_kernel.c      nb_kernel.h					    \
	$(KERNEL_SOURCES)



# mknb is special - we cant use a noinst_PROGRAMS target, since it might have 
# to be compiled with a non-MPI compiler whose files can be executed by make:

EXTRA_DIST       =         \
	mknb.c						\
	mknb_metacode.c		mknb_metacode.h		\
	mknb_common.c		mknb_common.h		\
	mknb_declarations.c	mknb_declarations.h	\
	mknb_outerloop.c	mknb_outerloop.h	\
	mknb_innerloop.c	mknb_innerloop.h	\
	mknb_interactions.c	mknb_interactions.h	




BUILD_COMPILE = $(CC_FOR_BUILD) $(DEFS) $(DEFAULT_INCLUDES)



mknb.o:			mknb.c
			$(BUILD_COMPILE) -c -o $@ $<
mknb_common.o:		mknb_common.c
			$(BUILD_COMPILE) -c -o $@ $<
mknb_metacode.o:	mknb_metacode.c
			$(BUILD_COMPILE) -c -o $@ $<
mknb_declarations.o:	mknb_declarations.c
			$(BUILD_COMPILE) -c -o $@ $<
mknb_outerloop.o:	mknb_outerloop.c
			$(BUILD_COMPILE) -c -o $@ $<
mknb_innerloop.o:	mknb_innerloop.c
			$(BUILD_COMPILE) -c -o $@ $<
mknb_interactions.o:	mknb_interactions.c
			$(BUILD_COMPILE) -c -o $@ $<

MKNB_OBJ = \
	mknb.o mknb_metacode.o mknb_common.o mknb_declarations.o \
	mknb_outerloop.o mknb_innerloop.o mknb_interactions.o

mknb:			$(MKNB_OBJ)
			$(BUILD_COMPILE) -o $@ $(MKNB_OBJ)


if GMX_DOUBLE
  DOUBLEFLAG = -double
endif

if GMX_THREAD
  THREADFLAG = -threads
endif

if GMX_SOFTWARE_SQRT
  SQRTFLAG   = -software_invsqrt
endif

if GMX_PREFETCH_FORCES
  PREFETCHFLAG = -prefetch_forces
endif

if GMX_FORTRAN
  FORTRANFLAG = -fortran
endif

MKNB_FLAGS = $(DOUBLEFLAG) $(THREADFLAG) $(SQRTFLAG) \
             $(PREFETCHFLAG) $(FORTRANFLAG)

#
# The kernels should only be built once, instead of once per
# source file when using parallel builds (make -j N). To solve this
# the mknb program writes out a stamp file after the kernels, which
# we use as a dummy dependency.
#
kernel-stamp:			mknb Makefile
				$(RM) kernel-stamp
				./mknb $(MKNB_FLAGS)

$(KERNEL_SOURCES): kernel-stamp


CLEANFILES     = *.la *~ \\\#* mknb nb_kernel*_f.f nb_kernel*_c.c log kernel-stamp 

