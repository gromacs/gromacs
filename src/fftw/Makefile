#!gmake
#
#	@(#) Makefile 1.154 07/24/97
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Makefile for libfftw
#
ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/makef/Makefile.def
#
LIB	= $(LIBDIR)/libfftw.$(GMXCPU).a
SONAME	= libfftw.$(GMXCPU).so
# Or, on primitive systems with 8.3 filenames:
NOTW_PREFIX=fn_
TWID_PREFIX=ftw_
NOTWI_PREFIX=fni_
TWIDI_PREFIX=ftwi_

NOTW_BLOCKS=$(NOTW_PREFIX)1.o $(NOTW_PREFIX)2.o $(NOTW_PREFIX)3.o     \
	$(NOTW_PREFIX)4.o $(NOTW_PREFIX)5.o $(NOTW_PREFIX)6.o         \
	$(NOTW_PREFIX)7.o $(NOTW_PREFIX)8.o $(NOTW_PREFIX)9.o         \
	$(NOTW_PREFIX)10.o $(NOTW_PREFIX)11.o $(NOTW_PREFIX)12.o      \
	$(NOTW_PREFIX)13.o $(NOTW_PREFIX)14.o $(NOTW_PREFIX)15.o      \
	$(NOTW_PREFIX)16.o $(NOTW_PREFIX)32.o $(NOTW_PREFIX)64.o      

TWID_BLOCKS=$(TWID_PREFIX)2.o $(TWID_PREFIX)3.o $(TWID_PREFIX)4.o     \
	$(TWID_PREFIX)5.o $(TWID_PREFIX)6.o $(TWID_PREFIX)7.o         \
	$(TWID_PREFIX)8.o $(TWID_PREFIX)9.o $(TWID_PREFIX)10.o        \
	$(TWID_PREFIX)16.o $(TWID_PREFIX)32.o $(TWID_PREFIX)64.o

NOTWI_BLOCKS=$(NOTWI_PREFIX)1.o $(NOTWI_PREFIX)2.o $(NOTWI_PREFIX)3.o \
	$(NOTWI_PREFIX)4.o $(NOTWI_PREFIX)5.o $(NOTWI_PREFIX)6.o      \
	$(NOTWI_PREFIX)7.o $(NOTWI_PREFIX)8.o $(NOTWI_PREFIX)9.o      \
	$(NOTWI_PREFIX)10.o $(NOTWI_PREFIX)11.o $(NOTWI_PREFIX)12.o   \
	$(NOTWI_PREFIX)13.o $(NOTWI_PREFIX)14.o $(NOTWI_PREFIX)15.o   \
	$(NOTWI_PREFIX)16.o $(NOTWI_PREFIX)32.o $(NOTWI_PREFIX)64.o

TWIDI_BLOCKS=$(TWIDI_PREFIX)2.o $(TWIDI_PREFIX)3.o                   \
	$(TWIDI_PREFIX)4.o $(TWIDI_PREFIX)5.o $(TWIDI_PREFIX)6.o    \
	$(TWIDI_PREFIX)7.o $(TWIDI_PREFIX)8.o $(TWIDI_PREFIX)9.o    \
	$(TWIDI_PREFIX)10.o $(TWIDI_PREFIX)16.o $(TWIDI_PREFIX)32.o \
        $(TWIDI_PREFIX)64.o

OTHEROBJ = 	timer.o 	config.o 	planner.o 	twiddle.o 	\
		executor.o 	generic.o 	fftwnd.o 	malloc.o 	\
		wisdom.o 	rader.o		putils.o	wisdomio.o

MPIOBJ	=	fftwnd_mpi.o	transpose_mpi.o	TOMS_transpose.o
# Apparently outdated files
# common_io.o	naive.o 

BLOCKS = $(NOTW_BLOCKS) $(NOTWI_BLOCKS) $(TWID_BLOCKS) $(TWIDI_BLOCKS) 

OBJ = $(BLOCKS) $(OTHEROBJ)

#ifeq ($(USE_MPI),yes)
#OBJ		+= $(MPIOBJ)
#endif

#
include		$(STDTARGET)
#

#ifeq ($(USE_MPI),yes)
#test_transpose_mpi:	test_transpose_mpi.o	$(OBJ)
#			$(LD) -o $@ $^ $(SYSLIBS)
#
#test_fftwnd_mpi:	test_fftwnd_mpi.o	$(OBJ)			
#			$(LD) -o $@ $^ $(SYSLIBS)
#
#time_fftwnd_mpi:	time_fftwnd_mpi.o	$(OBJ)
#			$(LD) -o $@ $^ $(SYSLIBS)
#endif

include Make.dep
endif
