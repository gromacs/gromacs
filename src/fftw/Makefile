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
FCR_PREFIX=fcr_
FHB_PREFIX=fhb_
FRC_PREFIX=frc_
FHF_PREFIX=fhf_
NOTW_PREFIX=fn_
TWID_PREFIX=ftw_
NOTWI_PREFIX=fni_
TWIDI_PREFIX=ftwi_

FCR_BLOCKS=$(FCR_PREFIX)1.o $(FCR_PREFIX)2.o $(FCR_PREFIX)3.o     \
	$(FCR_PREFIX)4.o $(FCR_PREFIX)5.o $(FCR_PREFIX)6.o         \
	$(FCR_PREFIX)7.o $(FCR_PREFIX)8.o $(FCR_PREFIX)9.o         \
	$(FCR_PREFIX)10.o $(FCR_PREFIX)11.o $(FCR_PREFIX)12.o      \
	$(FCR_PREFIX)13.o $(FCR_PREFIX)14.o $(FCR_PREFIX)15.o      \
	$(FCR_PREFIX)16.o $(FCR_PREFIX)32.o $(FCR_PREFIX)64.o      \
	$(FCR_PREFIX)128.o      

FHB_BLOCKS=$(FHB_PREFIX)2.o $(FHB_PREFIX)3.o $(FHB_PREFIX)4.o     \
	$(FHB_PREFIX)5.o $(FHB_PREFIX)6.o $(FHB_PREFIX)7.o         \
	$(FHB_PREFIX)8.o $(FHB_PREFIX)9.o $(FHB_PREFIX)10.o        \
	$(FHB_PREFIX)16.o $(FHB_PREFIX)32.o 

FRC_BLOCKS=$(FRC_PREFIX)1.o $(FRC_PREFIX)2.o $(FRC_PREFIX)3.o \
	$(FRC_PREFIX)4.o $(FRC_PREFIX)5.o $(FRC_PREFIX)6.o      \
	$(FRC_PREFIX)7.o $(FRC_PREFIX)8.o $(FRC_PREFIX)9.o      \
	$(FRC_PREFIX)10.o $(FRC_PREFIX)11.o $(FRC_PREFIX)12.o   \
	$(FRC_PREFIX)13.o $(FRC_PREFIX)14.o $(FRC_PREFIX)15.o   \
	$(FRC_PREFIX)16.o $(FRC_PREFIX)32.o $(FRC_PREFIX)64.o	\
	$(FRC_PREFIX)128.o

FHF_BLOCKS=$(FHF_PREFIX)2.o $(FHF_PREFIX)3.o                   \
	$(FHF_PREFIX)4.o $(FHF_PREFIX)5.o $(FHF_PREFIX)6.o    \
	$(FHF_PREFIX)7.o $(FHF_PREFIX)8.o $(FHF_PREFIX)9.o    \
	$(FHF_PREFIX)10.o $(FHF_PREFIX)16.o $(FHF_PREFIX)32.o

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


OTHEROBJ = 	rconfig.o       rexec.o         rexec2.o        rfftwnd.o \
		rgeneric.o	rplanner.o	fftwnd.o	malloc.o  \
		putils.o	twiddle.o	wisdom.o	planner.o \
		config.o	generic.o	executor.o	timer.o	\
		rader.o		wisdomio.o	

MPIOBJ	=	rfftwnd_mpi.o	transpose_mpi.o	TOMS_transpose.o sched.o \
		fftwnd_mpi.o

THREADOBJ =     executor_threads.o  rexec_threads.o   rfftwnd_threads.o   \
		fftw_threads.o	    fftwnd_threads.o	rexec2_threads.o  

MPI_THREADOBJ = rfftwnd_mpi_threads.o
# threaded mpi is not a part of fftw, but our hack 
# this should be useful e.g. on SMP nodes on IBM SP2

BLOCKS = $(FCR_BLOCKS) $(FHF_BLOCKS) $(FHB_BLOCKS) $(FRC_BLOCKS) \
	$(NOTW_BLOCKS) $(NOTWI_BLOCKS) $(TWID_BLOCKS) $(TWIDI_BLOCKS)

OBJ = $(BLOCKS) $(OTHEROBJ)

ifeq ($(USE_MPI),yes)
OBJ		+= $(MPIOBJ)
ifeq ($(USE_THREADS),yes)
OBJ		+= $(MPI_THREADOBJ)
endif
endif

ifeq ($(USE_THREADS),yes)
OBJ		+= $(THREADOBJ)
endif

#
include		$(STDTARGET)
#

ifeq ($(USE_MPI),yes)
test_transpose_mpi:	test_transpose_mpi.o     test_main.o	$(OBJ)
			$(LD) -o $@ $^ $(SYSLIBS)

rfftwnd_mpi_test:	rfftwnd_mpi_test.o     test_main.o	$(OBJ)			
			$(LD) -o $@ $^ $(SYSLIBS)
endif

ifeq ($(USE_THREADS),yes)
rfftw_threads_test:	rfftw_threads_test.o     test_main.o	$(OBJ)
			$(LD) -o $@ $^ $(SYSLIBS)

fftw_threads_test:	fftw_threads_test.o	     test_main.o $(OBJ)			
			$(LD) -o $@ $^ $(SYSLIBS)
endif

testprog:		testprog.o $(OBJ)
			mpcc_r -o $@ $^ -lm
rfftw_mpi_threads_test: rfftw_mpi_threads_test.o test_main.o   $(OBJ)
			mpcc_r -o $@ $^ -lm
include Make.dep
endif
