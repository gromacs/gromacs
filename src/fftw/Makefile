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

FCR_BLOCKS= fcr_1.o fcr_2.o fcr_3.o fcr_4.o fcr_5.o fcr_6.o       \
            fcr_7.o fcr_8.o fcr_9.o fcr_10.o fcr_11.o fcr_12.o    \
	    fcr_13.o fcr_14.o fcr_15.o fcr_16.o fcr_32.o fcr_64.o \
	    fcr_128.o      	

FHB_BLOCKS= fhb_2.o fhb_3.o fhb_4.o fhb_5.o fhb_6.o fhb_7.o       \
	    fhb_8.o fhb_9.o fhb_10.o fhb_16.o fhb_32.o 

FRC_BLOCKS= frc_1.o frc_2.o frc_3.o frc_4.o frc_5.o frc_6.o       \
	    frc_7.o frc_8.o frc_9.o frc_10.o frc_11.o frc_12.o    \
	    frc_13.o frc_14.o frc_15.o frc_16.o frc_32.o frc_64.o \
	    frc_128.o

FHF_BLOCKS= fhf_2.o fhf_3.o fhf_4.o fhf_5.o fhf_6.o fhf_7.o       \
	    fhf_8.o fhf_9.o fhf_10.o fhf_16.o fhf_32.o

FN_BLOCKS=  fn_1.o fn_2.o fn_3.o fn_4.o fn_5.o fn_6.o fn_7.o      \
	    fn_8.o fn_9.o fn_10.o fn_11.o fn_12.o fn_13.o fn_14.o \
            fn_15.o fn_16.o fn_32.o fn_64.o      

TWID_BLOCKS= ftw_2.o ftw_3.o ftw_4.o ftw_5.o ftw_6.o ftw_7.o      \
	     ftw_8.o ftw_9.o ftw_10.o ftw_16.o ftw_32.o ftw_64.o

NOTWI_BLOCKS= fni_1.o fni_2.o fni_3.o fni_4.o fni_5.o fni_6.o     \
	      fni_7.o fni_8.o fni_9.o fni_10.o fni_11.o fni_12.o  \
	      fni_13.o fni_14.o fni_15.o fni_16.o fni_32.o fni_64.o

TWIDI_BLOCKS= ftwi_2.o ftwi_3.o ftwi_4.o ftwi_5.o ftwi_6.o        \
	      ftwi_7.o ftwi_8.o ftwi_9.o ftwi_10.o ftwi_16.o      \
              ftwi_32.o ftwi_64.o

OTHEROBJ = 	rconfig.o       rexec.o         rexec2.o        rfftwnd.o \
		rgeneric.o	rplanner.o	fftwnd.o	malloc.o  \
		putils.o	twiddle.o	wisdom.o	planner.o \
		config.o	generic.o	executor.o	timer.o	  \
		rader.o		wisdomio.o	

MPIOBJ	=	rfftwnd_mpi.o	transpose_mpi.o	TOMS_transpose.o sched.o  \
		fftwnd_mpi.o

THREADOBJ =     executor_threads.o  rexec_threads.o   rfftwnd_threads.o   \
		fftw_threads.o	    fftwnd_threads.o	rexec2_threads.o  

MPI_THREADOBJ = rfftwnd_mpi_threads.o
# threaded mpi is not a part of fftw, but our hack 
# this should be useful e.g. on SMP nodes on IBM SP2

BLOCKS = $(FCR_BLOCKS) $(FHF_BLOCKS) $(FHB_BLOCKS) $(FRC_BLOCKS) \
	 $(FN_BLOCKS) $(NOTWI_BLOCKS) $(TWID_BLOCKS) $(TWIDI_BLOCKS)

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

#ifeq ($(USE_MPI),yes)
#test_transpose_mpi:	test_transpose_mpi.o     test_main.o	$(OBJ)
#			$(LD) -o $@ $^ $(SYSLIBS)
#
#rfftwnd_mpi_test:	rfftwnd_mpi_test.o     test_main.o	$(OBJ)			
#			$(LD) -o $@ $^ $(SYSLIBS)
#endif
#
#ifeq ($(USE_THREADS),yes)
#rfftw_threads_test:	rfftw_threads_test.o     test_main.o	$(OBJ)
#			$(LD) -o $@ $^ $(SYSLIBS)
#
#fftw_threads_test:	fftw_threads_test.o	     test_main.o $(OBJ)			
#			$(LD) -o $@ $^ $(SYSLIBS)
#endif

include Make.dep
endif
