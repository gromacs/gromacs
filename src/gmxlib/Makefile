#!gmake
#
#	$Id$
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990-1994, Groningen University
#
#	Makefile for gmxlib
#
ifeq ($(GMXHOME),)
all:
	@echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/makef/Makefile.def
#
CPUMAKEF	= ../makef/Makefile.$(GMXCPU)
M4 		= m4
M4DIR		= .
M4FILE		= $(M4DIR)/$(GMXCPU).m4
M4MAC		= $(M4DIR)/gmxmac.m4
LIB		= $(LIBDIR)/libgmx.$(GMXCPU).a
XTCLIB		= $(LIBDIR)/libxtcf.a
MKDEP_PRE	= $(MAKE) fnbf.c
MKDEP_POST	= $(RM) fnbf.c
SONAME		= libgmx.$(GMXCPU).so
TMPFILES        = mkinl   callf77.c libxdrf.c axp_asm.s \
		  inner.h innerf.f innerc.c fnbf.c callf77A.c 
OBJ		= \
	3dview.o	atomprop.o	block_tx.o	bondfree.o	\
	buffer.o	calcgrid.o	calch.o		callf77A.o	\
	confio.o	copyrite.o	disre.o		do_fit.o	\
	enxio.o		ewald_util.o	fatal.o		ffscanf.o	\
	filenm.o	futil.o		gbutil.o	fnbf.o		\
	gmxfio.o	ifunc.o		index.o		cinvsqrtdata.o	\
	crecipdata.o    invblock.o	javaio.o	macros.o	\
	main.o		maths.o		matio.o		memdump.o	\
	mshift.o	mvdata.o	mvxvf.o		names.o		\
	network.o	nrama.o		nrjac.o		nrnb.o		\
	pargs.o		pbc.o		pdbio.o		princ.o		\
	rando.o		random.o	rbin.o		rdgroup.o	\
	readinp.o	replace.o	rmpbc.o		shift_util.o	\
	sortwater.o	smalloc.o	stat.o		statutil.o	\
	strdb.o		string2.o	symtab.o	tpxio.o		\
	trnio.o		trxio.o		txtdump.o	typedefs.o	\
	viewit.o	wgms.o		wman.o		writeps.o	\
	xdrd.o		xtcio.o		xvgr.o         

MKINL_OBJ = \
	mkinl.o	        	mkinl_declarations.o	mkinl_outerloop.o \
	mkinl_innerloop.o   	mkinl_calcdist.o	mkinl_invsqrt.o \
	mkinl_recip.o           mkinl_interactions.o	metacode.o 

ifeq ($(USEF77),yes)
OBJ		+= innerf.o
ifeq ($(USE_THREADS),yes)
OBJ             += threadsync.o # f77 calls a C syncing routine
endif # no threads
else # not f77
OBJ		+= innerc.o
endif

ifeq ($(HAVE_MOTIF),yes)
OBJ		+= mgmx.o widget.o
endif

ifeq ($(USE_XDR),yes)
OBJ 		+= libxdrf.o ftocstr.o
else
OBJ		+= dumxdrf.o
endif

ifneq (,$(findstring USE_SSE_AND_3DNOW,$(SYSDEFS)))
ifeq ($(GMXCPU),lnx)
OBJ 		+= x86cpu.o x86asm.o inner_3dnow.o inner_sse.o
endif
endif

ifneq (,$(findstring USE_AXP_ASM,$(SYSDEFS)))
ifeq ($(GMXCPU),dec)
OBJ 		+= axp_asm.o
endif
ifeq ($(GMXCPU),ded)
OBJ 		+= axp_asm.o
endif
ifeq ($(GMXCPU),alx)
OBJ 		+= axp_asm.o
endif
ifeq ($(GMXCPU),ald)
OBJ 		+= axp_asm.o
endif
endif

# Some CPU dependencies...
ifeq ($(GMXCPU),cm5)
OBJ		+= cm5io.o
else
ifeq ($(USE_PVM3),yes)
OBJ		+= pvmio.o 
else
ifeq ($(USE_MPI),yes)
OBJ		+= mpiio.o
else
OBJ		+= libnet.o
endif
endif
endif


#
include		$(STDTARGET)
#

libxdrf.c:	libxdrf.m4 $(M4FILE) $(CPUMAKEF)
		$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c

callf77.c:	callf77.m4 $(M4FILE) $(CPUMAKEF)
		$(M4) $(M4FILE) callf77.m4 > callf77.c

callf77A.c:	callf77A.m4 $(M4FILE) $(CPUMAKEF)
		$(M4) $(M4FILE) callf77A.m4 > callf77A.c

$(XTCLIB)(libxdrf.o):	libxdrf.c
		$(CC) $(C-FLAGS) -c libxdrf.c
		$(AR) $(ARFLAGS) $(LIB) libxdrf.o

$(XTCLIB):	$(XTCLIB)(ftocstr.o)	$(XTCLIB)(libxdrf.o)	\
		$(XTCLIB)(xtciof.o)

xtclib:		$(XTCLIB)

testfft:	testfft.o
		$(LD) -o $@ $^ $(SYSLIBS) -lcomplib.sgimath
testtab:	testtab.o
		$(LD) -o $@ $^ $(SYSLIBS)

$(LIB)(libxdrf.o):	libxdrf.c
		$(CC) $(C-FLAGS) -c libxdrf.c
		$(AR) $(ARFLAGS) $(LIB) libxdrf.o

$(LIB)(callf77.o):	callf77.c $(M4FILE)
		$(CC) $(C-FLAGS) -c callf77.c
		$(AR) $(ARFLAGS) $(LIB) callf77.o

repfirst:	repfirst.o
		$(LD) -o $@ $^ $(SYSLIBS)

testlr:		testlr.o	ewald.o
		$(LD) -p -o $@ $^  $(SYSLIBS) $(FFTLIBS) $(SYSLIBS)

optghat:	optghat.o	ewald.o	\
		nrutil.o	powell.o	linmin.o	\
		brent.o		f1dim.o		mnbrak.o
		$(LD) -o $@ $^  $(SYSLIBS) $(FFTLIBS)

lfft:		lfft.o
		$(LD) -o $@ $^ $(SYSLIBS)
fitter:		fitter.o
		$(LD) -o $@ $^ $(SYSLIBS) $(FFTLIBS)
mgmxtest:	mgmxtest.o
		$(LD) -o $@ $^  $(SYSLIBS)
calcrf:		calcrf.o
		$(LD) -o $@ $^  $(SYSLIBS)


$(LIB)(callf77A.o):	callf77A.c 
		$(CC) $(C-FLAGS) -c callf77A.c
		$(AR) $(ARFLAGS) $(LIB) callf77A.o


# Innerloop related stuff
$(LIB)(shakef.o):	$(CPUMAKEF)

fnbf.c:			fnbf.m4 $(M4MAC) $(M4FILE) $(CPUMAKEF)
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) fnbf.m4 > fnbf.c

$(LIB)(fnbf.o):		fnbf.c inner.h $(M4MAC) $(M4FILE) $(CPUMAKEF)
			$(CC) $(C-FLAGS) -c fnbf.c
			$(AR) $(ARFLAGS) $(LIB) fnbf.o

x86asm.o:		x86asm.s
			@echo "If your build bails out here, it is most probably because"
			@echo "you are using a version of nasm which does not support the"
			@echo "*EXTENDED* 3dnow instruction set. Vanilla nasm 0.98 has to"
			@echo "be patched to enable this. You'll find both these patches and"
			@echo "a precompiled version under software on the gromacs homepage."
ifeq ($(GMXCPU),win)
			nasm -f win32 -o x86asm.o x86asm.s
else
			nasm -f elf -o x86asm.o x86asm.s
endif

inner_3dnow.o:		inner_3dnow.s
			@echo "If your build bails out here, it is most probably because"
			@echo "you are using a version of nasm which does not support the"
			@echo "*EXTENDED* 3dnow instruction set. Vanilla nasm 0.98 has to"
			@echo "be patched to enable this. You'll find both these patches and"
			@echo "a precompiled version under software on the gromacs homepage."
ifeq ($(GMXCPU),win)
			nasm -f win32 -o inner_3dnow.o inner_3dnow.s
else
			nasm -f elf -o inner_3dnow.o inner_3dnow.s
endif
inner_sse.o:		inner_sse.s
ifeq ($(GMXCPU),win)
			nasm -f win32 -o inner_sse.o inner_sse.s
else
			nasm -f elf -o inner_sse.o inner_sse.s
endif

axp_asm.o:		axp_asm.S
			$(CC) -E axp_asm.S > axp_asm.s
			$(CC) -c axp_asm.s -o axp_asm.o

ifeq ($(CC_MKINL),)
CC_MKINL	= $(CC) $(C-FLAGS)
endif
ifeq ($(LD_MKINL),)
LD_MKINL	= $(LD)
endif

ctabdata.o:		ctabdata.c $(CPUMAKEF)

mkinl.o:		mkinl.c $(CPUMAKEF)
			$(CC_MKINL) -c mkinl.c

mkinl_declarations.o:	mkinl_declarations.c  $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_declarations.c

mkinl_outerloop.o:	mkinl_outerloop.c  $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_outerloop.c

mkinl_innerloop.o:	mkinl_innerloop.c  $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_innerloop.c

mkinl_calcdist.o:	mkinl_calcdist.c  $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_calcdist.c

mkinl_invsqrt.o:	mkinl_invsqrt.c $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_invsqrt.c

mkinl_recip.o:		mkinl_recip.c $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_recip.c

mkinl_interactions.o:	mkinl_interactions.c  $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_interactions.c

metacode.o:		metacode.c $(CPUMAKEF)
			$(CC_MKINL) -c metacode.c

mkinl:			$(MKINL_OBJ) $(CPUMAKEF)
			$(LD_MKINL) -o $@ $(MKINL_OBJ)

inner.h:		inner.m4 $(CPUMAKEF)
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h

$(LIB)(innerf.o):	mkinl inner.h $(CPUMAKEF)
			./mkinl fortran
			$(F77) $(FFLAGS) -c innerf.f
			$(AR) $(ARFLAGS) $(LIB) innerf.o

$(LIB)(innerc.o):	mkinl inner.h Makefile $(CPUMAKEF)
			./mkinl c
			$(CC) $(SYSDEFS) $(C-FLAGS) -c innerc.c
			$(AR) $(ARFLAGS) $(LIB) innerc.o

innerf.f:		mkinl Makefile $(CPUMAKEF)
			./mkinl fortran

innerc.c:		mkinl inner.h Makefile $(CPUMAKEF)
			./mkinl c

tst:			tst.o
			$(FLD) -o tst tst.o $(SYSLIBS)

$(GMXHOME)/src/makef/Makefile.def:
	@echo "Wrong GROMACS environment, Please source GMXRC in this source tree"

include	Make.dep
endif
