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
	echo "No GROMACS environment, Please source your GMXRC file"
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
OBJ		= \
	3dview.o	atomprop.o	block_tx.o	bondfree.o	\
	buffer.o	calch.o		callf77A.o	confio.o	\
	copyrite.o	disre.o		do_fit.o	enxio.o		\
	ewald_util.o	fatal.o		ffscanf.o	filenm.o	\
	futil.o		gbutil.o	fnbf.o		gmxfio.o	\
	ifunc.o		index.o		cinvsqrtdata.o	\
	invblock.o	javaio.o	macros.o	\
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
	mkinl_interactions.o	metacode.o 

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
OBJ 		+= x86cpu.o x86cpu_asm.o
endif
endif

ifneq (,$(findstring USE_AXP_ASM,$(SYSDEFS)))
ifeq ($(GMXCPU),dec)
OBJ 		+= axp_asm.o
endif
ifeq ($(GMXCPU),ded)
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

$(XTCLIB)(libxdrf.o):	libxdrf.m4 $(M4FILE)
		$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c
		$(CC) $(C-FLAGS) -c libxdrf.c
		$(AR) $(ARFLAGS) $(LIB) libxdrf.o
		$(RM) libxdrf.c libxdrf.o
$(XTCLIB):	$(XTCLIB)(ftocstr.o)	$(XTCLIB)(libxdrf.o)	\
		$(XTCLIB)(xtciof.o)

xtclib:		$(XTCLIB)

testfft:	testfft.o
		$(LD) -o $@ $^ $(SYSLIBS) -lcomplib.sgimath
testtab:	testtab.o
		$(LD) -o $@ $^ $(SYSLIBS)

$(LIB)(libxdrf.o):	libxdrf.m4 $(M4FILE)
		$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c
		$(CC) $(C-FLAGS) -c libxdrf.c
		$(AR) $(ARFLAGS) $(LIB) libxdrf.o
		$(RM) libxdrf.c libxdrf.o

$(LIB)(callf77.o):	callf77.m4 $(M4FILE)
		$(M4) $(M4FILE) callf77.m4  > callf77.c
		$(CC) $(C-FLAGS) -c callf77.c
		$(AR) $(ARFLAGS) $(LIB) callf77.o
		$(RM) callf77.c callf77.o

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
compnl:		compnl.o
		$(LD) -o $@ $^  $(SYSLIBS)
calcrf:		calcrf.o
		$(LD) -o $@ $^  $(SYSLIBS)

$(LIB)(callf77A.o):	callf77A.m4 $(M4FILE)
		$(M4) $(M4FILE) callf77A.m4  > callf77A.c
		$(CC) $(C-FLAGS) -c callf77A.c
		$(AR) $(ARFLAGS) $(LIB) callf77A.o
		$(RM) callf77A.c callf77A.o

# Innerloop related stuff
$(LIB)(shakef.o):	$(CPUMAKEF)
$(LIB)(fnbf.o):		fnbf.m4 inner.m4 $(M4MAC) $(M4FILE) $(CPUMAKEF)
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) fnbf.m4 > fnbf.c
			$(CC) $(C-FLAGS) -c fnbf.c
			$(AR) $(ARFLAGS) $(LIB) fnbf.o
			$(RM) fnbf.c inner.h fnbf.o

fnbf.c:			fnbf.m4 inner.m4 $(M4MAC) $(M4FILE) $(CPUMAKEF) 
			 $(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) fnbf.m4 > fnbf.c 

x86cpu_asm.o:		x86cpu_asm.s
ifeq ($(GMXCPU),win)
			nasm -f win32 -o x86cpu_asm.o x86cpu_asm.s
else
			nasm -f elf -o x86cpu_asm.o x86cpu_asm.s
endif

axp_asm.o:		axp_asm.S
			$(CC) -E axp_asm.S > axp_asm.s
			$(CC) -c axp_asm.s -o axp_asm.o
			$(RM) axp_asm.s

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

mkinl_interactions.o:	mkinl_interactions.c  $(CPUMAKEF)
			$(CC_MKINL) -c mkinl_interactions.c


metacode.o:		metacode.c $(CPUMAKEF)
			$(CC_MKINL) -c metacode.c

mkinl:			$(MKINL_OBJ) $(CPUMAKEF)
			$(LD_MKINL) -o $@ $(MKINL_OBJ)

$(LIB)(innerf.o):	mkinl inner.m4 $(CPUMAKEF)
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
			./mkinl fortran
			$(F77) $(FFLAGS) -c innerf.f
			$(AR) $(ARFLAGS) $(LIB) innerf.o
			$(RM) innerf.f innerf.o mkinl

$(LIB)(innerc.o):	mkinl inner.m4 $(CPUMAKEF)
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
			./mkinl c
			$(RM) mkinl
			$(CC) $(SYSDEFS) $(C-FLAGS) -c innerc.c
			$(AR) $(ARFLAGS) $(LIB) innerc.o
			$(RM) innerc.c inner.h innerc.o  mkinl

innerf.f:		mkinl Makefile $(CPUMAKEF)
			./mkinl fortran
			$(RM) mkinl.o mkinl
innerc.c:		mkinl inner.m4 Makefile $(CPUMAKEF)
			./mkinl c
			$(RM) mkinl mkinl.o

tst:			tst.o
			$(FLD) -o tst tst.o $(SYSLIBS)

$(GMXHOME)/src/makef/Makefile.def:
	echo "Wrong GROMACS environment, Please source GMXRC in this source tree"

include	Make.dep
endif
