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
SONAME		= libgmx.$(GMXCPU).so
OBJ		= \
	3dview.o	atomprop.o	block_tx.o	bondfree.o	\
	buffer.o	calch.o		callf77A.o	confio.o	\
	copyrite.o	disre.o		do_fit.o	enxio.o		\
	ewald_util.o	fatal.o		ffscanf.o	filenm.o	\
	futil.o		gbutil.o	fnbf.o		gmxfio.o	\
	ifunc.o		index.o		\
	invblock.o	javaio.o	lutab.o		macros.o	\
	main.o		maths.o		matio.o		memdump.o	\
	mshift.o	mvdata.o	mvxvf.o		names.o		\
	network.o	nrama.o		nrjac.o		nrnb.o		\
	pargs.o		pbc.o		pdbio.o		princ.o		\
	rando.o		random.o	rbin.o		rdgroup.o	\
	readinp.o	replace.o	rmpbc.o		shift_util.o	\
	sortwater.o	smalloc.o	stat.o		statutil.o	\
	strdb.o		string2.o	symtab.o	tpxio.o		\
	trnio.o		trxio.o		txtdump.o	typedefs.o	\
	wgms.o		wman.o		writeps.o	xdrd.o		\
	xtcio.o		xvgr.o

ifeq ($(USEF77),yes)
OBJ		+= innerf.o
else
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
$(LIB)(lutab.o):	$(CPUMAKEF)
$(LIB)(shakef.o):	$(CPUMAKEF)
$(LIB)(fnbf.o):		fnbf.m4 inner.m4 $(M4MAC) $(M4FILE) $(CPUMAKEF)
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) fnbf.m4 > fnbf.c
			$(CC) $(C-FLAGS) -c fnbf.c
			$(AR) $(ARFLAGS) $(LIB) fnbf.o
			$(RM) fnbf.c inner.h fnbf.o

fnbf.c:			fnbf.m4 inner.m4 $(M4MAC) $(M4FILE) $(CPUMAKEF) 
			 $(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) fnbf.m4 > fnbf.c 
ifeq ($(CC_MKINL),)
CC_MKINL	= $(CC) $(C-FLAGS)
endif
ifeq ($(LD_MKINL),)
LD_MKINL	= $(LD)
endif

mkinl.o:		mkinl.c $(CPUMAKEF)
			$(CC_MKINL) -g -c $(SYSDEFS) mkinl.c 

mkinl:			mkinl.o
			$(LD_MKINL) -o mkinl mkinl.o

$(LIB)(innerf.o):	mkinl.c $(CPUMAKEF)
			$(CC_MKINL) -g -c $(SYSDEFS) mkinl.c 
			$(LD_MKINL) -o mkinl mkinl.o
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
			./mkinl 0 1
			$(F77) $(FFLAGS) -c innerf.f
			$(AR) $(ARFLAGS) $(LIB) innerf.o
			$(RM) innerf.f innerf.o mkinl.o mkinl

$(LIB)(innerc.o):	mkinl.c inner.m4 $(CPUMAKEF)
			$(CC_MKINL) -g -c $(SYSDEFS) mkinl.c 
			$(LD_MKINL) -o mkinl mkinl.o
			$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
			./mkinl 1 1
			$(RM) mkinl
			$(CC) $(C-FLAGS) -c innerc.c
			$(AR) $(ARFLAGS) $(LIB) innerc.o
			$(RM) innerc.c inner.h innerc.o mkinl.o mkinl

innerf.f:		mkinl Makefile $(CPUMAKEF)
			./mkinl 0 1
			$(RM) mkinl.o mkinl

innerc.c:		mkinl inner.m4 Makefile $(CPUMAKEF)
			./mkinl 1 1
			$(RM) mkinl mkinl.o

tst:			tst.o
			$(FLD) -o tst tst.o $(SYSLIBS)
include	Make.dep
endif
