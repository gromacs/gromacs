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
SONAME		= libgmx.$(GMXCPU).so
OBJ		= \
		3dview.o	atomprop.o	block_tx.o	\
		bondfree.o	buffer.o	calch.o		\
		confio.o	copyrite.o	\
		disre.o		do_fit.o	enxio.o		\
		fatal.o		ffscanf.o	filenm.o	\
		futil.o		gbutil.o	genhydro.o	\
		fnbf.o		gmxfio.o	h_db.o		\
		ifunc.o		index.o		ghat.o		\
		invblock.o	javaio.o	lincsutl.o	\
		lutab.o		macros.o	main.o		\
		maths.o		matio.o		lrutil.o	\
		memdump.o	mshift.o	mvdata.o	\
		mvxvf.o		names.o		network.o	\
		nrama.o		nrjac.o		nrnb.o		\
		pargs.o		pbc.o		pdbio.o		\
		princ.o		\
		rando.o		random.o	rbin.o		\
		rdgroup.o	readinp.o	replace.o	\
		rmpbc.o		smalloc.o	stat.o		\
		statutil.o	strdb.o		string2.o	\
		symtab.o	\
		tpxio.o		trnio.o		trxio.o		\
		txtdump.o	typedefs.o	wgms.o		\
		wman.o		writeps.o	\
		xdrd.o		xtcio.o		xvgr.o

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

EXEC_GENTAB	= gentab

# Some CPU dependencies...
ifeq ($(GMXCPU),860)
OBJ		+= idtio.o inout.o delay.o  \
		   comlib.o semlib.o synclib.o  led.o 
EXEC_GENTAB	= (setenv SPC_POLL 1000000 ; s860 -w 1 gentab)

else

ifeq ($(GMXCPU),cm5)
OBJ		+= libsys.o cm5io.o
else
ifeq ($(GMXCPU),amb)
OBJ		+= ambio.o 	ambaled.o	delay.o
EXEC_GENTAB	= ambasrv -c 0 gentab
else
ifeq ($(GMXCPU),am8)
OBJ		+= ambio.o 	delay.o 	ambaled.o
EXEC_GENTAB	= ambasrv -c 0 gentab
else
ifeq ($(GMXCPU),a23)
OBJ		+= ambio.o 	delay.o 	ambaled.o
EXEC_GENTAB	= ambasrv -c 0 gentab
else
ifeq ($(USE_PVM3),yes)
OBJ		+= libsys.o 	pvmio.o 
else
ifeq ($(USE_MPI),yes)
OBJ		+= libsys.o 	mpiio.o
else
OBJ		+= libsys.o	libnet.o
endif
endif
endif
endif
endif
endif
endif

#
include		$(STDTARGET)
#

testfft:	testfft.o
		$(LD) -o $@ $^ $(SYSLIBS) -lcomplib.sgimath
testtab:	testtab.o
		$(LD) -o $@ $^ $(SYSLIBS)

$(LIB)(libxdrf.o):	libxdrf.m4 $(M4FILE)
		$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c
		$(CC) $(C-FLAGS) -c libxdrf.c
		$(AR) $(ARFLAGS) $(LIB) libxdrf.o
		$(RM) libxdrf.c libxdrf.o

prfn:		prfn.o
		$(LD) -o $@ $^ $(SYSLIBS)

repfirst:	repfirst.o
		$(LD) -o $@ $^ $(SYSLIBS)

testlr:		testlr.o	ewald.o
		$(LD) -p -o $@ $^  $(SYSLIBS) $(FFTLIBS) $(SYSLIBS)

optghat:	optghat.o	ewald.o	\
		nrutil.o	powell.o	linmin.o	\
		brent.o		f1dim.o		mnbrak.o
		$(LD) -o $@ $^  $(SYSLIBS) $(FFTLIBS)

gentab:		gentab.o
		$(LD) -o $@ $^  $(SYSLIBS) 

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

# Innerloop related stuff
$(LIB)(lutab.o):	$(CPUMAKEF)
$(LIB)(shakef.o):	$(CPUMAKEF)
fnbf.c:		fnbf.m4 $(M4MAC) $(M4FILE) $(CPUMAKEF)
		m4 $(SYSDEFS) $(M4FILE) $(M4MAC) fnbf.m4 > $@
$(LIB)(fnbf.o):		fnbf.c inner.h $(CPUMAKEF)
mkinl.o:	$(CPUMAKEF)
mkinl:		mkinl.o 
		$(LD) -o $@ mkinl.o
inner.h:	inner.m4 $(M4MAC) $(CPUMAKEF)
		$(M4) $(SYSDEFS) $(M4FILE) $(M4MAC) inner.m4 > inner.h
innerc.c:	mkinl inner.h Makefile $(CPUMAKEF)
		./mkinl 1 1
innerf.f:	mkinl inner.h Makefile $(CPUMAKEF)
		./mkinl 0 1

include	Make.dep
endif


