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
M4 		= m4
M4FILE		= $(GMXCPU).m4
LIB		= $(LIBDIR)/libmd.$(GMXCPU).a
SONAME		= libmd.$(GMXCPU).so
OBJ		= \
		calcmu.o	calcvir.o	\
		congrad.o	coupling.o	edsam.o		\
		fftgrid.o	fnbf.o		force.o		\
		init.o		lrutil.o	\
		md.o		mdatom.o	nm.o		\
		ns.o		nsb.o		nsgrid.o	\
		poisson.o	pppm.o		\
		psgather.o	pssolve.o	psspread.o	\
		runner.o	sim_util.o	splittop.o	\
		steep.o		tables.o	tgroup.o	\
		update.o	vcm.o		wnblist.o

#
include		$(STDTARGET)
#

testfft:	testfft.o
		$(LD) -o $@ $^ $(SYSLIBS) -lcomplib.sgimath
testtab:	testtab.o
		$(LD) -o $@ $^ $(SYSLIBS)
tstsqrtc:	tstsqrtc.o
		$(LD) -o $@ $^ $(SYSLIBS)
tstsqrtf:	tstsqrtf.o inloopf.o
		$(F77) -o $@ $^ -lm
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
$(LIB)(inloopf.o):	inloopf.m4 
		/lib/cpp $(SYSDEFS) inloopf.m4 | grep -v '#' > inloopf.f
		$(F77) $(FFLAGS) -c inloopf.f
		$(AR) $(ARFLAGS) $(LIB) inloopf.o
		$(RM) inloopf.f inloopf.o

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

include	Make.dep
endif

