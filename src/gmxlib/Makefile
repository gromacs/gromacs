#
#	$Id$
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990-1994, Groningen University
#
#	Makefile for gmxlib
#
include		$(GMXHOME)/src/makef/Makefile.def
#
M4 		= m4
M4FILE		= $(GMXCPU).m4
LIB		= $(LIBDIR)/libgmx.$(GMXCPU).a
SONAME		= $(LIBDIR)/libgmx.$(GMXCPU).so
OBJ		= \
		3dview.o	block_tx.o	bondfree.o	\
		buffer.o	calch.o		calcmu.o	\
		calcvir.o	callf77.o	clincs.o	\
		clincsld.o	confio.o	copyrite.o	\
		coupling.o	csettle.o	cshake.o	\
		disre.o		do_fit.o	dummies.o	\
		ebin.o		edsam.o		enxio.o		\
		fatal.o		ffscanf.o	fftgrid.o	\
		filenm.o	fnbf.o		force.o		\
		futil.o		gbutil.o	genhydro.o	\
		ghat.o		gmxfio.o	h_db.o		\
		ifunc.o		index.o		init.o		\
		inloopc.o	invblock.o	javaio.o	\
		lrutil.o	lutab.o		macros.o	\
		main.o		mass.o		maths.o		\
		matio.o		mdatom.o	mdebin.o	\
		memdump.o	mshift.o	mvdata.o	\
		mvxvf.o		names.o		network.o	\
		nrama.o		nrjac.o		nrnb.o		\
		ns.o		nsb.o		nsgrid.o	\
		pargs.o		pbc.o		pdbio.o		\
		poisson.o	pppm.o		princ.o		\
		psgather.o	pssolve.o	psspread.o	\
		rando.o		random.o	rbin.o		\
		rdgroup.o	readinp.o	replace.o	\
		rmpbc.o		shakef.o	sim_util.o	\
		smalloc.o	splittop.o	stat.o		\
		statutil.o	strdb.o		string2.o	\
		symtab.o	tables.o	tgroup.o	\
		tpxio.o		trnio.o		trxio.o		\
		txtdump.o	typedefs.o	update.o	\
		vcm.o		vdw.o		wgms.o		\
		wman.o		wnblist.o	writeps.o	\
		xdrd.o		xtcio.o		xvgr.o

ifeq ($(USEF77),yes)
OBJ		+= inloopf.o	fshake.o	fsettle.o	\
		   flincs.o	flincsld.o	fql77.o
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
#		$(RM) inloopf.f inloopf.o

libxtcf.a:	xtciof.o	libxdrf.o	ftocstr.o	#xdrf.o
		$(AR) $@ $^
		$(RANLIB) $@
		$(RM) $(LIBDIR)/$@
		-chmod $(LIBMODE) $@
		mv $@ $(LIBDIR)

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

include	Make.dep
