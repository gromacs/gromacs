#
#	@(#) Makefile 1.154 07/24/97
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
M4FILE		= $(CPU).m4
LIBRARY		= libgmx
OBJ		= 	\
		3dview.o	block_tx.o	bondfree.o	\
		buffer.o	calch.o		calcmu.o	\
		calcvir.o	callf77.o	clincs.o	\
		clincsld.o	confio.o	\
		copyrite.o	coupling.o	\
		csettle.o	cshake.o	disre.o		\
		do_fit.o	dummies.o	ebin.o		\
		edsam.o		enxio.o		fatal.o		\
		ffscanf.o	fftgrid.o	filenm.o	\
		fnbf.o		force.o		futil.o		\
		gbutil.o	genhydro.o	ghat.o		\
		gmxfio.o	h_db.o		ifunc.o		\
		index.o		init.o		inloopc.o	\
		invblock.o	javaio.o	lrutil.o	\
		lutab.o		macros.o	main.o		\
		mass.o						\
		maths.o		matio.o		mdatom.o	\
		mdebin.o	memdump.o	mshift.o	\
		mvdata.o	mvxvf.o		names.o		\
		network.o	nrama.o		nrjac.o		\
		nrnb.o		ns.o		nsb.o		\
		nsgrid.o	pargs.o		pbc.o		\
		pdbio.o		pppm.o		princ.o		\
		rando.o		random.o	rbin.o		\
		rdgroup.o	readinp.o	replace.o	\
		rmpbc.o		shakef.o	\
		sim_util.o	smalloc.o	splittop.o	\
		stat.o		statutil.o	\
		strdb.o		string2.o	symtab.o	\
		tables.o	tgroup.o	tpxio.o		\
		trnio.o		trxio.o		txtdump.o	\
		typedefs.o	update.o	vcm.o		\
		vdw.o		wgms.o		wman.o		\
		wnblist.o	writeps.o	xdrd.o		\
		xtcio.o		xvgr.o

ifeq ($(USEF77),yes)
OBJ		+= inloopf.o	fshake.o	fsettle.o	\
		   flincs.o	flincsld.o	fql77.o
endif

ifeq ($(USE_XDR),yes)
OBJ 		+= libxdrf.o ftocstr.o
else
OBJ		+= dumxdrf.o
endif

EXEC_GENTAB	= gentab

# Some CPU dependencies...
ifeq ($(CPU),860)
OBJ		+= idtio.o inout.o delay.o  \
		   comlib.o semlib.o synclib.o  led.o 
EXEC_GENTAB	= (setenv SPC_POLL 1000000 ; s860 -w 1 gentab)

else

ifeq ($(CPU),cm5)
OBJ		+= libsys.o cm5io.o
else
ifeq ($(CPU),amb)
OBJ		+= ambio.o 	ambaled.o	delay.o
EXEC_GENTAB	= ambasrv -c 0 gentab
else
ifeq ($(CPU),am8)
OBJ		+= ambio.o 	delay.o 	ambaled.o
EXEC_GENTAB	= ambasrv -c 0 gentab
else
ifeq ($(CPU),a23)
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
testtab:	testtab.o	tables.o	lrutil.o
		$(LD) -o $@ $^ $(SYSLIBS)
tstsqrtc:	tstsqrtc.o
		$(LD) -o $@ $^ $(SYSLIBS)
tstsqrtf:	tstsqrtf.o inloopf.o
		$(F77) -o $@ $^ -lm
libxdrf.o:	libxdrf.m4 $(M4FILE)
#		$(RM) libxdrf.c
		$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c
		$(CC) $(C-FLAGS) -c libxdrf.c
		$(RM) libxdrf.c
callf77.o:	callf77.m4 $(M4FILE)
		$(M4) $(M4FILE) callf77.m4  > callf77.c
		$(CC) $(C-FLAGS) -c callf77.c
		$(RM) callf77.c

inloopf.o:	inloopf.m4 
		$(CC) -E $(SYSDEFS) inloopf.m4 | grep -v '#' > inloopf.f
		$(F77) $(FFLAGS) -c inloopf.f
		$(RM) inloopf.f

libxtcf.a:	xtciof.o	libxdrf.o	ftocstr.o
		$(AR) $@ $^
		$(RANLIB) $@
		$(RM) $(LIBDIR)/$@
		-chmod $(LIBMODE) $@
		mv $@ $(LIBDIR)

repfirst:	repfirst.o
		$(LD) -o $@ $^ $(SYSLIBS)

testlr:		testlr.o	ewald.o		pppm.o	\
		ghat.o		lrutil.o	fftgrid.o	poisson.o
		$(LD) -o $@ $^  $(SYSLIBS) $(FFTLIBS)

testwc:		testwc.o	coupling.o
		$(LD) -o $@ $^  $(SYSLIBS)

gentab:		gentab.o
		$(LD) -o $@ $^  $(SYSLIBS) 

lfft:		lfft.o
		$(LD) -o $@ $^ $(SYSLIBS)
fitter:		fitter.o
		$(LD) -o $@ $^ $(SYSLIBS) $(FFTLIBS)

include	Make.dep

