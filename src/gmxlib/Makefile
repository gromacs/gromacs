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
OBJ		= force.o 	ns.o 		nsgrid.o     	\
		  lutab.o	calch.o		calcvir.o 	\
		  random.o 	pbc.o  		wnblist.o 	\
		  bondfree.o 	nrnb.o 		maths.o   	\
		  nsb.o		mshift.o	ifunc.o		\
		  fnbf.o	inloopc.o	typedefs.o	\
		  disre.o	init.o 		stat.o 		\
		  mvdata.o	steep.o		sim_util.o	\
		  update.o 	tgroup.o	dummies.o	\
                  block_tx.o	mvxvf.o		runner.o	\
		  ebin.o	rbin.o		vcm.o		\
		  mdebin.o	splittop.o	mdatom.o	\
		  shakef.o	cshake.o	csettle.o	\
		  clincs.o	clincsld.o	calcmu.o	\
		  filenm.o 	txtdump.o 	\
		  confio.o 	rdgroup.o 	\
		  names.o 	statutil.o	vdw.o		\
		  xvgr.o	readinp.o	\
		  wman.o	replace.o	pdbio.o  	\
		  xtcio.o	trxio.o		xdrd.o		\
		  h_db.o	genhydro.o	pargs.o		\
		  fatal.o 	copyrite.o 	memdump.o 	\
		  smalloc.o 	string2.o 	gbutil.o	\
		  macros.o 	ffscanf.o 	futil.o 	\
		  symtab.o 	invblock.o 	fftgrid.o	\
		  strdb.o 	3dview.o 	rando.o		\
		  buffer.o	main.o		network.o 	\
		  nrama.o	do_fit.o	writeps.o	\
		  rmpbc.o	matio.o 	nrjac.o   	\
		  wgms.o 	princ.o		pppm.o		\
		  lrutil.o	ghat.o		tables.o	\
		  edsam.o	congrad.o	javaio.o        \
		  index.o       \
		  gmxfio.o	trnio.o		tpxio.o		\
		  enxio.o
		  
ifeq ($(USEF77),yes)
OBJ		+= inloopf.o	fshake.o	fsettle.o	\
		   flincs.o	flincsld.o	doconv.o
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
		ghat.o		lrutil.o	fftgrid.o
		$(LD) -o $@ $^  $(SYSLIBS) 
gentab:		gentab.o
		$(LD) -o $@ $^  $(SYSLIBS) 
lfft:		lfft.o
		$(LD) -o $@ $^ $(SYSLIBS)

include	Make.dep

