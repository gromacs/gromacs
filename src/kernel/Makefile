#!gmake
#
#	@(#) Makefile 1.112 4/15/97
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990-1997, Groningen University
#
# 	Makefile for Gromacs kernel
#
ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/makef/Makefile.def
#
# Common object files
PGOBJ		= toputil.o 	topdirs.o 	topexcl.o
MDLIBS		= $(LIBDIR)/libgmx.$(GMXCPU).a	\
		  $(LIBDIR)/libmd.$(GMXCPU).a	\
		  $(LIBDIR)/libfftw.$(GMXCPU).a
# Objects for grompp
GOBJ		= \
		topio.o 	toppush.o 	\
		topcat.o 	topshake.o 	convparm.o 	\
		tomorse.o	sorting.o 	splitter.o	\
		readir.o	dum_parm.o	calcgrid.o	\
		$(PGOBJ)
# Objects for pdb2gmx
POBJ		= \
		pdb2top.o 	gen_ad.o 	\
		pgutil.o	resall.o	ter_db.o	\
		hizzie.o	xlate.o		specbond.o	\
		add_par.o	gen_dum.o	$(PGOBJ)	
ifneq ($(PARALLEL_ONLY_MDRUN),yes)
PROGRAMS	= \
		pdb2gmx 	grompp		tpbconv		\
		mdrun		nmrun		gmxcheck	\
		gmxdump		luck
else
PROGRAMS	= mdrun
endif

include		$(STDTARGET)

all:		bin

grompp:		grompp.o $(GOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

tpbconv:	$(PGOBJ)	tpbconv.o	readir.o
		$(LD) -o $@ $^ $(SYSLIBS)

pdb2gmx:	pdb2gmx.o $(POBJ) 
		$(LD) -o $@ $^ $(SYSLIBS) 

mdrun:		mdrun.o $(MDLIBS)
		$(LD) $(RUNLDOPT) -o $@ mdrun.o $(RUNSYSLIBS)

mk_ghat:	mk_ghat.o
		$(LD) -o $@ $^ $(SYSLIBS) $(FFTLIBS) $(SYSLIBS)

nmrun:		nmrun.o
		$(LD) $(RUNLDOPT) -o $@ nmrun.o $(RUNSYSLIBS)

luck:		luck.o
		$(LD) $^ -o $@ $(SYSLIBS)  

gmxdump:	gmxdump.o
		$(LD) -o $@ $^ $(SYSLIBS)

gmxcheck:	gmxcheck.o	tpbcmp.o
		$(LD) -o $@ $^ $(SYSLIBS)
tpbcmp.o:	tpbcmp.h
gmxcheck.o:	tpbcmp.h

x2top:		x2top.o	$(POBJ) $(PGOBJ) toppush.o nm2type.o
		$(LD) -o $@ $^ $(SYSLIBS) -lm

#
#	Testing stuff
#
test:		test.o
		$(LD) -o $@ $^ $(SYSLIBS)

resall:		resall.o
		$(LD) -o $@ $^ $(SYSLIBS)

ter_db:		ter_db.o h_db.o resall.o pgutil.o $(GOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

opls2rtp:	opls2rtp.o	$(POBJ) toputil.o topdirs.o topexcl.o
		$(LD) -o $@ $^ $(SYSLIBS)

testtop:	testtop.o Makefile
		$(LD) testtop.o $(SYSLIBS) -o $@ 

wb:		wb.o Makefile
		$(LD) wb.o -o $@ $(SYSLIBS) 

ftest:		ftest.o 
		$(LD) ftest.o $(SYSLIBS) -o $@ 

rbin:		rbin.o
		$(LD) $^  $(SYSLIBS) -o $@ 

ebin:		ebin.o
		$(LD) $^  $(SYSLIBS) -o $@ 

align:		align.o
		$(LD) $^  $(SYSLIBS) -o $@ 

ifeq ($(GMXCPU),amb)
tgroup.o:	tgroup.c
		$(CCCCC) $^
steep.o:	steep.c
		$(CCCCC) $^
mvxvf.o:	mvxvf.c
		$(CCCCC) $^
endif

ifeq ($(GMXCPU),am8)
tgroup.o:	tgroup.c
		$(CCCCC) $^
steep.o:	steep.c
		$(CCCCC) $^
mvxvf.o:	mvxvf.c
		$(CCCCC) $^
endif

hallo:		hallo.o
		$(LD) $^ -o $@ $(SYSLIBS)

testf:		testf.o
		$(LD) -o $@ $^ -lm

libdump.a:	dumpargs.o
		$(AR) rv libdump.a $^		
testcomp:	testcomp.o libdump.a
		$(LD) -L. -o $@ testcomp.o -ldump $(SYSLIBS)
#topcat.o:	topcat.c
#		$(CC) CFLAGS=$(LOWOPT) -c topcat.c

include	Make.dep
endif
