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

# libraries for mdrun
MDLIBS	= $(LIBDIR)/libgmx.$(GMXCPU).a	\
	  $(LIBDIR)/libmd.$(GMXCPU).a	\
	  $(LIBDIR)/libfftw.$(GMXCPU).a

# some Common Objects:
# Objects for grompp, pdb2gmx, tpbconv and protonate
PGTRXOBJ = toputil.o	topdirs.o

# Objects for grompp, pdb2gmx, tpbconv and x2top
PGTXOBJ	= add_par.o	topexcl.o	$(PGTRXOBJ)

# Objects for pdb2gmx, protonate and x2top
PRXOBJ	= pgutil.o	resall.o	hackblock.o	$(PGTXOBJ)

# Objects for pdb2gmx and protonate
PROBJ	= ter_db.o	h_db.o		genhydro.o	$(PRXOBJ)

# Objects for pdb2gmx and x2top
PXOBJ	= pdb2top.o 	gen_ad.o 	gen_dum.o	$(PRXOBJ)

# Objects for grompp and tpbconv
GTOBJ	= readir.o	$(PGTROBJ)	$(PGTXOBJ)

# program-specific Objects
# Objects for grompp
GOBJ	= \
	topio.o 	toppush.o 	topcat.o 	topshake.o 	\
	convparm.o 	tomorse.o	sorting.o 	splitter.o	\
	dum_parm.o	calcgrid.o	$(GTOBJ)

# Objects for pdb2gmx
POBJ	= \
	hizzie.o	xlate.o		specbond.o	\
	$(PROBJ)	$(PXOBJ)

# Objects for tpbconv
TOBJ	= $(GTOBJ)

# Objects for protonate
ROBJ	= $(PROBJ)	hackblock.o

# Objects for x2top
XOBJ	= $(PXOBJ) toppush.o nm2type.o

ifneq ($(PARALLEL_ONLY_MDRUN),yes)
PROGRAMS = \
	pdb2gmx 	protonate	grompp		tpbconv		\
	mdrun		nmrun		gmxcheck	gmxdump		\
	luck		x2top
else
PROGRAMS = mdrun
endif

include		$(STDTARGET)

$(GMXHOME)/src/makef/Makefile.def:
	echo "Wrong GROMACS environment, Please source GMXRC in this source tree"

all:		bin

grompp:		grompp.o $(GOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

tpbconv:	tpbconv.o $(TOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

pdb2gmx:	pdb2gmx.o $(POBJ) 
		$(LD) -o $@ $^ $(SYSLIBS) 

protonate:	protonate.o $(ROBJ)
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

x2top:		x2top.o	$(XOBJ)
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
