#
#	@(#) Makefile 1.112 4/15/97
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990-1997, Groningen University
#
# 	Makefile for Gromacs kernel
#
include		$(GMXHOME)/src/makef/Makefile.def
#
# Empty library
LIB		= $(LIBDIR)/libmd.$(CPU).a
# Objects for mdrun and nmrun
RUNOBJ		= runner.o 	md.o		steep.o 	\
		  congrad.o	nm.o	
OBJ		= $(RUNOBJ)
# Common object files
PGOBJ		= toputil.o 	topdirs.o 	topexcl.o	
# Objects for grompp
GOBJ		= grompp.o	topio.o 	toppush.o 	\
		  topcat.o 	topshake.o 	convparm.o 	\
		  tomorse.o	sorting.o 	splitter.o	\
		  readir.o	$(PGOBJ)
# Objects for pdb2gmx
POBJ		= pdb2gmx.o	pdb2top.o 	gen_ad.o 	\
		  pgutil.o	resall.o	ter_db.o	\
		  hizzie.o	xlate.o		specbond.o	\
		  gen_dum.o	add_par.o	$(PGOBJ)	

#
ifneq ($(CPU),amb)
ifneq ($(CPU),860)
PROGRAMS	= pdb2gmx 	grompp		tpbconv		\
		  nmrun		gmxcheck	luck	mk_ghat
endif
endif
PROGRAMS	+= mdrun	gmxdump
#
include		$(STDTARGET)
#

all:		bin

grompp:		$(GOBJ)
		$(LD) -o $@ $^ $(SYSLIBS)

tpbconv:	$(PGOBJ)	tpbconv.o	readir.o
		$(LD) -o $@ $^ $(SYSLIBS)

pdb2gmx:	$(POBJ) 
		$(LD) -o $@ $^ $(SYSLIBS) 

mdrun:		mdrun.o $(RUNOBJ)
		$(LD) $(RUNLDOPT) -o $@ $^ $(LINKLIB) $(FFTLIBS) $(SYSLIBS)

mk_ghat:	mk_ghat.o
		$(LD) -o $@ $^ $(SYSLIBS) $(FFTLIBS) $(SYSLIBS)

nmrun:		nmrun.o $(RUNOBJ)
		$(LD) $(RUNLDOPT) -o $@ $^ $(SYSLIBS) $(FFTLIBS) $(SYSLIBS)

luck:		luck.o
		$(LD) $^ -o $@ $(SYSLIBS)  

gmxdump:	gmxdump.o
		$(LD) -o $@ $^ $(SYSLIBS)

gmxcheck:	gmxcheck.o	tpbcmp.o
		$(LD) -o $@ $^ $(SYSLIBS)
tpbcmp.o:	tpbcmp.h
gmxcheck.o:	tpbcmp.h

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

ifeq ($(CPU),amb)
tgroup.o:	tgroup.c
		$(CCCCC) $^
steep.o:	steep.c
		$(CCCCC) $^
mvxvf.o:	mvxvf.c
		$(CCCCC) $^
endif

ifeq ($(CPU),am8)
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

include	Make.dep
