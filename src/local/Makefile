#
#	@(#) Makefile 1.120 10/21/97
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
LIB		= 
XOBJ		= do_gct.o init_sh.o gctio.o glaasje.o ionize.o relax_sh.o
PROGRAMS	= my_dssp	  xmdrun	total		do_shift
#
include		$(STDTARGET)
#

all:		bin

xmdrun:		xmdrun.o $(XOBJ) calcfdev.o
		$(LD) -o $@ $^ $(RUNSYSLIBS)

copyrgt:	copyrgt.o1;
		$(LD) -o $@ $^ $(SYSLIBS) 

my_dssp:	my_dssp.o	dsspcore.o
		$(LD) -o $@ $^ $(SYSLIBS) 

hrefify:	hrefify.o
		$(LD) -o $@ $^ $(SYSLIBS) 

total:		total.o
		$(FLD) -o $@ $^ $(SYSLIBS)

do_shift:	do_shift.o
		$(LD) -o $@ $^ $(SYSLIBS)
testpt:		testpt.o ionize.o
		$(LD) -o $@ $^ $(SYSLIBS) -lm
anaf:		anaf.o
		$(LD) -o $@ $^ $(SYSLIBS) -lm

include	Make.dep
