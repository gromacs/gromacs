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
XOBJ		= do_gct.o init_sh.o runner.o gctio.o update.o glaasje.o
PROGRAMS	= xmdrun	my_dssp	copyrgt	total	do_shift hrefify
#
include		$(STDTARGET)
#

all:		bin

xmdrun:		xmdrun.o $(XOBJ)
		$(LD) -o $@ $^ $(SYSLIBS) 

copyrgt:	copyrgt.o
		$(LD) -o $@ $^ $(SYSLIBS) 

my_dssp:	my_dssp.o	dsspcore.o
		$(LD) -o $@ $^ $(SYSLIBS) 

hrefify:	hrefify.o
		$(LD) -o $@ $^ $(SYSLIBS) 

total:		total.o
		$(LD) -o $@ $^ $(SYSLIBS)

do_shift:	do_shift.o
		$(LD) -o $@ $^ $(SYSLIBS)

runner.c:	$(GMXHOME)/src/kernel/SCCS/s.runner.c
		$(GET) $^

include	Make.dep
