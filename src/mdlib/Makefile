#!gmake
#
#	$Id$
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990-1994, Groningen University
#
#	Makefile for mdlib
#
ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/makef/Makefile.def
#
M4 		= m4
M4FILE		= $(GMXHOME)/src/gmxlib/$(GMXCPU).m4
LIB		= $(LIBDIR)/libmd.$(GMXCPU).a
SONAME		= libmd.$(GMXCPU).so
OBJ		= \
		calcmu.o	calcvir.o	callf77.o	\
		congrad.o	constr.o	coupling.o	\
		dummies.o	ebin.o		edsam.o		\
		ewald.o         fftgrid.o	force.o  	\
		ghat.o		init.o		md.o		\
		mdatom.o	mdebin.o	nm.o		\
		ns.o		nsb.o		nsgrid.o	\
		pme.o		poisson.o	pppm.o		\
		psgather.o	pull.o 		pullinit.o	\
		pullio.o	pullutil.o	\
		pssolve.o	psspread.o	runner.o	\
		shakef.o	sim_util.o	splittop.o	\
		steep.o		tables.o	tgroup.o	\
		update.o	vcm.o		wnblist.o

ifeq ($(USEF77),yes)
OBJ		+= fshake.o	fsettle.o	flincs.o
else
OBJ		+= cshake.o	csettle.o	clincs.o
endif

#
include		$(STDTARGET)
#

$(GMXHOME)/src/makef/Makefile.def:
	echo "Wrong GROMACS environment, Please source GMXRC in this source tree"

# This is to force rebuild when the Makefile.CPU has changed
$(LIB)(update.o):	$(CPUMAKEF)
$(LIB)(ns.o):		$(CPUMAKEF)

$(LIB)(callf77.o):	callf77.m4 $(M4FILE)
		$(M4) $(M4FILE) callf77.m4  > callf77.c
		$(CC) $(C-FLAGS) -c callf77.c
		$(AR) $(ARFLAGS) $(LIB) callf77.o
		$(RM) callf77.c callf77.o
timefft:	timefft.o
		$(LD) $(RUNLDOPT) -o $@ $^ $(RUNSYSLIBS)
include	Make.dep
endif

