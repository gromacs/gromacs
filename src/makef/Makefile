#!gmake
#
#	@(#) Makefile 1.25 11/19/95
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990, 1991, 1992, Groningen University
#
#	Makefile for gromacs
#
ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include		$(GMXHOME)/src/Makefile.def
#
SHELL		=/bin/csh
CPUMAKEF	=Makefile.$(GMXCPU)
OBJMAKEF	=Makefile.obj
GROMTAZ		=gromacs.taz
SRCDIRS		= $(GROMLIB) $(notdir $(INCDIR))
EXCLUDE		=core ',*' '*log*' '*~' '\#*\#' '*.o'
IGNORE		=\( $(foreach EXCL,$(EXCLUDE),\( \! -name $(EXCL) \) )

obj:	$(ROOT)/obj/Makefile

template:
	echo "#" > $(CPUMAKEF)
	echo "#	@(#) Makefile 1.25 11/19/95" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "#	GROMACS - Groningen Machine for Chemical Simulation" >> $(CPUMAKEF)
	echo "#	Copyright (c) 1990, 1991, 1992, Groningen University" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "#	Makefile for gromacs on a" $(GMXCPU) "processor" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "# compiler, linker & archiver" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "CC		= " >> $(CPUMAKEF)
	echo "CFLAGS		= " >> $(CPUMAKEF)
	echo "LD		= " >> $(CPUMAKEF)
	echo "LDFLAGS		= " >> $(CPUMAKEF)
	echo "SYSLIBS		= " >> $(CPUMAKEF)
	echo "AR		= " >> $(CPUMAKEF)
	echo "RANLIB		= " >> $(CPUMAKEF)

.DEFAULT:	
	$(MAKE) -f $(OBJMAKEF) OBJDIR=$(SRCDIR) $@
endif
