#
#	@(#) Makefile 1.25 11/19/95
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990, 1991, 1992, Groningen University
#
#	Makefile for gromacs
#
include		$(GMXHOME)/src/Makefile.def
#
SHELL		=/bin/csh
CPUMAKEF	=Makefile.$(CPU)
OBJMAKEF	=Makefile.obj
GROMTAZ		=gromacs.taz
SRCDIRS		= $(GROMLIB) $(notdir $(INCDIR))
EXCLUDE		=core ',*' '*log*' '*~' '\#*\#' '*.o'
IGNORE		=\( $(foreach EXCL,$(EXCLUDE),\( \! -name $(EXCL) \) )
TAZFILES	=Makefile Makefile.def Makefile.std Makefile.obj \
		 $(foreach CPU,$(CPUS),Makefile.$(CPU)) \
                 $(foreach DIR,$(SRCDIRS),`find $(DIR) $(IGNORE) -print `)

obj:	$(ROOT)/obj/Makefile

taz:	source
	tar -cf - $(TAZFILES) | compress -c > $(GROMTAZ)

template:
	echo "#" > $(CPUMAKEF)
	echo "#	@(#) Makefile 1.25 11/19/95" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "#	GROMACS - Groningen Machine for Chemical Simulation" >> $(CPUMAKEF)
	echo "#	Copyright (c) 1990, 1991, 1992, Groningen University" >> $(CPUMAKEF)
	echo "#" >> $(CPUMAKEF)
	echo "#	Makefile for gromacs on a" $(CPU) "processor" >> $(CPUMAKEF)
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
