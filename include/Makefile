#
#	@(#) Makefile 1.85 9/30/97
#
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990, 1991, 1992, Groningen University
#
#	Makefile for the Gromacs include file directory
#
SHELL		= /bin/csh
CPP		= /lib/cpp -I.

MAGIC_H		= magic.h
MAGIC_H_CR	= magic.h.bak

all:		$(MAGIC_H)

$(MAGIC_H):	typedefs.h Makefile
		echo "/* Generated from" `what typedefs.h | awk '{ if (FNR==2) print }'` "*/" > $(MAGIC_H)
		echo "" >> $(MAGIC_H)
		echo '#ifndef GROMACS_MAGIC' >> $(MAGIC_H)
		echo '#define GROMACS_MAGIC '`\
		$(CPP) -P typedefs.h | tr -d "\ \	\\\
" | sum | awk ' { printf("0x%x",$$1) } '` >> $(MAGIC_H)
		echo '#endif' >> $(MAGIC_H)
		copyrgt $(MAGIC_H)
		echo "y" | $(RM) $(MAGIC_H_CR)

clean:		
		$(RM) *~ #*# 

diff:	
		$(RM) diff[12].list
		/bin/ls *.h *.bm > diff1.list
		echo $(HEADERS) | tr "\ " "\012" | sort > diff2.list
		diff diff1.list diff2.list
		$(RM) diff[12].list

typedefs.h:	types/drblock.h\
		types/forcerec.h\
		types/idef.h\
		types/mdatom.h\
		types/nsborder.h\
		types/atoms.h\
		types/energy.h\
		types/graph.h\
		types/ifunc.h\
		types/nblist.h\
		types/simple.h\
		types/block.h\
		types/enums.h\
		types/grid.h\
		types/inputrec.h\
		types/nbslist.h\
		types/symtab.h\
		types/commrec.h\
		types/filenm.h\
		types/group.h\
		types/ishift.h\
		types/nrnb.h\
		types/topology.h\
		types/group.h \
		types/parm.h\
		types/matrix.h\
		types/edsams.h


