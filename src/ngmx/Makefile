#!gmake
#
# Makefile for ngmx library
#
ifeq ($(GMXHOME),)
all:
	echo "No GROMACS environment, Please source your GMXRC file"
else
include 	$(GMXHOME)/src/makef/Makefile.def

XOBJ	=	xutil.o		xmb.o		x11.o		\
		xdlghi.o   	fgrid.o    	xdlgitem.o 	\
		popup.o 	pulldown.o	xdlg.o		\
		manager.o 	buttons.o 	nmol.o 		\
		nleg.o 		dialogs.o	logo.o 		\
		filter.o	molps.o
ifneq ($(PARALLEL_ONLY_MDRUN),yes)
PROGRAMS=	ngmx 		xrama		highway
endif
SLIBS	=	$(XLIBS) $(SYSLIBS)

include 	$(STDTARGET)

$(GMXHOME)/src/makef/Makefile.def:
	echo "Wrong GROMACS environment, Please source GMXRC in this source tree"

ngmx:		ngmx.o $(XOBJ)
		$(LD) -o $@ $^  $(SLIBS) 

xrama:		xrama.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

gmxlogo:	gmxlogo.o logo.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

highway:	highway.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

dlg:		dlg.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

scrollw:	scrollw.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

showcol:	showcol.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

include	Make.dep
endif
