#
# Makefile for ngmx library
#
include 	$(GMXHOME)/src/makef/Makefile.def
DEFS	=
LIB	=
XOBJ	=	xutil.o		xmb.o		x11.o		\
		xdlghi.o   	fgrid.o    	xdlgitem.o 	\
		popup.o 	pulldown.o	xdlg.o
OBJ	=	manager.o 	buttons.o 	nmol.o 		\
		nleg.o 		dialogs.o	logo.o 		\
		filter.o	molps.o 	$(XOBJ)
POBJ	=	ngmx.o		xrama.o		gmxlogo.o
PROGRAMS=	ngmx 		xrama		highway
SLIBS	=	$(XLIBS) $(SYSLIBS)

include 	$(STDTARGET)

ngmx:		$(OBJ) ngmx.o
		$(LD) -o $@ $^  $(SLIBS) 

xrama:		xrama.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

gmxlogo:	gmxlogo.o logo.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

highway:	highway.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

dlg:		dlg.o $(OBJ)
		$(LD) -o $@ $^	$(SLIBS)

scrollw:	scrollw.o $(OBJ)
		$(LD) -o $@ $^	$(SLIBS)

showcol:	showcol.o $(OBJ)
		$(LD) -o $@ $^	$(SLIBS)

include	Make.dep
