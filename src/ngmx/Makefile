#
# Makefile for ngmx library
#
include 	$(GMXHOME)/src/makef/Makefile.def
DEFS	=
XOBJ	=	xutil.o		xmb.o		x11.o		\
		xdlghi.o   	fgrid.o    	xdlgitem.o 	\
		popup.o 	pulldown.o	xdlg.o
NOBJ	=	manager.o 	buttons.o 	nmol.o 		\
		nleg.o 		dialogs.o	logo.o 		\
		filter.o	molps.o 	$(XOBJ)
POBJ	=	ngmx.o		xrama.o		gmxlogo.o
PROGRAMS=	ngmx 		xrama		highway
SLIBS	=	$(XLIBS) $(SYSLIBS)

include 	$(STDTARGET)

ngmx:		$(NOBJ) ngmx.o
		$(LD) -o $@ $^  $(SLIBS) 

xrama:		xrama.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

gmxlogo:	gmxlogo.o logo.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

highway:	highway.o $(XOBJ)
		$(LD) -o $@ $^	$(SLIBS)

dlg:		dlg.o $(NOBJ)
		$(LD) -o $@ $^	$(SLIBS)

scrollw:	scrollw.o $(NOBJ)
		$(LD) -o $@ $^	$(SLIBS)

showcol:	showcol.o $(NOBJ)
		$(LD) -o $@ $^	$(SLIBS)

include	Make.dep
