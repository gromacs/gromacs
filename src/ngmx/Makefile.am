## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la

# Ngmx requires X11 - nothing is built if it doesn't exist

if USE_X11

bin_PROGRAMS     = ngmx xrama highway
EXTRA_PROGRAMS   = gmxlogo dlg scrollw showcol

endif

XSRC	=	xutil.c		xmb.c		x11.c		\
		xdlghi.c   	fgrid.c    	xdlgitem.c 	\
		popup.c 	pulldown.c	xdlg.c		\
		manager.c 	buttons.c 	nmol.c 		\
		nleg.c 		dialogs.c	logo.c 		\
		filter.c	molps.c	        Xstuff.h	\
		dialogs.h	logo.h		molps.h		\
		nleg.h		nmol.h		pulldown.h	\
		xdlg.h		xdlgitem.h	xutil.h		\
		buttons.h	fgrid.h		manager.h	\
		nener.h		nload.h		popup.h		\
		x11.h		xdlghi.h	xmb.h		

EXTRA_DIST = 	alert.bm 	gromacs.bm 	play.bm 	\
		rewind.bm 	stop_ani.bm 	ff.bm 		\
		info.bm 	rama.bm 	stop.bm


ngmx_SOURCES    = ngmx.c $(XSRC)
xrama_SOURCES   = xrama.c $(XSRC)
highway_SOURCES = highway.c $(XSRC) 
gmxlogo_SOURCES = gmxlogo.c logo.c $(XSRC)
dlg_SOURCES     = dlg.c $(XSRC)
scrollw_SOURCES = scrollw.c $(XSRC)
showcol_SOURCES = showcol.c $(XSRC)

CLEANFILES   = 	*~ \\\#*
