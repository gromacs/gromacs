## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la

# Ngmx requires X11 - nothing is built if it doesn't exist

if USE_X11

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\" $(X_CFLAGS)
AM_LDFLAGS = $(X_LIBS) $(X_PRE_LIBS) -lX11 $(X_EXTRA_LIBS)
bin_PROGRAMS     = ngmx g_xrama 
EXTRA_PROGRAMS   = test_ngmx_dialog  g_showcol

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
g_xrama_SOURCES   = g_xrama.c $(XSRC)
g_showcol_SOURCES = g_showcol.c $(XSRC)
test_ngmx_dialog_SOURCES     = test_ngmx_dialog.c $(XSRC)

CLEANFILES   = 	*~ \\\#*
