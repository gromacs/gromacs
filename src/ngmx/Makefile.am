## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

include $(srcdir)/../Makefile.inc

# This is necessary for VPATH to work, 
# but it can probably be done in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include
LDFLAGS = @LDFLAGS@ -L${top_builddir}/src/gmxlib -L${top_builddir}/src/mdlib 
LDADD = -lmdXXX_SUFFIX_XXX -lgmxXXX_SUFFIX_XXX

# Ngmx requires X - nothing is built if it doesn't exist

if USE_X11 

bin_PROGRAMS     = 	ngmxXXX_SUFFIX_XXX  	xramaXXX_SUFFIX_XXX	\
			highwayXXX_SUFFIX_XXX

EXTRA_PROGRAMS   = 	gmxlogoXXX_SUFFIX_XXX  dlgXXX_SUFFIX_XXX 	\
			scrollwXXX_SUFFIX_XXX showcolXXX_SUFFIX_XXX

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
		x11.h		xdlghi.h	xmb.h		\
		alert.bm	gromacs.bm	play.bm		\
		rewind.bm	stop_ani.bm	ff.bm		\
		info.bm		rama.bm		stop.bm

ngmxXXX_SUFFIX_XXX_SOURCES    = ngmx.c $(XSRC)

xramaXXX_SUFFIX_XXX_SOURCES   = xrama.c $(XSRC)

highwayXXX_SUFFIX_XXX_SOURCES = highway.c $(XSRC) 

gmxlogoXXX_SUFFIX_XXX_SOURCES = gmxlogo.c logo.c $(XSRC)

dlgXXX_SUFFIX_XXX_SOURCES     = dlg.c $(XSRC)

scrollwXXX_SUFFIX_XXX_SOURCES = scrollw.c $(XSRC)

showcolXXX_SUFFIX_XXX_SOURCES = showcol.c $(XSRC)

endif

# clean things explicitly, since the target names might have changed
CLEANFILES   = 	${bin_PROGRAMS}	${EXTRA_PROGRAMS}	\
		*_mpi *_d *~ \\\#*