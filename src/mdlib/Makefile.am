## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

include $(srcdir)/../Makefile.inc

# This is necessary for VPATH to work, 
# but it can probably be done in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include

# Dependencies and extra objects should come before the target definition
libmdXXX_SUFFIX_XXX_la_LIBADD         = @MDLIB_COND_OBJ@	
libmdXXX_SUFFIX_XXX_la_DEPENDENCIES   = @MDLIB_COND_OBJ@	



# Finally, the library definition

lib_LTLIBRARIES = libmdXXX_SUFFIX_XXX.la

libmdXXX_SUFFIX_XXX_la_SOURCES = \
	calcmu.c	calcvir.c			\
	congrad.c	constr.c	coupling.c	\
	dummies.c	ebin.c		edsam.c		\
	ewald.c         fftgrid.c	force.c  	\
	ghat.c		init.c		md.c		\
	mdatom.c	mdebin.c	nm.c		\
	ns.c		nsb.c		nsgrid.c	\
	pme.c		poisson.c	pppm.c		\
	psgather.c	pull.c 		pullinit.c	\
	pullio.c	pullutil.c	pssolve.c	\
	psspread.c	runner.c	shakef.c	\
	sim_util.c	splittop.c	steep.c		\
	tables.c	tgroup.c	update.c	\
	vcm.c		wnblist.c	poisson.h	\
	splittop.h	wnblist.h

EXTRA_libmdXXX_SUFFIX_XXX_la_SOURCES = \
	cshake.c	csettle.c	clincs.c	\
	fshake.f	fsettle.f	flincs.f	\
	fshaked.f	fsettled.f	flincsd.f

# clean things explicitly, since the target names might have changed
CLEANFILES     = ${lib_LIBRARIES} *_d.a *_mpi.a *~ \\\#*



