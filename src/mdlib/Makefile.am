## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

include $(srcdir)/../Makefile.inc

# This is necessary for VPATH to work, 
# but it can probably be done in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include

lib_LIBRARIES = libmdXXX_SUFFIX_XXX.a

libmdXXX_SUFFIX_XXX_a_SOURCES = \
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

EXTRA_libmdXXX_SUFFIX_XXX_a_SOURCES = \
	cshake.c	csettle.c	clincs.c	\
	fshake.f	fsettle.f	flincs.f	\
	fshaked.f	fsettled.f	flincsd.f

if USE_DOUBLE
  inner_f77_obj	= fshaked.o	fsettled.o	flincsd.o
else
  inner_f77_obj	= fshake.o	fsettle.o	flincs.o
endif

inner_c_obj	= cshake.o	csettle.o	clincs.o


libmdXXX_SUFFIX_XXX_a_LIBADD       = 	@INNER_F77_OBJ@		@INNER_C_OBJ@
libmdXXX_SUFFIX_XXX_a_DEPENDENCIES = 	@INNER_F77_OBJ@		@INNER_C_OBJ@

# clean things explicitly, since the target names might have changed
CLEANFILES     = ${lib_LIBRARIES} *_d.a *_mpi.a *~ \\\#*



