## Process this file with automake to produce Makefile.in
#
# Note that Makefile is generated automatically from Makefile.in,
# which is automatically generated from Makefile.am

include $(srcdir)/../Makefile.inc
# 
# This is necessary for VPATH builds (and thus distcheck) to work, 
# but it can probably be included in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include

# produce the gmx library
#

lib_LTLIBRARIES = libgmxXXX_SUFFIX_XXX.la

libgmxXXX_SUFFIX_XXX_la_SOURCES = \
	3dview.c	atomprop.c	block_tx.c	bondfree.c	\
	buffer.c	calcgrid.c	calch.c		inner.h		\
	confio.c	copyrite.c	disre.c		do_fit.c	\
	enxio.c		ewald_util.c	fatal.c		ffscanf.c	\
	filenm.c	futil.c		gbutil.c	fnbf.c		\
	gmxfio.c	ifunc.c		index.c		cinvsqrtdata.c	\
	crecipdata.c    invblock.c	javaio.c	macros.c	\
	main.c		maths.c		matio.c		memdump.c	\
	mshift.c	mvdata.c	mvxvf.c		names.c		\
	network.c	nrama.c		nrjac.c		nrnb.c		\
	pargs.c		pbc.c		pdbio.c		princ.c		\
	rando.c		random.c	rbin.c		rdgroup.c	\
	readinp.c	replace.c	rmpbc.c		shift_util.c	\
	sortwater.c	smalloc.c	stat.c		statutil.c	\
	strdb.c		string2.c	symtab.c	tpxio.c		\
	trnio.c		trxio.c		txtdump.c	typedefs.c	\
	viewit.c	wgms.c		wman.c		writeps.c	\
	xdrd.c		xtcio.c		xvgr.c   	replace.h	\
	x86_cpu.c	

EXTRA_libgmxXXX_SUFFIX_XXX_la_SOURCES = \
	mpiio.c		libnet.c	x86_3dnow.S	f77_wrappers.c  \
	libxdrf.c	ftocstr.c	dumxdrf.c	mgmx.c		\
	widget.c	widget.h	x86_cpuid.S	axp_asm.s	\
	x86_sse.S	

if USE_FORTRAN
BUILT_SOURCES = innerf.f
else
BUILT_SOURCES = innerc.c
endif

#
# NB: The contents of GMXLIB_COND_OBJ is defined in
#     the main configure.in script file
#

libgmxXXX_SUFFIX_XXX_la_LIBADD = @GMXLIB_COND_OBJ@	

libgmxXXX_SUFFIX_XXX_la_DEPENDENCIES = @GMXLIB_COND_OBJ@	

# The inner loops
innerf.f:	mkinl 
		./mkinl	fortran	

innerc.c:	mkinl
		./mkinl c

# The compaq compiler is really stupid and thinks a .S file is some
# kind of object. Instead it preprocesses the .s files.
# To make this work with both gcc and the compaq compilers, we
# do a workaround to make gcc preprocess the .s file instead:

if GNU_CC
axp_asm.lo:	axp_asm.s
	$(COMPILE) -x assembler-with-cpp -c $(srcdir)/axp_asm.s
endif

# Mkinl is special - we cant use a noinst_PROGRAMS target, since it
# might have to be compiled with a special non-MPI compiler whose files
# can be executed on the build machine (i.e., not for the target host):

EXTRA_DIST 	 = \
	mkinl.c			mkinl_declarations.c	mkinl_outerloop.c  \
	mkinl_innerloop.c	mkinl_calcdist.c	mkinl_invsqrt.c	   \
	mkinl_recip.c		mkinl_interactions.c	metacode.c	   \
	mkinl.h			mkinl_fortrandata.h

metacode.o:	metacode.c
		$(BUILD_COMPILE) -c $(srcdir)/metacode.c
mkinl.o:	mkinl.c
		$(BUILD_COMPILE) -c $(srcdir)/mkinl.c
mkinl_declarations.o: mkinl_declarations.c
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_declarations.c
mkinl_outerloop.o: mkinl_outerloop.c
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_outerloop.c
mkinl_innerloop.o: mkinl_innerloop.c	
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_innerloop.c
mkinl_calcdist.o: mkinl_calcdist.c	
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_calcdist.c	
mkinl_invsqrt.o: mkinl_invsqrt.c
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_invsqrt.c
mkinl_recip.o:	mkinl_recip.c
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_recip.c
mkinl_interactions.o: mkinl_interactions.c
		$(BUILD_COMPILE) -c $(srcdir)/mkinl_interactions.c

MKINL_OBJ	= \
		mkinl.o mkinl_declarations.o mkinl_outerloop.o mkinl_innerloop.o \
		mkinl_calcdist.o mkinl_invsqrt.o mkinl_recip.o mkinl_interactions.o \
		metacode.o

mkinl:		$(MKINL_OBJ)
		$(BUILD_COMPILE) -o $@ $(MKINL_OBJ)

# clean things explicitly, since the target names might have changed
CLEANFILES     = ${lib_LTLIBRARIES} *_d.la *_mpi.la *~ \\\#* innerc.c innerf.f mkinl

