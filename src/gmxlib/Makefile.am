
## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

lib_LTLIBRARIES = libgmx@LIBSUFFIX@.la

libgmx@LIBSUFFIX@_la_LIBADD = @GMXLIB_LIBOBJS@
libgmx@LIBSUFFIX@_la_DEPENDENCIES = @GMXLIB_LIBOBJS@
libgmx@LIBSUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

libgmx@LIBSUFFIX@_la_SOURCES = \
	3dview.c	atomprop.c	block_tx.c	bondfree.c	\
	buffer.c	calcgrid.c	calch.c		inner.h		\
	confio.c	copyrite.c	disre.c		do_fit.c	\
	enxio.c		ewald_util.c	fatal.c		ffscanf.c	\
	filenm.c	futil.c		gbutil.c	fnbf.c		\
	gmxfio.c	ifunc.c		index.c		cinvsqrtdata.c	\
	crecipdata.c    invblock.c	macros.c	\
	main.c		maths.c		matio.c		memdump.c	\
	mshift.c	mvdata.c	mvxvf.c		names.c		\
	network.c	nrama.c		nrjac.c		nrnb.c		\
	orires.c	\
	pargs.c		pbc.c		pdbio.c		princ.c		\
	rando.c		random.c	rbin.c		rdgroup.c	\
	readinp.c	replace.c	rmpbc.c		shift_util.c	\
	sortwater.c	smalloc.c	stat.c		statutil.c	\
	strdb.c		string2.c	symtab.c	tpxio.c		\
	trnio.c		trxio.c		txtdump.c	typedefs.c	\
	viewit.c	wgms.c		wman.c		writeps.c	\
	xdrd.c		xtcio.c		xvgr.c   	replace.h	\
	detectcpu.c	libxdrf.c	vec.c		

EXTRA_libgmx@LIBSUFFIX@_la_SOURCES = \
	x86_sse.s	x86_sse2.s	x86_3dnow.s  	inner_altivec.c \
	mgmx.c		widget.c	widget.h	x86_cpuid.s	\
	axp_asm.s	f77_wrappers.c	

nodist_EXTRA_libgmx@LIBSUFFIX@_la_SOURCES = innerf.f innerc.c


# Compaq are fscking around with their compiler... some versions require
# .s and some .S. It's checked by configure. COMPAQ: FIX YOUR COMPILER!
if AXP_ASM_CAPITAL_S
axp_asm.lo: $(srcdir)/axp_asm.s
	ln -sf $(srcdir)/axp_asm.s axp_asm.S
	$(LTASCOMPILE) -c axp_asm.S
endif

# The inner loops are made with our program and not included in the dist.
innerf.f: mkinl 
	./mkinl	fortran	
innerc.c: mkinl
	./mkinl c


# mkinl is special - cant use a noinst_PROGRAMS target, since it might have 
# to be compiled with a non-MPI compiler whose files can be executed directly:

EXTRA_DIST 	 =			mkinl_declarations.c	\
		mkinl_outerloop.c	mkinl_innerloop.c	\
		mkinl_invsqrt.c		mkinl_calcdist.c	\
		mkinl_recip.c		mkinl_interactions.c	\
		metacode.c		mkinl.h			\
		mkinl_fortrandata.h	mkinl.c

BUILD_COMPILE = $(CC_FOR_BUILD) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)

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

# link the mpi library to non-mpi names if the latter are not present after install
install-exec-hook:
	libname="libgmx@LIBSUFFIX@"; \
	nompi="`echo $$libname | sed -e 's,_mpi,,'`"; \
	libdir="$(libdir)"; \
	if echo $$libname | grep mpi >/dev/null ; then \
	  (cd $$libdir && test -e $$libname.a -a ! -e $$nompi.a && $(LN_S) $$libname.a $$nompi.a ; exit 0); \
	  (cd $$libdir && test -e $$libname.so -a ! -e $$nompi.so && $(LN_S) $$libname.so $$nompi.so ; exit 0); \
	fi; 

# clean all libtool libraries, since the target names might have changed
CLEANFILES     = *.la *~ \\\#* innerc.c innerf.f mkinl axp_asm.S







