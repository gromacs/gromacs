## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

include $(srcdir)/../Makefile.inc
# 
# This is necessary for VPATH builds (and thus distcheck) to work, 
# but it can probably included in a nicer way...
INCLUDES = @INCLUDES@ -I$(top_srcdir)/src/include

# produce the gmx library
#

lib_LIBRARIES = libgmxXXX_SUFFIX_XXX.a

libgmxXXX_SUFFIX_XXX_a_SOURCES = \
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

EXTRA_libgmxXXX_SUFFIX_XXX_a_SOURCES = \
	mpiio.c		libnet.c	x86_3dnow.asm	f77_wrappers.c  \
	libxdrf.c	ftocstr.c	dumxdrf.c	mgmx.c		\
	widget.c	widget.h	x86_cpuid.asm	axp_asm.S	\
	x86_sse.asm	

if USE_F77
BUILT_SOURCES = innerf.f
else
BUILT_SOURCES = innerc.c
endif

mpi_obj 	= mpiio.o
libnet_obj 	= libnet.o
inner_f77_obj	= innerf.o	f77_wrappers.o
inner_c_obj	= innerc.o 
xdr_obj		= libxdrf.o ftocstr.o
noxdr_obj	= dumxdrf.o
motif_obj       = mgmx.o widget.o
x86_asm_obj     = x86_cpuid.o
sse_obj         = x86_sse.o
tdn_obj         = x86_3dnow.o 
axp_asm_obj     = axp_asm.o


libgmxXXX_SUFFIX_XXX_a_LIBADD = \
	@PAR_OBJ@ 	@INNER_F77_OBJ@	@INNER_C_OBJ@	\
	@AXP_ASM_OBJ@	@X86_ASM_OBJ@	@SSE_OBJ@	\
	@TDN_OBJ@	@MOTIF_OBJ@	@XDR_OBJ@	


libgmxXXX_SUFFIX_XXX_a_DEPENDENCIES = \
	@PAR_OBJ@ 	@INNER_F77_OBJ@	@INNER_C_OBJ@	\
	@AXP_ASM_OBJ@	@X86_ASM_OBJ@	@SSE_OBJ@	\
	@TDN_OBJ@	@MOTIF_OBJ@	@XDR_OBJ@	


# The inner loops
innerf.f:	mkinl 
		./mkinl	fortran	

innerc.c:	mkinl
		./mkinl c


noinst_PROGRAMS  = mkinl

mkinl_SOURCES 	 = \
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
CLEANFILES     = ${lib_LIBRARIES} *_d.a *_mpi.a *~ \\\#* innerc.c innerf.f

