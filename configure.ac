#######################################################################
# Process this file wth autoconf to produce a configure script.
#######################################################################
 
AC_PREREQ(2.50)
AC_INIT(gromacs, 4.0.7, [gmx-users@gromacs.org])
AC_CONFIG_SRCDIR(src/gmxlib/3dview.c)
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST

AM_INIT_AUTOMAKE(tar-ustar)
AC_PREFIX_DEFAULT(/usr/local/gromacs)

AM_CONFIG_HEADER(src/config.h)
dnl This is the version info according to the libtool versioning system.
dnl It does *not* correspond to the release number.
SHARED_VERSION_INFO="5:0:0"
AC_SUBST(SHARED_VERSION_INFO)

AC_DISABLE_SHARED


#######################################################################
# Simple options and makefile variables
#######################################################################


### Single/Double
AC_ARG_ENABLE(float,
	     [AC_HELP_STRING([--disable-float],
	                     [Use double instead of single precision])],, enable_float=yes)

AC_ARG_ENABLE(double,
	     [AC_HELP_STRING([--enable-double],
	                     [Same effect as --disable-float])],, enable_double=no)
if test "$enable_double" = "yes"; then
   enable_float=no;
fi
if test "$enable_float" = "no"; then
  AC_DEFINE(GMX_DOUBLE,,[Compile in double precision])
fi
AM_CONDITIONAL([GMX_DOUBLE],[test "$enable_float" = "no"])



### Fortran or not
AC_ARG_ENABLE(fortran,
             [AC_HELP_STRING([--enable-fortran],
	                     [Use Fortran kernels (default on sgi,ibm,sun,axp)])],, 
[case "${host_cpu}-${host_os}" in
  sparc*-solaris* | alpha*-* | rs6000*-aix* | powerpc*-aix* | mips*-irix*) enable_fortran=yes ;;
  *) enable_fortran=no ;;
 esac])
AM_CONDITIONAL([GMX_FORTRAN],[test "$enable_fortran" = "yes"])
if test "$enable_fortran" = "yes"; then
  AC_DEFINE(GMX_FORTRAN,,[Use Fortran kernels])
fi
# always use CC for linking:
AC_SUBST(F77LINK,"\$(LINK)")

### MPI or not
AC_ARG_ENABLE(mpi,
              [AC_HELP_STRING([--enable-mpi],
	                      [Compile for parallel runs using MPI])],,enable_mpi=no)


### MPI environment
AC_ARG_ENABLE(mpi-environment,
	      [AC_HELP_STRING([--enable-mpi-environment=VAR],
	                      [Only start parallel runs when VAR is set])],,
[case "${host_cpu}" in
  mips*-irix*) enable_mpi_environment="MPI_ENVIRONMENT" ;;
  *) enable_mpi_environment=no ;;
esac])
if test "$enable_mpi_environment" != "no"; then
  AC_DEFINE_UNQUOTED(GMX_CHECK_MPI_ENV,"$enable_mpi_environment",[If defined, only start MPI runs when this variable is set])
fi

### MPI_IN_PLACE
AC_ARG_ENABLE(mpi-in-place,
              [AC_HELP_STRING([--disable-mpi-in-place],
                              [Disable MPI_IN_PLACE. For MPI implementations that have broken support for it (some versions of MPICH, for example).])],, enable_mpi_in_place=yes)



AC_ARG_ENABLE(debug,
	      [AC_HELP_STRING([--enable-debug],
                              [Add -g to other command line arguments])],,enable_debug=no)

### IA32 assembly code
AC_ARG_ENABLE(ia32_3dnow,    
              [AC_HELP_STRING([--disable-ia32-3dnow],
                              [Don't build 3DNow! assembly loops on ia32])],,enable_ia32_3dnow=yes)
case "${host_cpu}-${host_vendor}" in
   i?86-*) ;;
   *) enable_ia32_3dnow=no ;;
esac


### IA32 assembly code
AC_ARG_ENABLE(ia32_sse,    
              [AC_HELP_STRING([--disable-ia32-sse],
                              [Don't build SSE/SSE2 assembly loops on ia32])],,enable_ia32_sse=yes)

case "${host_cpu}-${host_vendor}" in
   i?86-*) ;;
   *) enable_ia32_sse=no ;;
esac

### AMD64 assembly code
AC_ARG_ENABLE(x86_64_sse,
              [AC_HELP_STRING([--disable-x86-64-sse],
                              [Don't build SSE assembly loops on X86_64])],,enable_x86_64_sse=yes)
case "${host_cpu}-${host_vendor}" in
   x86_64-* | amd64-*) ;;
   i?86-apple) ;;
   *) enable_x86_64_sse=no ;;
esac


### Altivec on powerpc (motorola)
AC_ARG_ENABLE(ppc_altivec,
              [AC_HELP_STRING([--disable-ppc-altivec],
                              [Don't build Altivec loops on PowerPC])],,enable_ppc_altivec=undef)
# Power3/4/5 running AIX does not support altivec. We identify this with vendor=IBM.
# IBM powerpc 970 boxes run linux, and have vendor=unknown.
if test "$enable_ppc_altivec" = "undef"; then
  case "${host_cpu}-${host_vendor}" in
    ppc*-ibm* | powerpc*-ibm*) enable_ppc_altivec=no ;;
    ppc*-* | powerpc*-*) enable_ppc_altivec=yes ;; 
    *) enable_ppc_altivec=no ;;
  esac
fi



### ia64 assembly code
AC_ARG_ENABLE(ia64_asm,
	      [AC_HELP_STRING([--disable-ia64-asm],
	                      [Don't build assembly loops on ia64])],,enable_ia64_asm=yes)
case "${host_cpu}" in
   ia64) ;;
   *) enable_ia64_asm=no ;;
esac
if test "$enable_ia64_asm" = "yes"; then
  AC_DEFINE(GMX_IA64_ASM,,[Use ia64 assembly tuned for Itanium2])
  AC_DEFINE(DISABLE_WATER_NLIST,,[Turn off all water neighborlist optimization])
  AC_DEFINE(DISABLE_WATERWATER_NLIST,,[Turn off water-water neighborlist optimization only])
fi


AC_ARG_ENABLE(bluegene,
              [AC_HELP_STRING([--enable-bluegene],
                              [Build assembly intrinsics kernels on bluegene])],,enable_bluegene=no)
if test "$enable_bluegene" = "yes"; then
  AC_DEFINE(GMX_BLUEGENE,,[Use assembly intrinsics kernels for BlueGene])
fi

### Optimize for host cpu version
AC_ARG_ENABLE(cpu-optimization,     
              [AC_HELP_STRING([--disable-cpu-optimization],
	                      [No detection or tuning flags for cpu version])],, enable_cpu_optimization=yes)


### Do software 1/sqrt(x)                                                                                                            
AC_ARG_ENABLE(software-sqrt,
              [AC_HELP_STRING([--disable-software-sqrt],
                              [No software 1/sqrt (disabled on sgi,ibm,ia64)])],,
[case "${host_cpu}-${host_os}" in
  mips*-irix* | rs6000*-aix* | powerpc*-aix | ia64*-*) enable_software_sqrt=no ;;
  *) enable_software_sqrt=yes ;;
esac])
# We might turn off software sqrt on PowerPC, so the
# define and makefile conditional is set below that check.

### Do PowerPC 1/sqrt(x)
AC_ARG_ENABLE(ppc-sqrt,     
             [AC_HELP_STRING([--enable-ppc-sqrt],
                             [PowerPC 1/sqrt table lookup. Single precision]
                             [requires Power5 or later. By default we]
                             [use 2 iterations, but if you only run on]
			     [Power6,PPC440,PPC450 or later CPUs you can]
			     [set --enable-ppc-sqrt=1 for better speed.])],,
			     enable_ppc_sqrt=no)
if test "$enable_ppc_sqrt" != "no"; then
   if test "$enable_ppc_sqrt" == "1"; then
     ppc_sqrt_niter=1;
   else 
     ppc_sqrt_niter=2;
   fi
AC_DEFINE_UNQUOTED(GMX_POWERPC_SQRT,$ppc_sqrt_niter,[Use PowerPC 1/sqrt(x) table lookup (number of iterations)])
# Turn off generic software sqrt
   enable_software_sqrt="no";
fi
AM_CONDITIONAL([GMX_POWERPC_SQRT],[test "$enable_ppc_sqrt" != "no"])
AM_CONDITIONAL([GMX_POWERPC_SQRT_SINGLE_ITERATION],[test "$ppc_sqrt_niter" == "1"])


if test "$enable_software_sqrt" = "yes"; then
  AC_DEFINE(GMX_SOFTWARE_SQRT,,[Use the GROMACS software 1/sqrt(x)])
fi
AM_CONDITIONAL([GMX_SOFTWARE_SQRT],[test "$enable_software_sqrt" = "yes"])




# The four options below used to improve code scheduling, but the newest CPUs can usually rearrange
# instructions dynamically at runtime, so they might not make much difference now:

### Force prefetching
AC_ARG_ENABLE(prefetch,     
              [AC_HELP_STRING([--enable-prefetch-forces],
                              [Prefetch forces in innerloops])],, enable_prefetch=no)
# processing of prefetch choices is done just before the final output stage!
AM_CONDITIONAL([GMX_PREFETCH_FORCES],[test "$enable_prefetch" = "yes"])




AC_ARG_ENABLE(all-static,
              [AC_HELP_STRING([--enable-all-static],
                              [Make completely static binaries])],, enable_all_static=no)
# Dont add the -all-static flag until after the compiler test.


# FFT libraries - only use one...
AC_ARG_WITH(fft,
            [AC_HELP_STRING([--with-fft=@<:@fftw3/fftw2/mkl(>=6.0)/fftpack@:>@],
                            [FFT library to use. fftw3 is default, fftpack built in.])],,
                             with_fft=fftw3)

AC_ARG_ENABLE(fftw-measure,
            [AC_HELP_STRING([--disable-fftw-measure],
			    [Dont optimize FFTs. This produces binary exactly]
                            [reproducible runs (as the -reprod flag to mdrun),]
			    [and if you are using FFTW on x86 and specified its]
			    [--enable-sse/sse2 flag (NOT default in FFTW!) the]
			    [kernels used are close-to-optimal anyway.])],,
			     enable_fftw_measure=yes)



### Use external BLAS/LAPACK libraries if the user wants to.
###
AC_ARG_WITH(external_blas,
            [AC_HELP_STRING([--with-external-blas],
                            [Use system BLAS library (add to LIBS). Automatic on OS X.])],,
[
# default action - check if we are on OS X
case "${host_vendor}-${host_os}" in
  apple-darwin*) with_external_blas="yes" ; extra_LDFLAGS="-framework Accelerate" ;;
  *)             with_external_blas="no"  ;;
esac
])
AM_CONDITIONAL([GMX_EXTERNAL_BLAS], [test "$with_external_blas" = "yes"])


AC_ARG_WITH(external_lapack,
            [AC_HELP_STRING([--with-external-lapack],
                            [Use system LAPACK library (add to LIBS). Automatic on OS X.])],,
[
# default action - check if we are on OS X
case "${host_vendor}-${host_os}" in
  apple-darwin*) with_external_lapack="yes" ; extra_LDFLAGS="-framework Accelerate" ;;
  *)             with_external_lapack="no"  ;;
esac
])
AM_CONDITIONAL([GMX_EXTERNAL_LAPACK], [test "$with_external_lapack" = "yes"])




AC_ARG_WITH(qmmm_gaussian,
              [AC_HELP_STRING([--without-qmmm-gaussian],
                              [Use modified Gaussian0x for QM-MM (see website)])],,with_qmmm_gaussian=yes)
if test "$with_qmmm_gaussian" = "yes"; then
  AC_DEFINE(GMX_QMMM_GAUSSIAN,,[Use (modified) Gaussian0x for QM-MM calculations])
  MDLIB_LIBOBJS="$MDLIB_LIBOBJS qm_gaussian.lo"
fi


AC_ARG_WITH(qmmm_gamess,
              [AC_HELP_STRING([--without-qmmm-gamess],
                              [Use modified Gamess-UK for QM-MM (see website)])],,with_qmmm_gamess=no)
if test "$with_qmmm_gamess" = "yes"; then
  AC_DEFINE(GMX_QMMM_GAMESS,,[Use (modified) Gamess-UK for QM-MM calculations])
  MDLIB_LIBOBJS="$MDLIB_LIBOBJS qm_gamess.lo"
fi


AC_ARG_WITH(qmmm_mopac,
              [AC_HELP_STRING([--without-qmmm-mopac],
                              [Use modified Mopac 7 for QM-MM (see website)])],,with_qmmm_mopac=no)
if test "$with_qmmm_mopac" = "yes"; then
  AC_DEFINE(GMX_QMMM_MOPAC,,[Use (modified) Mopac 7 for QM-MM calculations])
  MDLIB_LIBOBJS="$MDLIB_LIBOBJS qm_mopac.lo"
fi
   

############################################################
# Add some debug info: Who is building, and on what machine?
############################################################
AC_DEFINE_UNQUOTED([BUILD_TIME],["`date`"],[Date and time for build])
AC_DEFINE_UNQUOTED([BUILD_USER],["`whoami`@`hostname`"],[User doing build])
AC_DEFINE_UNQUOTED([BUILD_MACHINE],["`uname -srm`"],[Hardware and OS version for build host])



############################################################
# Checks for programs
############################################################

# normally we use the default compiler, but on alpha/linux we try the compaq ones first
case "${host_cpu}-${host_os}" in
  alpha*-linux*)
    cc_names="ccc cc gcc"
    f77_names="fort f77 g77"
    ;;
  *) 
    cc_names="cc icc xlc gcc"
    f77_names="f77 ifort xlf xlf77 cf77 fl32 g77 fort77 f90 xlf90 pgf77 cf77 fort fort77 pgf90"
    ;;
esac
# the (non-mpi) c compiler, which is also used for programs executed during build stage
AC_PROG_CC( $cc_names )
# Check for environment variable CC_FOR_BUILD before overwriting it
if test -z "$CC_FOR_BUILD"; then
  CC_FOR_BUILD=$CC
fi
AC_SUBST(CC_FOR_BUILD) 

if test "$enable_fortran" = "yes"; then 
  AC_PROG_F77( $f77_names )
  if test -z "$F77"; then
    AC_MSG_ERROR([No fortran compiler found])
  fi
fi

AM_PROG_AS( $CC )

# if we are using mpi, also get an MPICC. We cannot set that in the PROG_CC macro
# above, since autoconf checks that the created file can be executed. This would
# fail on platforms where MPI executables can only be run through a batchqueue.

if test "$enable_mpi" = "yes"; then
  AC_CHECK_PROGS(MPICC, mpxlc mpicc mpcc hcc, $CC)
# now change the normal cc to the MPI one - see the comment above.
  CC=$MPICC
  AC_MSG_CHECKING([whether the MPI cc command works]) # be paranoid
  AC_TRY_LINK([#include <mpi.h>],[int argc; char **argv; MPI_Init(&argc,&argv);],
AC_MSG_RESULT([yes]),AC_MSG_ERROR([Cannot compile and link MPI code with $CC]))
  AC_DEFINE(GMX_MPI,,[Make a parallel version of GROMACS using MPI])

# on the cray xt3 we have to tell autoconf that we 
# are actually cross-compiling even if the architecture
# of host and target are the same. so we enforce 
# cross-compilation if __QK_USER__ is defined
  AC_MSG_CHECKING([for catamount])
  AC_TRY_COMPILE([
#if defined __QK_USER__
#else
#error not catamount
#endif
],,[
  AC_MSG_RESULT([yes])
  cross_compiling=yes 
  AC_DEFINE(GMX_CRAY_XT3,[],[Enable special hacks for Cray XT3])
  AC_DEFINE(GMX_NO_SYSTEM,[],[Ignore calls to system(3)])
  AC_DEFINE(GMX_NO_NICE,[],[Ignore calls to nice(3)])
  AC_DEFINE(GMX_BROKEN_CALLOC,[],[Don't use calloc(3)])
  extra_LDFLAGS="$extra_LDFLAGS -lgmalloc -lacml"
],[
  AC_MSG_RESULT([no])
])
# end of "$enable_mpi" = "yes"
fi

if test "$enable_threads" = "yes"; then
   AC_MSG_CHECKING([for pthreads])
   AC_MSG_RESULT([yes])
   AC_DEFINE([GMX_THREAD_PTHREAD],,[Use pthreads for Gromacs multithreading])
#  ACX_PTHREAD([],AC_MSG_ERROR([Cannot compile with pthreads; use --disable-threads])
fi



AH_TEMPLATE([F77_OR_C_FUNC],
            [Set to F77_FUNC(name,NAME) if Fortran used, otherwise 'name' for C.])
AH_TEMPLATE([F77_OR_C_FUNC_],
            [Set to F77_FUNC_(name,NAME) if Fortran used, otherwise 'name' for C.])



if test "$enable_fortran" = "yes"; then
  AC_F77_LIBRARY_LDFLAGS
  AC_F77_WRAPPERS
  AC_DEFINE([F77_OR_C_FUNC(name,NAME)],[F77_FUNC(name,NAME)])
  AC_DEFINE([F77_OR_C_FUNC_(name,NAME)],[F77_FUNC_(name,NAME)])
else
  AC_DEFINE([F77_OR_C_FUNC(name,NAME)],[name])
  AC_DEFINE([F77_OR_C_FUNC_(name,NAME)],[name])
fi



AC_PROG_CPP

if test "$enable_cpu_optimization" = "yes"; then
  ACX_DETECT_GMXCPU
else
  gmxcpu=""
fi
ACX_COMPILER_MAXOPT




#############
# Check integer endian
#############
AC_C_BIGENDIAN([AC_DEFINE([GMX_INTEGER_BIG_ENDIAN],,[Integer byte order is big endian.])],[],[
AC_MSG_ERROR([Cannot determine endian in compiled output])])
# Little endian implied if GMX_INTEGER_BIG_ENDIAN not defined in source.

#############
# Check integer sizes, set to zero if e.g. long long is not supported.
#############
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(long long int)
AC_CHECK_SIZEOF(off_t)

#############
# Make sure size_t can hold pointers.
# On some Nec machines, you reportedly need to specify a flag for size_t to be 64 bits wide.
#############
AC_MSG_CHECKING([that size_t can hold pointers])
AC_TRY_COMPILE([
#include <stdlib.h>
#include <stddef.h>
],[int chk[sizeof(size_t)-sizeof(int *)+1];],AC_MSG_RESULT([yes]),[
   AC_MSG_ERROR([size_t is smaller than int*. Check if your compiler needs a flag to make it 64 bit.])])



#############
# Check floating-point format and endian
#############
ACX_FLOAT_FORMAT


AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
AC_SYS_LARGEFILE
#
# Use our own version of AC_FUNC_FSEEKO from autoconf 2.63 instead
# of the broken one we might get from autoconf 2.61.
ACX_FUNC_FSEEKO_FIXED




############################################################################
# Checks for libraries.
############################################################################
#
# Don't add math library for intel compilers
if $CC -V 2>&1 | grep 'Intel Corporation' > /dev/null 2>&1; then
  AC_MSG_NOTICE([Using built-in math library with intel compiler])
else
  AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR([No math library found]))
fi

#####
# Disable MASS support for now.
#####


case "$with_fft" in
fftw2)
  if test "$enable_float" = "yes"; then
    sizeof_real=4
  else
    sizeof_real=8
  fi
  ACX_CHECK_FFTW2(fftw,$sizeof_real)
  ACX_CHECK_FFTW2(rfftw,$sizeof_real)

  case ${ac_fftw_savedprefix} in
   d) AC_DEFINE(FFTW2_NAME_DFFTW,,[Use the d prefix on fftw2 includes]) ;;
   s) AC_DEFINE(FFTW2_NAME_SFFTW,,[Use the s prefix on fftw2 includes]) ;;
   *) AC_DEFINE(FFTW2_NAME_FFTW,,[Dont use any prefix on fftw2 includes]) ;;
  esac

  MDLIB_LIBOBJS="$MDLIB_LIBOBJS gmx_fft_fftw2.lo"

  ;; #  end of fftw2 check

fftw3) # Much simpler check than fftw2
# Check for header <fftw3.h> AC_CHECK_HEADERS doesnt work, since we must
# use mpicc to get includes - cpp isnt always the same compiler.
AC_MSG_CHECKING([for fftw3.h])
AC_TRY_COMPILE([#include<fftw3.h>],,[
# ok, look for library file too
AC_MSG_RESULT(yes)
if test "$enable_float" = "yes"; then
  AC_CHECK_LIB([fftw3f],main,,AC_MSG_ERROR([Cannot find fftw3f library]))
else
  AC_CHECK_LIB([fftw3],main,,AC_MSG_ERROR([Cannot find fftw3 library]))
fi
  MDLIB_LIBOBJS="$MDLIB_LIBOBJS gmx_fft_fftw3.lo"
],[
# not ok, echo a warning
AC_MSG_ERROR(
[Cannot find the default external FFT library (fftw3).
Other alternatives are 'fftw2', or 'mkl' for Intel MKL.
You are STRONGLY recommended to use one of these - fftw is free.

Use CPPFLAGS and LDFLAGS if the library is installed in a 
non-standard location. (see FAQ at http://www.gromacs.org)
                                                          
If you dont care about performance you can also specify 'fftpack'
to use a slower set of FFTs built into Gromacs. 
(Just install FFTW3 unless you really know what you are doing).
])
])
  ;;

mkl*)
###########
# Intel Math Kernel Library version 6 and later.
##########
AC_MSG_CHECKING([for mkl_dfti.h])
AC_TRY_COMPILE([#include<mkl_dfti.h>],,AC_MSG_RESULT(yes),AC_MSG_ERROR([Cannot find mkl_dfti.h header from Intel Math Kernel Library>=6.0.]\
))
## Check for library
#  AC_CHECK_LIB([guide],main,,AC_MSG_ERROR([Cannot find libguide (Intel MKL)]))
  AC_CHECK_LIB([mkl],DftiComputeForward,,
    AC_MSG_ERROR([Cannot find Intel Math Kernel Library >= 6.0]))

  MDLIB_LIBOBJS="$MDLIB_LIBOBJS gmx_fft_mkl.lo"
  ;;

fftpack | no)
  AC_MSG_WARN([The built-in FFTPACK routines are slow.
Are you sure you don't want to use FFTW? It is free and much faster...])
  MDLIB_LIBOBJS="$MDLIB_LIBOBJS gmx_fft_fftpack.lo"
  ;;
*)
  AC_MSG_ERROR([Unknown FFT library])
  ;;
esac


if test "$enable_fftw_measure" != "yes"; then
   AC_DEFINE(GMX_DISABLE_FFTW_MEASURE,,[Do not optimize FFTW setups (not needed with SSE FFT kernels)])
fi

######
# check for xtc headers
case "${host_os}" in
cygwin* | mingw*)   # use our own xdr stuff for cygwin/mingw
  GMXLIB_LIBOBJS="$GMXLIB_LIBOBJS gmx_system_xdr.lo"
  ;;
*)
  AC_CHECK_HEADERS([rpc/rpc.h])
  AC_CHECK_HEADERS([rpc/xdr.h],,,
[#if HAVE_RPC_RPC_H
# include <rpc/rpc.h>
#endif
])
  AC_CHECK_LIB(nsl,xdr_float)
  AC_TRY_LINK([
#if HAVE_RPC_RPC_H
#include<rpc/rpc.h> 
#else
#error NO rpc.h header
#endif
#if HAVE_RPC_XDR_H
#include<rpc/xdr.h>
#else
#error No xdr.h header
#endif
],[ XDR *xd; float f; xdr_float(xd,&f);],[have_xdr="yes"])

  if test "$have_xdr" != "yes"; then
    AC_DEFINE(GMX_INTERNAL_XDR,,[Use our own instead of system XDR libraries])	
      GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} gmx_system_xdr.lo"
    AC_MSG_WARN([Couldn't find XDR headers and/or libraries - using our own])
  fi
  ;;
esac

#####
# Checks for additional and/or optional functions or libraries.
#AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_TYPE_SIGNAL
AC_TYPE_OFF_T
AC_FUNC_VPRINTF
AM_WITH_DMALLOC
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(strdup)
# check for bool (must come late, since CFLAGS/includes might affect it)
AC_CHECK_TYPES([bool])





#############
# Check integer sizes, set to zero if e.g. long long is not supported.
# (Note: this must come *after* type checks for the corresponding types!)
#############
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(long long int)
AC_CHECK_SIZEOF(off_t)
AC_CHECK_SIZEOF(void*)



# Apple OS X is a bit strange/braindead in that it always reports itself as i386, both for
# 32 and 64 bit builds (even with the 64 bit kernel!).
#
# To get the SSE support right we check for the size of void * that we detected above...
case "${host_vendor}-${host_os}" in
  apple-darwin*)
    enable_ia32_3dnow=no;
    if test $ac_cv_sizeof_voidp = 8; then
      enable_ia32_sse=no;
    else
      enable_x86_64_sse=no;
    fi
    ;;
  *) 
    ;;
esac




#####
# Look for X
AC_PATH_XTRA
if test "$no_x" != "yes"; then
  LIBS="$X_LIBS $LIBS $X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
  INCLUDES="$INCLUDES $X_CFLAGS"
  AM_CONDITIONAL(USE_X11,true)
else
  AM_CONDITIONAL(USE_X11,false)
fi

# GNU Scientific Library or not
AC_ARG_WITH(gsl,
 [  --with-gsl                    Link to the GNU scientific library, enables extra functions in analysis],,with_gsl=no)

######
if test "$with_gsl" = "yes"; then
  ac_save_LDFLAGS="$LDFLAGS"
  AC_CHECK_HEADERS([gsl/gsl_version.h],[
# header found, check for libraries 
  AC_CHECK_LIB(gslcblas,main,[
  LDFLAGS="$ac_save_LDFLAGS -lgslcblas"
# gslcblas library found, check main lib
  AC_CHECK_LIB(gsl,main)])])
#
fi


AC_ARG_WITH(xml,
 [  --without-xml                 do not link to the xml2 library, disallows the use of certain file formats],,with_xml=yes)
 
#########
# Check for libxml2, but it is optional, so dont stop
# if it isnt there.
if test "$with_xml" = "yes"; then
   ac_save_CPPFLAGS="$CPPFLAGS"
   CPPFLAGS="$ac_save_CPPFLAGS -I/usr/include/libxml2"
   AC_CHECK_HEADERS([libxml/parser.h],AC_CHECK_LIB(xml2,main))
   if test "$ac_cv_lib_xml2_main" = "yes"; then
     INCLUDES="$INCLUDES -I/usr/include/libxml2"
   fi
# restore CPPFLAGS
   CPPFLAGS="$ac_save_CPPFLAGS"
fi

#### 
# And collect the f77 libs to the linker - on Solaris, f77 stuff must go first,
# and we need -L/usr/ucblib to find libucb... (stupid OS)
# On e.g. ibm the mass libraries must come before -lm... 
case "${host_cpu}-${host_os}" in
  *-solaris*) 
    LIBS="$FLIBS $LIBS" 
    LDFLAGS="$LDFLAGS -L/usr/ucblib"
    ;;
  *) LIBS="$LIBS $FLIBS" ;;
esac


# add the Accelerate framework if enabled above.
# LDFLAGS might have been set in the CPU optimization, so we cant assign it directly before that.
LDFLAGS="$LDFLAGS $extra_LDFLAGS"






########################################################################
# Checks for header files and types
########################################################################

AC_CHECK_HEADERS([limits.h strings.h unistd.h],,AC_MSG_ERROR([Include headers not found]))
AC_CHECK_HEADERS(unistd.h)

#####
# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_TYPE_UID_T
AC_C_INLINE
    


# Test stuff for ia32 and x86_64 assembly, and add objects/libraries.
if test "$enable_ia32_sse" = "yes" -o "$enable_ia32_3dnow" = "yes" -o "$enable_x86_64_sse" = "yes"; then
  AC_MSG_CHECKING([whether your compiler can handle assembly files (*.s)])
#
# Use 64-bit registers (rsp/rbp) on x86_64
#
if test "$enable_x86_64_sse" = "yes"; then
cat > conftestasm.s << EOF
checkasm:
        mov     %rsp, %rbp
        ret
EOF
else
cat > conftestasm.s << EOF
checkasm:
        movl    %esp, %ebp
        ret
EOF
fi
  if AC_TRY_COMMAND($CC $CFLAGS -c conftestasm.s); then
    if test -f conftestasm.o; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([Upgrade your compiler (or disable assembly loops).])
    fi
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Upgrade your compiler (or disable assembly loops).])
  fi

# ok, cc understands assembly.
#
# If we are using double precision we need SSE2
  if test "$enable_float" = "no"; then
#   test SSE2 on ia32:
    if test "$enable_ia32_sse" = "yes"; then
      AC_MSG_CHECKING([whether as fully supports ia32 SSE2])
cat > conftest.s << EOF
checkasm:
        emms
        movapd 32(%esi,%edx,8),%xmm4
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC $CFLAGS -c conftest.s); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no]) 
        AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
      fi
      AC_DEFINE([GMX_IA32_SSE2],,[Double-precision SSE2 instructions on ia32])
    fi
#   test SSE2 on X86_64:
    if test "$enable_x86_64_sse" = "yes"; then
      AC_MSG_CHECKING([whether as fully supports X86_64 SSE2])
cat > conftest.s << EOF
checkasm:
        emms
	movapd 32(%rsi,%rdx,8),%xmm4
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC $CFLAGS -c conftest.s); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no]) 
        AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
      fi
      AC_DEFINE([GMX_X86_64_SSE2],,[Double-precision SSE2 instructions on X86_64])
    fi
  else 
# end of double precision testing, now do single.
    if test "$enable_ia32_sse" = "yes"; then
      AC_MSG_CHECKING([whether as fully supports ia32 SSE])
cat > conftest.s << EOF
checkasm:
        emms
	movaps 32(%esi,%edx,8),%xmm4
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC $CFLAGS -c conftest.s); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no]) 
        AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
      fi
      AC_DEFINE([GMX_IA32_SSE],,[Single-precision SSE instructions on ia32])
    fi
#   test SSE2 on X86_64:
    if test "$enable_x86_64_sse" = "yes"; then
      AC_MSG_CHECKING([whether as fully supports X86_64 SSE instructions])
cat > conftest.s << EOF
checkasm:
        emms
        movaps 32(%rsi,%rdx,8),%xmm4
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC $CFLAGS -c conftest.s); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no]) 
        AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
      fi
      AC_DEFINE([GMX_X86_64_SSE],,[Single-precision SSE instructions on X86_64])
    fi
    if test "$enable_ia32_3dnow" = "yes"; then
      AC_MSG_CHECKING([whether as fully supports ia32 3DNow! instructions])
cat > conftest.s << EOF
checkasm:
        emms
        pswapd %mm0, %mm0
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC $CFLAGS -c conftest.s); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no]) 
        AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
      fi
      AC_DEFINE([GMX_IA32_3DNOW],,[Single-precision 3DNow! instructions on ia32])
    fi
  fi  
fi
        

# Test stuff for ia64 assembly
#
# Nothing to test right now.


  

#
# Test altivec support.
#
if test "$enable_ppc_altivec" = "yes"; then
  if test "$enable_float" = "no"; then
    AC_MSG_WARN([PowerPC Altivec loops can only be used in single precision - disabling])
    enable_ppc_altivec=no
  else 
    # On FSF (vanilla) GCC we must include altivec.h, and on apple we 
    # should NOT do it after using -faltivec. Go figure...
    case "${host_os}" in
    darwin*)                            
        # do nothing on Apple systems
       ;;
    *)
        AC_CHECK_HEADERS([altivec.h])
       ;;
    esac

    AC_MSG_CHECKING([whether the compiler supports altivec extensions])
    AC_TRY_COMPILE([],[
#ifndef __VEC__
choke_me
#endif
],[
    AC_MSG_RESULT([yes])
###    GMX_LIBOBJS="${GMX_LIBOBJS} inner_altivec.lo"
    AC_DEFINE([GMX_PPC_ALTIVEC],,[Use PowerPC Altivec inner loops])],[
      AC_MSG_RESULT([no])
      AC_MSG_WARN([no altivec support found - On recent PowerPC CPUs]
[you can get a huge performance gain by using an altivec-enabled compiler.]
[On OS X, you need GCC>=3.3, or the IBM compilers. If you are using standard]
[GCC, you need at least version 4.0])
      enable_ppc_altivec=no])
  fi  
fi



##
# Add the appropriate assembly loops 
#
AM_CONDITIONAL([GMX_IA32_SSE],[test "$enable_ia32_sse" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_IA32_SSE2],[test "$enable_ia32_sse" = "yes" -a "$enable_float" = "no"])
AM_CONDITIONAL([GMX_X86_64_SSE],[test "$enable_x86_64_sse" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_X86_64_SSE2],[test "$enable_x86_64_sse" = "yes" -a "$enable_float" = "no"])
AM_CONDITIONAL([GMX_IA32_3DNOW],[test "$enable_ia32_3dnow" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_PPC_ALTIVEC],[test "$enable_ppc_altivec" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_IA64_ASM_SINGLE],[test "$enable_ia64_asm" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_IA64_ASM_DOUBLE],[test "$enable_ia64_asm" = "yes" -a "$enable_float" = "no"])
AM_CONDITIONAL([GMX_BLUEGENE],[test "$enable_bluegene" = "yes"])


#############
# Check integer endian
#############
AC_C_BIGENDIAN([AC_DEFINE([GMX_INTEGER_BIG_ENDIAN],,[Integer byte order is big endian.])],[],[
AC_MSG_ERROR([Cannot determine endian in compiled output])])
# Little endian implied if GMX_INTEGER_BIG_ENDIAN not defined in source.



case "${host_cpu}" in
   i?86) 
         # Check if inline assembly works
         AC_MSG_CHECKING(if the compiler supports gcc inline assembly)
      AC_TRY_LINK([],[float f;int i; asm("fld %1\nfistpl %0\n" : "=m" (*&i) : "f" (f));],[AC_MSG_RESULT(yes)
AC_DEFINE(GMX_X86_GCC_INLINE_ASM,,[Enable x86 gcc inline assembly])],[AC_MSG_RESULT(no)])
         # Check for MSVC inline assembly
     AC_MSG_CHECKING(if the compiler supports MSVC inline assembly)
      AC_TRY_LINK([],[float f;int i; _asm { fld f } ; _asm { fistpl i };],[AC_MSG_RESULT(yes)
AC_DEFINE(GMX_X86_MSVC_INLINE_ASM,,[Enable x86 MSVC inline assembly])],[AC_MSG_RESULT(no)])
        ;;              
   *) ;;
esac






########################################################################
# Final output stage
########################################################################
AS="$CC"
AC_SUBST(ASFLAGS)
AC_SUBST(AS)
AC_SUBST(INCLUDES)          # should be automatic, but doesnt seem to be?
AC_SUBST(GMXLIB_LIBOBJS)
AC_SUBST(MDLIB_LIBOBJS)
AC_SUBST(GMXANA_LIBOBJS)

# To tell libraries apart, we use four versions of the name. For instance, libgmx:
# libgmx		single prec, no mpi.
# libgmx_d		double prec, no mpi.
# libgmx_mpi		single prec, mpi.
# libgmx_mpi_d		double prec, mpi.
# However, the non-suffixed names are linked to the _mpi and/or _d suffixed
# ones upon installation if that is the only version available.

if test "$enable_mpi" = "yes"; then
  LIBSUFFIX="_mpi"
fi
if test "$enable_float" = "no"; then
  LIBSUFFIX="${LIBSUFFIX}_d"
fi
AC_SUBST(LIBSUFFIX) 

# Unless the user has explicitly requested a prefix/suffix/transform, we
# use _d on double precision files. Don't add anything for mpi, but at the
# end we tell the user it is possible and smart to do in some cases!
if test "$program_transform_name" = "s,x,x,"; then
  name_transform_provided=no;
  if test "$enable_float" = "no"; then
    program_transform_name="s,\$\$,_d,"
  fi
else
  name_transform_provided=yes;
fi


# Use a variable for RPM - this way it can be 
# overridden with make RPM=rpm3 rpm
RPM=rpm
AC_SUBST(RPM)


if test "$enable_all_static" = "yes"; then
  LDFLAGS="$LDFLAGS -all-static"
fi

# we have asm source, so this is necessary for automake 1.6
CCAS=$CC
CCASFLAGS=$CFLAGS
 
AC_SUBST(CCAS)
AC_SUBST(CCASFLAGS)

AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ src/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/gmx_blas/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/gmx_lapack/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia32_3dnow/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia32_sse/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia32_sse2/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_x86_64_sse/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_x86_64_sse2/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ppc_altivec/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia64_single/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia64_double/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_bluegene/Makefile ])
AC_CONFIG_FILES([ include/Makefile ])
AC_CONFIG_FILES([ include/types/Makefile ])
AC_CONFIG_FILES([ src/mdlib/Makefile ])
AC_CONFIG_FILES([ src/kernel/Makefile ])
AC_CONFIG_FILES([ src/tools/Makefile ])
AC_CONFIG_FILES([ src/ngmx/Makefile ])
AC_CONFIG_FILES([ src/contrib/Makefile ])
AC_CONFIG_FILES([ scripts/Makefile ])
AC_CONFIG_FILES([ admin/Makefile ])
AC_CONFIG_FILES([ share/Makefile ])
AC_CONFIG_FILES([ share/tutor/Makefile ])
AC_CONFIG_FILES([ share/tutor/gmxdemo/Makefile ])
AC_CONFIG_FILES([ share/tutor/nmr1/Makefile ])
AC_CONFIG_FILES([ share/tutor/nmr2/Makefile ])
AC_CONFIG_FILES([ share/tutor/water/Makefile ])
AC_CONFIG_FILES([ share/tutor/mixed/Makefile ])
AC_CONFIG_FILES([ share/tutor/methanol/Makefile ])
AC_CONFIG_FILES([ share/tutor/speptide/Makefile ])
AC_CONFIG_FILES([ share/template/Makefile ])
AC_CONFIG_FILES([ share/top/Makefile ])
AC_CONFIG_FILES([ share/html/Makefile ])
AC_CONFIG_FILES([ share/html/images/Makefile ])
AC_CONFIG_FILES([ share/html/online/Makefile ])
AC_CONFIG_FILES([ man/Makefile man/man1/Makefile ])

AC_OUTPUT







#########################################################
# Echo some important info, to avoid stupid mistakes
#

if test "$enable_float" = "no" -a "$name_transform_provided" = "no"; then
  echo ""
  echo "* You are compiling a double precision version of Gromacs -"
  echo "  program names will be suffixed with _d to avoid overwriting single"
  echo "  precision files. You can override it with --program-suffix"
fi

if test "$name_transform_provided" = "no" -a "$enable_mpi" = "yes"; then
  echo ""
  echo "* Seems you are compiling with MPI support. You can install the MPI-"
  echo "  enabled programs with suffixed names to have both MPI and non-MPI"
  echo "  versions. This is useful e.g. on supercomputers where you usually"
  echo "  cannot run MPI-linked programs on the login node."
  echo "  Set a suffix with e.g. --program-suffix=_mpi (or _mpi_d for double)."
  echo "  You only need MPI for mdrun, so if you already have non-MPI stuff"
  echo "  installed you can issue make mdrun; make install-mdrun."
fi

if test "$enable_shared" = "no"; then
  echo ""
  echo "* On most platforms you can save 10X space with dynamic libraries, although"
  echo "  the binaries might be less portable. Why not try --enable-shared ?"
fi  

echo ""

























