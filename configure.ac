#######################################################################
# Process this file wth autoconf to produce a configure script.
#######################################################################
 
AC_PREREQ(2.50)
AC_INIT(GROMACS, 3.3_beta_20050823, gmx-users@gromacs.org)
AC_CONFIG_SRCDIR(src/gmxlib/3dview.c)
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST

AM_INIT_AUTOMAKE(gromacs, 3.3_beta_20050823)
AC_PREFIX_DEFAULT(/usr/local/gromacs)

AM_CONFIG_HEADER(src/config.h)
dnl This is the version info according to the libtool versioning system.
dnl It does *not* correspond to the release number.
SHARED_VERSION_INFO="3:0:0"
AC_SUBST(SHARED_VERSION_INFO)


#######################################################################
# Simple options and makefile variables
#######################################################################


### Single/Double
AC_ARG_ENABLE(float,
 [  --disable-float               use double instead of single precision],, enable_float=yes)

AC_ARG_ENABLE(double,
 [  --enable-double               same effect as --disable-float],, enable_double=no)
if test "$enable_double" = "yes"; then
   enable_float=no;
fi
if test "$enable_float" = "no"; then
  AC_DEFINE(GMX_DOUBLE,,[Compile in double precision])
fi
AM_CONDITIONAL([GMX_DOUBLE],[test "$enable_float" = "no"])



### Fortran or not
AC_ARG_ENABLE(fortran,
 [  --enable-fortran              use fortran (default on sgi,ibm,sun,axp)],, 
[case "${host_cpu}-${host_os}" in
  sparc*-solaris* | alpha*-* | rs6000*-aix* | powerpc*-aix* | mips*-irix*) enable_fortran=yes ;;
  *) enable_fortran=no ;;
 esac])
AM_CONDITIONAL([GMX_FORTRAN],[test "$enable_fortran" = "yes"])
if test "$enable_fortran" = "yes"; then
  AC_DEFINE(GMX_FORTRAN,,[Use Fortran innerloops])
fi
# always use CC for linking:
AC_SUBST(F77LINK,"\$(LINK)")


### MPI or not
AC_ARG_ENABLE(mpi,
 [  --enable-mpi                  compile for parallel runs using MPI],,enable_mpi=no)

### Threads
AC_ARG_ENABLE(threads,
 [  --disable-threads             don't try to use multithreading],,enable_threads=no)
AM_CONDITIONAL([GMX_THREAD],[test "$enable_threads" = "yes"])


### MPI environment
AC_ARG_ENABLE(mpi-environment,
 [  --enable-mpi-environment=VAR  only start parallel runs when VAR is set],,
[case "${host_cpu}" in
  mips*-irix*) enable_mpi_environment="MPI_ENVIRONMENT" ;;
  *) enable_mpi_environment=no ;;
esac])
if test "$enable_mpi_environment" != "no"; then
  AC_DEFINE_UNQUOTED(GMX_CHECK_MPI_ENV,"$enable_mpi_environment",[If defined, only start MPI runs when this variable is set])
fi

### IA32 assembly code
AC_ARG_ENABLE(ia32_3dnow,    
 [  --disable-ia32-3dnow          don't build 3DNow! assembly loops on ia32],,enable_ia32_3dnow=yes)
case "${host_cpu}" in
   i?86) ;;     
   *) enable_ia32_3dnow=no ;;
esac


### IA32 assembly code
AC_ARG_ENABLE(ia32_sse,    
 [  --disable-ia32-sse            don't build SSE/SSE2 assembly loops on ia32],,enable_ia32_sse=yes)
case "${host_cpu}" in
   i?86) ;;     
   *) enable_ia32_sse=no ;;
esac

### AMD64 assembly code
AC_ARG_ENABLE(x86_64_sse,
 [  --disable-x86-64-sse           don't build SSE assembly loops on X86_64],,enable_x86_64_sse=yes)
case "${host_cpu}" in
   x86_64 | amd64) ;;
   *) enable_x86_64_sse=no ;;
esac


### Altivec on powerpc (motorola)
AC_ARG_ENABLE(ppc_altivec,
 [  --disable-ppc-altivec         don't build Altivec loops on PowerPC],,enable_ppc_altivec=undef)
#We will probably see ppc970 boxes from ibm with altivec support soon, so make it possible to override...
if test "$enable_ppc_altivec" = "undef"; then
  case "${host_cpu}-${host_vendor}" in
    ppc*-ibm* | powerpc*-ibm*) enable_ppc_altivec=no ;;
    ppc*-* | powerpc*-*) enable_ppc_altivec=yes ;; 
    *) enable_ppc_altivec=no ;;
  esac
fi


### ia64 assembly code
AC_ARG_ENABLE(ia64_asm,
 [  --disable-ia64-asm            don't build assembly loops on ia64],,enable_ia64_asm=yes)
case "${host_cpu}" in
   ia64) ;;
   *) enable_ia64_asm=no ;;
esac
if test "$enable_ia64_asm" = "yes"; then
  AC_DEFINE(GMX_IA64_ASM,,[Use ia64 assembly tuned for Itanium2])
  AC_DEFINE(DISABLE_WATER_NLIST,,[Turn off all water neighborlist optimization])
  AC_DEFINE(DISABLE_WATERWATER_NLIST,,[Turn off water-water neighborlist optimization only])
fi



### Optimize for host cpu version
AC_ARG_ENABLE(cpu-optimization,     
 [  --disable-cpu-optimization    no detection or tuning flags for cpu version],, enable_cpu_optimization=yes)



### Do software 1/sqrt(x)
AC_ARG_ENABLE(software-sqrt,     
 [  --disable-software-sqrt       no software 1/sqrt (disabled on sgi,ibm,ia64)],,
[case "${host_cpu}-${host_os}" in
  mips*-irix* | rs6000*-aix* | powerpc*-aix | ia64*-*) enable_software_sqrt=no ;;
  *) enable_software_sqrt=yes ;;
esac])
if test "$enable_software_sqrt" = "yes"; then
  AC_DEFINE(GMX_SOFTWARE_SQRT,,[Use the GROMACS software 1/sqrt(x)])
fi
AM_CONDITIONAL([GMX_SOFTWARE_SQRT],[test "$enable_software_sqrt" = "yes"])


# The four options below used to improve code scheduling, but the newest CPUs can usually rearrange
# instructions dynamically at runtime, so they might not make much difference now:

### Force prefetching
AC_ARG_ENABLE(prefetch,     
 [  --enable-prefetch-forces      prefetch forces in innerloops],, enable_prefetch=no)
# processing of prefetch choices is done just before the final output stage!
AM_CONDITIONAL([GMX_PREFETCH_FORCES],[test "$enable_prefetch" = "yes"])




AC_ARG_ENABLE(all-static,
 [  --enable-all-static           make completely static binaries],, enable_all_static=no)
# Dont add the -all-static flag until after the compiler test.



### Compile without FFTW - disabling it means you can't use PME or PPPM!
### N.B.: This is only intended as a last resort, or for tests - don't use it!
AC_ARG_WITH(fftw,
 [  --without-fftw          if you turn this off you cannot use PME/PPPM],, with_fftw=yes)
if test "$with_fftw" = "no"; then
  AC_DEFINE(GMX_WITHOUT_FFTW,,[Disable all FFTW (and thus PME/PPPM) support])
fi


### Use external BLAS/LAPACK libraries if the user wants to.
###
AC_ARG_WITH(system_blas,
 [  --with-system-blas      use system BLAS (add to LIBS) instead of built-in],, with_system_blas=no)
AM_CONDITIONAL([GMX_SYSTEM_BLAS], [test "with_system_blas" = "yes"])


AC_ARG_WITH(system_lapack,
 [  --with-system-lapack    use system LAPACK (add to LIBS) instead of built-in],, with_system_lapack=no)
AM_CONDITIONAL([GMX_SYSTEM_LAPACK], [test "with_system_lapack" = "yes"])


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
    cc_names="cc xlc icc gcc"
    f77_names="f77 xlf ifort xlf77 cf77 fl32 g77 fort77 f90 xlf90 pgf77 cf77 fort fort77 pgf90"
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
# Check floating-point format and endian
#############
ACX_FLOAT_FORMAT


# MPI builds on aix break with shared libs, so they are disabled by default.
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO




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
# IBM provides a MASS library with vectorized and optimized math stuff,
# and an xlopt library which partly overlaps this. They must come before
# libm in the link list, thus the test goes after that for libm!
if test "${host_vendor}" = "ibm"; then
  AC_CHECK_LIB(xlopt,main)
  AC_CHECK_LIB(mass,main)  
  vectmass=no
  case "$gmxcpu" in
    power4*)
	# first check for a power4 lib, then power3, power2.
      AC_SEARCH_LIBS(vsrsqrt,massvp4 massvp3 massvp2 massv,[vectmass=yes]) ;;
    power3*)
      AC_SEARCH_LIBS(vsrsqrt,massvp3 massvp2 massv,[vectmass=yes]) ;;
    power2*)
      AC_SEARCH_LIBS(vsrsqrt,massvp2 massv,[vectmass=yes]) ;;
    *)
      AC_CHECK_LIB(massv,vsrsqrt,[vectmass=yes LIBS="$LIBS -lmassv"]) ;;
  esac
  if test "$vectmass" = "yes"; then
    AC_DEFINE(HAVE_LIBMASSV_ANY,,[Use vectorized routines in the IBM MASS library])
    if test "$enable_vectorized_recip" = "not_selected"; then
      AC_MSG_NOTICE([vector MASS lib present - vectorizing 1/x inner loops])
      enable_vectorized_recip="yes"
    fi
    if test "$enable_vectorized_sqrt" = "not_selected"; then
      AC_MSG_NOTICE([vector MASS lib present - vectorizing 1/sqrt inner loops])
      enable_vectorized_sqrt="yes"
    fi    
  else
    echo "******************************************************************************"
    echo "*Didn't find the IBM MASS library. It's not necessary, but can provide 20%   *"
    echo "*performance improvement by loop vectorization. Add the library search path  *"
    echo "*with -L in the LDFLAGS variable. If you cannot find MASS, download it from  *"
    echo "*http://www.research.ibm.com/mass                                            *"
    echo "******************************************************************************"
  fi
fi 

if test "$with_fftw" = "yes"; then
  if test "$enable_float" = "yes"; then
    sizeof_real=4
  else 
    sizeof_real=8
  fi

  ACX_CHECK_FFTW(fftw,$sizeof_real)
  ACX_CHECK_FFTW(rfftw,$sizeof_real)

  if test "$enable_mpi" = "yes"; then
    ACX_CHECK_FFTW(fftw_mpi,$sizeof_real)
    ACX_CHECK_FFTW(rfftw_mpi,$sizeof_real)
  fi

  case ${ac_fftw_savedprefix} in
   d) AC_DEFINE(FFTW_NAME_DFFTW,,[Use the d prefix on fftw includes]) ;;
   s) AC_DEFINE(FFTW_NAME_SFFTW,,[Use the s prefix on fftw includes]) ;;
   *) AC_DEFINE(FFTW_NAME_FFTW,,[Dont use any prefix on fftw includes]) ;;
  esac
else
  echo "*******************************************************************"
  echo "*             WARNING! You have disabled FFTW support!            *"
  echo "* You will not be able to use PME (Particle Mesh Ewald) summation *"
  echo "* in your simulations. Unless you are absolutely sure you won't   *"
  echo "* ever need this we suggest you install FFTW instead of compiling *"
  echo "* a version of GROMACS that will be crippled!                     *"
  echo "*******************************************************************"
fi # end of fftw check

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
AC_FUNC_VPRINTF
AM_WITH_DMALLOC
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(strdup)

#####
# Look for X
AC_PATH_XTRA
# try to find motif headers and libraries
AC_FIND_MOTIF
use_motif=no
if test "$motif_libraries" != "no"; then
  if test "$motif_includes" != "no"; then
    GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} mgmx.lo widget.lo"
    use_motif=yes
  else
    AC_MSG_NOTICE([Both libraries and includes are needed for Motif - disabling.])
    motif_includes=no
  fi
fi
if test "$no_x" != "yes"; then
  LIBS="$X_LIBS $LIBS $X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
  INCLUDES="$INCLUDES $X_CFLAGS $MOTIF_CFLAGS" 
  AM_CONDITIONAL(USE_X11,true)
else
  AM_CONDITIONAL(USE_X11,false)
fi

#########
# Check for libxml2, but it is optional, so dont stop
# if it isnt there.
#AC_CHECK_LIB(xml2,main)
# This should really only be done when the library is found. Doesn't hurt though.
#INCLUDES="$INCLUDES -I/usr/include/libxml2"


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
  if AC_TRY_COMMAND($CC -c conftestasm.s); then
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
# But does it understand the lock cmpxchg instruction (needed for atomicity on SMP)?
  AC_MSG_CHECKING([whether as understands 'lock cmpxchg'])
  if test "$enable_x86_64_sse" = "yes"; then
    cat > conftest.s << EOF
.intel_syntax noprefix 
checkasm:
        lock cmpxchg [[rsi]], ebx 
        emms
        ret
EOF
  else
    cat > conftest.s << EOF
.intel_syntax noprefix
checkasm:
        lock cmpxchg [[esi]], ebx
        emms
        ret
EOF
  fi
    if AC_TRY_COMMAND($CC -c conftest.s); then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no]) 
      AC_MSG_ERROR([Upgrade to binutils >= 2.12.1 (or disable assembly loops).])
    fi
#
# If we are using double precision we need SSE2
  if test "$enable_float" = "no"; then
#   test SSE2 on ia32:
    if test "$enable_ia32_sse" = "yes"; then
      AC_MSG_CHECKING([whether as fully supports ia32 SSE2])
cat > conftest.s << EOF
.intel_syntax noprefix 
checkasm:
        emms
        movupd xmm0,[[esp]]
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC -c conftest.s); then
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
.intel_syntax noprefix 
checkasm:
        emms
        movupd xmm0,[[rsp]]
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC -c conftest.s); then
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
.intel_syntax noprefix 
checkasm:
        emms
        movups xmm0,[[esp]]
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC -c conftest.s); then
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
.intel_syntax noprefix 
checkasm:
        emms
        movups xmm0,[[rsp]]
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC -c conftest.s); then
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
.intel_syntax noprefix 
checkasm:
        emms
        pswapd mm0, mm0
        emms
        ret
EOF
      if AC_TRY_COMMAND($CC -c conftest.s); then
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
[Download the latest devtools from Apple, or gcc-3.3 from ftp.gnu.org])
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

echo ""

























