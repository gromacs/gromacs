#######################################################################
# Process this file wth autoconf to produce a configure script.
#######################################################################

AC_PREREQ(2.50)
AC_INIT(GROMACS, 3.1, gmx-users@gromacs.org)
AC_CONFIG_SRCDIR(src/gmxlib/3dview.c)
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST

AM_INIT_AUTOMAKE(gromacs, 3.1)
AC_PREFIX_DEFAULT(/usr/local/gromacs)
AM_CONFIG_HEADER(src/config.h)
dnl This is the version info according to the libtool versioning system.
dnl It does *not* correspond to the release number.
SHARED_VERSION_INFO="2:0:0"
AC_SUBST(SHARED_VERSION_INFO)


#######################################################################
# Simple options and makefile variables
#######################################################################


### Single/Double
AC_ARG_ENABLE(float,
 [  --disable-float               use double instead of single precision],, enable_float=yes)
if test "$enable_float" = "no"; then
  AC_DEFINE(DOUBLE,,[Compile in double precision])
  AM_CONDITIONAL(DOUBLE,true)
else
  AM_CONDITIONAL(DOUBLE,false)
fi

### Fortran or not
AC_ARG_ENABLE(fortran,
 [  --enable-fortran              use fortran (default on sgi,ibm,sun,axp)],, 
[case "${host_cpu}-${host_os}" in
  sparc*-solaris* | alpha*-* | rs6000*-aix* | powerpc*-aix* | mips*-irix*) enable_fortran=yes ;;
  *) enable_fortran=no ;;
 esac])
if test "$enable_fortran" = "yes"; then
  AC_DEFINE(USE_FORTRAN,,[Use Fortran for innerloops and some other core stuff])
  GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} innerf.lo f77_wrappers.lo"
  if test "$enable_float" = "yes"; then
    MDLIB_LIBOBJS="${MDLIB_LIBOBJS} flincs.lo fsettle.lo fshake.lo"
  else
    MDLIB_LIBOBJS="${MDLIB_LIBOBJS} flincsd.lo fsettled.lo fshaked.lo"
  fi
else
  GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} innerc.lo"
  MDLIB_LIBOBJS="${MDLIB_LIBOBJS} clincs.lo csettle.lo"
fi
# always use CC for linking:
AC_SUBST(F77LINK,"\$(LINK)")


### MPI or not
AC_ARG_ENABLE(mpi,
 [  --enable-mpi                  compile for parallel runs using MPI],,enable_mpi=no)


### MPI environment
AC_ARG_ENABLE(mpi-environment,
 [  --enable-mpi-environment=VAR  only start parallel runs when VAR is set],,
[case "${host_cpu}" in
  mips*-irix*) enable_mpi_environment="MPI_ENVIRONMENT" ;;
  *) enable_mpi_environment=no ;;
esac])
if test "$enable_mpi_environment" != "no"; then
  AC_DEFINE_UNQUOTED(CHECK_MPI_ENV,"$enable_mpi_environment",[If defined, only start MPI runs when this variable is set])
fi


### X86 assembly code
AC_ARG_ENABLE(x86_asm,     
 [  --disable-x86-asm             don't build assembly loops on x86],,enable_x86_asm=yes)
case "${host_cpu}" in
   i?86) ;;	
   *) enable_x86_asm=no ;;
esac

### Altivec on powerpc (motorola)
AC_ARG_ENABLE(ppc_altivec,
 [  --disable-ppc-altivec         don't build Altivec loops on PowerPC],,enable_ppc_altivec=yes)
case "${host_cpu}-${host_vendor}" in
   ppc*-ibm* | powerpc*-ibm*) enable_ppc_altivec=no ;;
   ppc*-* | powerpc*-*) ;; 
   *) enable_ppc_altivec=no ;;
esac

AC_ARG_ENABLE(axp_asm,
 [  --enable-axp-asm              include vectorized sqrt routines for alpha],,enable_axp_asm=no)


### Optimize for host cpu version
AC_ARG_ENABLE(cpu-optimization,     
 [  --disable-cpu-optimization    no detection or tuning flags for cpu version],, enable_cpu_optimization=yes)


### Vector machine inner loops
AC_ARG_ENABLE(vector,     
 [  --enable-vector               create inner loops for a vector machine],, enable_vector=no)
if test "$enable_vector" = "yes"; then
  AC_DEFINE(USE_VECTOR,,[Optimize for a vector architecture])
fi


### Simplewater hack for SGI - disabled unrolling is the same as simplewater=yes
AC_ARG_ENABLE(waterloop-unrolling,     
 [  --disable-waterloop-unrolling expand the loops instead (hack for sgi)],, 
[case "${host_cpu}-${host_os}" in
  mips*-irix*) enable_waterloop_unrolling=yes ;;
  *) enable_waterloop_unrolling=no ;;
esac])
if test "$enable_waterloops_unrolling" = "no"; then
  AC_DEFINE(SIMPLEWATER,,[Expand the water loops to three inner loops])
fi


### Turn off the special water-water loops on SGI...
AC_ARG_ENABLE(waterwater-loops,     
 [  --disable-waterwater-loops    turn off double unrolled loops (hack for sgi)],,
[case "${host_cpu}-${host_os}" in
  mips*-irix*) enable_waterwater_loops=no ;;
  *) enable_waterwater_loops=yes ;;
esac])
if test "$enable_waterwater_loops" = "no"; then
  AC_DEFINE(DISABLE_WATERWATER_LOOPS,,[Don't use special loops for water-water interactions])
fi


### Do software 1/x
AC_ARG_ENABLE(software-recip,     
 [  --enable-software-recip       perform software 1/x],, enable_software_recip=no)
if test "$enable_software_recip" = "yes"; then
  AC_DEFINE(SOFTWARE_RECIP,,[Use the GROMACS software 1/x])
fi


### Do software 1/sqrt(x)
AC_ARG_ENABLE(software-sqrt,     
 [  --disable-software-sqrt       no software 1/sqrt (disabled on sgi,rs6000)],,
[case "${host_cpu}-${host_os}" in
  mips*-irix* | rs6000*-aix*) enable_software_sqrt=no ;;
  *) enable_software_sqrt=yes ;;
esac])
if test "$enable_software_sqrt" = "yes"; then
  AC_DEFINE(SOFTWARE_SQRT,,[Use the GROMACS software 1/sqrt(x)])
fi


### Vectorize 1/x
AC_ARG_ENABLE(vectorized-recip,     
 [  --enable-vectorized-recip     vectorize LJ-only inner loops],, enable_vectorized_recip=no)
if test "$enable_vectorized_recip" = "yes"; then
  AC_DEFINE(VECTORIZE_RECIP,,[Vectorize the reciprocal calculation])
else 
  if test "$enable_vectorized_recip" != "no"; then
    AC_MSG_ERROR([Vectorized reciprocal is only meaningful for normal LJ-only loops])
  fi
fi

### Vectorize 1/sqrt(x)
AC_ARG_ENABLE(vectorized-sqrt,     
 [  --enable-vectorized-sqrt=A,B  vectorize selected inner loops (default none)
                                choices: normal,solvent,water,water-water],, enable_vectorized_sqrt=not_selected)
# processing of vectorized_sqrt choices is done just before the final output stage!


# The four options below used to improve code scheduling, but the newest CPUs can usually rearrange
# instructions dynamically at runtime, so they might not make much difference now:

### Coordinate prefetching
AC_ARG_ENABLE(prefetch-x,     
 [  --enable-prefetch-x=A,B       prefetch coordinates in selected innerloops],, enable_prefetch_x=no)
# processing of prefetch-x choices is done just before the final output stage!


### Force prefetching
AC_ARG_ENABLE(prefetch-f,     
 [  --enable-prefetch-f=A,B       prefetch forces in selected innerloops],, enable_prefetch_f=no)
# processing of prefetch-f choices is done just before the final output stage!


### First load all coords in a water-water iteration, then square all pairs, then add the sums to get r2.
AC_ARG_ENABLE(hide-square-latency,     
 [  --enable-hide-square-latency  load coords to cache in advance],, enable_hide_square_latency=no)
if test "$enable_hide_square_latency" = "yes"; then
  AC_DEFINE(DECREASE_SQUARE_LATENCY,,[Try to get coordinates to cache before using them])
fi


### Do the table lookups for all interactions between two water molecules before using the results
AC_ARG_ENABLE(hide-table-latency,     
 [  --enable-hide-table-latency   load table data to cache in advance],, enable_hide_table_latency=no)
if test "$enable_hide_table_latency" = "yes"; then
  AC_DEFINE(DECREASE_LOOKUP_LATENCY,,[Try to get table data to cache before using it])
fi



### Compile without FFTW - disabling it means you can't use PME or PPPM!
### N.B.: This is only intended as a last resort, or for tests - don't use it!
AC_ARG_WITH(fftw,
 [  --without-fftw          if you turn this off you cannot use PME/PPPM],, with_fftw=yes)
if test "$with_fftw" = "no"; then
  AC_DEFINE(WITHOUT_FFTW,,[Disable all FFTW (and thus PME/PPPM) support])
fi




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
    cc_names="cc xlc gcc"
    f77_names="f77 xlf xlf77 cf77 fl32 g77 fort77 f90 xlf90 pgf77 cf77 fort fort77 pgf90"
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

# This is a libtool hack. We never use c++, but libtool supports it and
# thus there are macros that require it. For some reason, these macros are
# always expanded by the libtool scripts, unless the are already present 
# in the file. Thus, we include them here but never run them...
if test 1 = 0; then
  AC_PROG_CXX
  AC_PROG_CXXCPP
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
  AC_DEFINE(USE_MPI,,[Make a parallel version of GROMACS using MPI])
fi

if test "$enable_fortran" = "yes"; then
  AC_F77_LIBRARY_LDFLAGS
  AC_F77_WRAPPERS
fi

AC_PROG_CPP

if test "$enable_cpu_optimization" = "yes"; then
  ACX_DETECT_GMXCPU
else
  gmxcpu=""
fi
ACX_COMPILER_MAXOPT

# Check x86 asm prerequisites and the capabilities of as
if test "$enable_x86_asm" = "yes"; then
  AC_MSG_CHECKING([whether your compiler can handle assembly files (*.s)])
cat > conftestasm.s << EOF
checkasm:
	movl	%esp, %ebp
	ret
EOF
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

# ok, gcc understands assembly, but can it handle SSE/3DNow/SSE2 instructions?

  if test "$enable_float" = "no"; then
    AC_MSG_CHECKING([whether as fully supports intel syntax SSE2])
cat > conftest.s << EOF
.intel_syntax noprefix 
checkasm:
	emms
	movupd xmm0,[[checkasm]]
	emms
	ret
EOF
    if AC_TRY_COMMAND($CC -c conftest.s); then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no]) 
      AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
    fi
    GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} x86_cpuid.lo x86_sse2.lo"
    AC_DEFINE([USE_X86_SSE2],,[Use X86 SSE2 multimedia assembly code])
  else 
    AC_MSG_CHECKING([whether as fully supports intel syntax SSE & 3DNow])
cat > conftest.s << EOF
.intel_syntax noprefix 
checkasm:
	emms
	pswapd mm0,mm0
	movups xmm0,[[checkasm]]
	emms
	ret
EOF
    if AC_TRY_COMMAND($CC -c conftest.s); then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no]) 
      AC_MSG_ERROR([Upgrade to a more recent binutils (or disable assembly loops).])
    fi
    GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} x86_cpuid.lo x86_sse.lo x86_3dnow.lo"
    AC_DEFINE([USE_X86_SSE_AND_3DNOW],,[Use X86 SSE/3DNow multimedia assembly code])
  fi
fi

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
    GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} inner_altivec.lo"
    AC_DEFINE([USE_PPC_ALTIVEC],,[Use PowerPC Altivec inner loops])],[
      AC_MSG_RESULT([no])
      AC_MSG_WARN([no altivec support found - On G4 and later Motorola PowerPC]
[CPUs you can get a HUGE performance gain by using an altivec-enabled compiler.]
[Check out www.altivec.org, or use the development version of gcc-3.1])
      enable_ppc_altivec=no])
  fi  
  # the gcc instruction scheduler destroys our handcoded altivec.
  if test "$enable_ppc_altivec" = "yes"; then
    ACX_CHECK_CC_FLAGS(-fno-schedule-insns,fno_schedule_insns,CFLAGS="$CFLAGS -fno-schedule-insns")
  fi
fi

if test "$enable_axp_asm" = "yes"; then
  case "${host_cpu}" in
    alpha*) 
      AC_DEFINE([USE_AXP_ASM],,[Use Alpha assembly for vectorized 1/sqrt])
      cat > conftest.s << EOF
#ifdef __ELF__
.section .rodata
#else			
	 .rdata
#endif	
EOF
      if AC_TRY_COMMAND($CC -c conftest.s); then
        AM_CONDITIONAL(AXP_ASM_CAPITAL_S,false)
      else
        AM_CONDITIONAL(AXP_ASM_CAPITAL_S,true)
      fi
      GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} axp_asm.lo"
       ;;
    *) ;;
  esac
else
  AM_CONDITIONAL(AXP_ASM_CAPITAL_S,false)
fi

AC_PATH_PROG(IDENT,ident,no)
if test "$IDENT" != "no"; then
  # seems as if we have the ident program, but does the
  # compiler support it?
  AC_MSG_CHECKING([whether the compiler supports ident])	
  AC_TRY_CPP([#ident  "@(#) file.h 1.1 12/16/92"],[
    AC_MSG_RESULT([yes]) 
    AC_DEFINE(HAVE_IDENT,,[Does this machine have the ident command])],
    AC_MSG_RESULT([no]))
fi

# MPI builds on aix break with shared libs.
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
AC_CHECK_HEADERS(rpc/rpc.h rpc/xdr.h,
[# check for xtc libs if we found headers - on solaris the xdr stuff is in libnsl
AC_CHECK_LIB(nsl,xdr_float)
AC_TRY_LINK([#include<rpc/rpc.h> 
#include<rpc/xdr.h>],[ XDR *xd; float f; xdr_float(xd,&f);],
[have_xdr=yes],AC_MSG_WARN([RPC/XDR libraries not found - disabling xtc]))
],AC_MSG_WARN([RPC/XDR include headers not found - disabling xtc]))

if test "$have_xdr" = "yes"; then
  AC_DEFINE(USE_XDR,,[Use portable trajectory routines])
else
  echo "*****************************************************************************"
  echo "*WARNING - XDR not found. This cripples GROMACS significantly - you cannot  *"
  echo "*read/write portable trajectories. No choice on windows, but locate the RPC *"
  echo "*header files and libraries if you run UNIX - they do exist!                *"
  echo "*****************************************************************************"
fi

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
    





case "${host_cpu}" in
   i?86) 
      # disable x86trunc if inlined assembly doesnt work 
      AC_MSG_CHECKING(if the C compiler supports inline assembly)
      AC_TRY_LINK([],[float f;int i; asm("fld %1\nfistpl %0\n" : "=m" (*&i) : "f" (f));],[AC_MSG_RESULT(yes)
AC_DEFINE(USE_X86TRUNC,,[Disable inline assembly truncation for x86 cpu])],[AC_MSG_RESULT(no)])
        ;;	
   *) enable_x86_asm=no ;;
esac






########################################################################
# Process vectorization and prefetch options for inner loops
########################################################################
# check which loops to vectorize
list_of_vectorized_sqrt=""
if test "$enable_vectorized_sqrt" != "no"; then
  if test "$enable_vectorized_sqrt" = "yes" -o "$enable_vectorized_sqrt" = "all" ; then
    enable_vectorized_sqrt="normal,solvent,water,water-water"
  fi
  if echo $enable_vectorized_sqrt | grep "normal" >/dev/null; then
    list_of_vectorized_sqrt="normal "
    AC_DEFINE(VECTORIZE_INVSQRT,,[Vectorize 1/sqrt(x) in normal loops])
  fi
  if echo $enable_vectorized_sqrt | grep "solvent" >/dev/null; then
    list_of_vectorized_sqrt="${list_of_vectorized_sqrt}solvent "
    AC_DEFINE(VECTORIZE_INVSQRT_S,,[Vectorize 1/sqrt(x) in solvent loops])
  fi
  if echo $enable_vectorized_sqrt | sed 's/water-water//' | grep "water" >/dev/null; then
    list_of_vectorized_sqrt="${list_of_vectorized_sqrt}water "
    AC_DEFINE(VECTORIZE_INVSQRT_W,,[Vectorize 1/sqrt(x) in water loops])
  fi
  if echo $enable_vectorized_sqrt | grep "water-water" >/dev/null; then
    list_of_vectorized_sqrt="${list_of_vectorized_sqrt}water-water "
    AC_DEFINE(VECTORIZE_INVSQRT_WW,,[Vectorize 1/sqrt(x) in water-water loops])
  fi
fi
if test -z "$list_of_vectorized_sqrt"; then
  list_of_vectorized_sqrt="no"
fi


# check loops to prefetch coordinates in
list_of_prefetch_x=""
if test "$enable_prefetch_x" != "no"; then
  if test "$enable_prefetch_x" = "yes" -o "$enable_prefetch_x" = "all"; then
    enable_prefetch_x="normal,solvent,water,water-water"
  fi
  if echo $enable_prefetch_x | grep "normal" >/dev/null; then
    list_of_prefetch_x="normal "
    AC_DEFINE(PREFETCH_X,,[Prefetch coordinates in normal loops])
  fi
  if echo $enable_prefetch_x | grep "solvent" >/dev/null; then
    list_of_prefetch_x="${list_of_prefetch_x}solvent "
    AC_DEFINE(PREFETCH_X_S,,[Prefetch coordinates in solvent loops])
  fi
  if echo $enable_prefetch_x | sed 's/water-water//' | grep "water" >/dev/null; then
    list_of_prefetch_x="${list_of_prefetch_x}water "
    AC_DEFINE(PREFETCH_X_W,,[Prefetch coordinates in water loops])
  fi
  if echo $enable_prefetch_x | grep "water-water" >/dev/null; then
    list_of_prefetch_x="${list_of_prefetch_x}water-water "
    AC_DEFINE(PREFETCH_X_WW,,[Prefetch coordinates in water-water loops])
  fi
fi
if test -z "$list_of_prefetch_x"; then
  list_of_prefetch_x="no"
fi


# check loops to prefetch forces in
list_of_prefetch_f=""
if test "$enable_prefetch_f" != "no"; then
  if test "$enable_prefetch_f" = "yes" -o "$enable_prefetch_f" = "all"; then
    enable_prefetch_f="normal,solvent,water,water-water"
  fi
  if echo $enable_prefetch_f | grep "normal" >/dev/null; then
    list_of_prefetch_f="normal "
    AC_DEFINE(PREFETCH_F,,[Prefetch coordinates in normal loops])
  fi
  if echo $enable_prefetch_f | grep "solvent" >/dev/null; then
    list_of_prefetch_f="${list_of_prefetch_f}solvent "
    AC_DEFINE(PREFETCH_F_S,,[Prefetch coordinates in solvent loops])
  fi
  if echo $enable_prefetch_f | sed 's/water-water//' | grep "water"; then
    list_of_prefetch_f="${list_of_prefetch_f}water "
    AC_DEFINE(PREFETCH_F_W,,[Prefetch coordinates in water loops])
  fi
  if echo $enable_prefetch_f | grep "water-water" >/dev/null; then
    list_of_prefetch_f="${list_of_prefetch_f}water-water "
    AC_DEFINE(PREFETCH_F_WW,,[Prefetch coordinates in water-water loops])
  fi
fi
if test -z "$list_of_prefetch_f"; then
  list_of_prefetch_f="no"
fi








########################################################################
# Final output stage
########################################################################
AS="$CC"
AC_SUBST(ASFLAGS)
AC_SUBST(AS)
AC_SUBST(INCLUDES)          # should be automatic, but doesnt seem to be?
AC_SUBST(GMXLIB_LIBOBJS)
AC_SUBST(MDLIB_LIBOBJS)

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

# Set exec-prefix from the architecture and cpu, unless it was 
# explicitly given by the user.
if test "$exec_prefix" = "NONE"; then
  arch_exec_prefix=yes;	
  if test -n "$gmxcpu"; then
    exec_prefix="\${prefix}/${host}/${gmxcpu}"
  else
    exec_prefix="\${prefix}/${host}"
  fi
fi
AC_SUBST(exec_prefix)



#If the prefix doesn't contain "gromacs" or "gmx" and the datadir
#wasn't explicitly specified we add gromacs to the end of datadir:
if test "$datadir" = "\${prefix}/share" -a -z "`echo ${prefix} | grep gromacs`"; then
  if test "${prefix}" != "NONE" -a -z "`echo ${prefix} | grep gromacs`" -a -z "`echo ${prefix} | grep gmx`"; then
    datadir="\${prefix}/share/gromacs"
  fi
fi
AC_SUBST(datadir)

#If the prefix doesn't contain "gromacs" or "gmx" and the includedir
#wasn't explicitly specified we add gromacs to the end of includedir:
if test "$includedir" = "\${prefix}/include" -a -z "`echo ${prefix} | grep gromacs`"; then
  if test "${prefix}" != "NONE" -a -z "`echo ${prefix} | grep gromacs`" -a -z "`echo ${prefix} | grep gmx`"; then
    includedir="\${prefix}/include/gromacs"
  fi
fi
AC_SUBST(includedir)

# automake has started to put the man pages under share, and says it
# should be so according to FHS, but this is still debated so we keep
# the old setup unless the user says different.
if test "$mandir" = "\${prefix}/share/man"; then
  mandir="\${prefix}/man"
  AC_SUBST(mandir)
fi

# Use a variable for RPM - this way it can be 
# overridden with make RPM=rpm3 rpm
RPM=rpm
AC_SUBST(RPM)

AC_CONFIG_FILES([Makefile src/Makefile src/gmxlib/Makefile 
include/Makefile include/types/Makefile
src/mdlib/Makefile src/kernel/Makefile src/tools/Makefile 
src/ngmx/Makefile src/contrib/Makefile scripts/Makefile admin/Makefile 
share/Makefile share/tutor/Makefile share/tutor/gmxdemo/Makefile 
share/tutor/nmr1/Makefile share/tutor/nmr2/Makefile 
share/tutor/water/Makefile share/tutor/speptide/Makefile 
share/template/Makefile share/top/Makefile share/html/Makefile 
share/html/images/Makefile share/html/online/Makefile man/Makefile 
man/man1/Makefile])

AC_OUTPUT 







#########################################################
# Echo some important info, to avoid stupid mistakes
#

echo ""
echo "GROMACS is ready to compile."

eval "show_path=`echo ${exec_prefix} | sed -e s,NONE,${ac_default_prefix},`"
eval "show_datadir=`echo ${datadir} | sed -e s,NONE,${ac_default_prefix},`"
echo ""
echo "* Binaries and libraries for this host will be installed in"
echo "    ${show_path}/bin" 
echo "    ${show_path}/lib" 
echo "  and shared data will placed in"
echo "    ${show_datadir}"
echo "  (You can set this with --prefix, --exec-prefix and --datadir)"


if test "$enable_float" = "no" -a "$name_transform_provided" = "no"; then
  echo ""
  echo "* You are compiling a double precision version of the package -"
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

if test "$enable_shared" = "yes"; then
  echo ""
  echo "* To save space, you might want to use shared libraries."
  echo "  This should work on almost any platform, but since there are"
  echo "  a few exceptions (MPI builds on aix) it is not used by default."
  echo "  Turn it on with --enable-shared."
fi
echo ""

























