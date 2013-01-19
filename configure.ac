#######################################################################
# Process this file wth autoconf to produce a configure script.
#######################################################################
 
AC_PREREQ(2.50)

# The "-dev" suffix is important to keep because it makes possible to distinguish 
# between a build from official release and a build from git release branch on a 
# machine with no git. 
#
# NOTE: when releasing the "-dev" suffix needs to be stripped off!
AC_INIT(gromacs, 4.5.6-dev, [gmx-users@gromacs.org])
AC_CONFIG_SRCDIR(src/gmxlib/3dview.c)
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST

AM_INIT_AUTOMAKE(tar-ustar)
AC_PREFIX_DEFAULT(/usr/local/gromacs)

AM_CONFIG_HEADER(src/config.h)
dnl This is the version info according to the libtool versioning system.
dnl It does *not* correspond to the release number.
SHARED_VERSION_INFO="6:0:0"
AC_SUBST(SHARED_VERSION_INFO)

AC_ENABLE_SHARED
test "$enable_mpi" = "yes" && AC_DISABLE_SHARED
PKG_CFLAGS=""

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
  PKG_CFLAGS="$PKG_CFLAGS -DGMX_DOUBLE"
fi
AM_CONDITIONAL([GMX_DOUBLE],[test "$enable_float" = "no"])



### Fortran or not
AC_ARG_ENABLE(fortran,
             [AC_HELP_STRING([--enable-fortran],
	                     [Use Fortran kernels])],, 
[case "${host_cpu}-${host_os}" in
  sparc*-solaris* | alpha*-* | rs6000*-aix* | powerpc*-aix* | mips*-irix*) enable_fortran=yes ;;
  *) enable_fortran=no ;;
 esac])
# always use CC for linking:
AC_SUBST(F77LINK,"\$(LINK)")

#On some systems GNU libtool's la files are more confusing for libtool than helpful due to the 
#hard coded paths/libs in the files. And for better linking we have pkg-config files (.pc).
AC_ARG_ENABLE(la-files,
	[AS_HELP_STRING([--disable-la-files],[Do NOT install GNU libtool's la files])],,
	[enable_la_files=yes])
AM_CONDITIONAL(NO_LA_FILES,[test "$enable_la_files" = no])




### MPI or not
AC_ARG_ENABLE(mpi,
              [AC_HELP_STRING([--enable-mpi],
	                      [Compile for parallel runs using MPI]
                              [(instead of threads)])],,enable_mpi=no)


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
                              [Disable MPI_IN_PLACE. For MPI implementations]
                              [that have broken support for it ]
                              [(some versions of MPICH, for example).])],, enable_mpi_in_place=yes)



AC_ARG_ENABLE(debug,
	      [AC_HELP_STRING([--enable-debug],
                              [Add -g to other command line arguments])],,enable_debug=no)

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

###################################
## Detect buggy version of gcc
###################################
AC_ARG_ENABLE(gcc41_check,
              [AC_HELP_STRING([--disable-gcc41-check],
                              [Disable the check for buggy gcc 4.1])],,enable_gcc41_check=yes)



### ia64 assembly code
AC_ARG_ENABLE(ia64_asm,
	      [AC_HELP_STRING([--enable-ia64-asm],
	                      [Build assembly loops on ia64])],,enable_ia64_asm=no)
case "${host_cpu}" in
   ia64) ;;
   *) enable_ia64_asm=no ;;
esac


# IBM Power6-specific optimization
AC_ARG_ENABLE(power6,
             [AC_HELP_STRING([--enable-power6],
                             [Use IBM Power6-specific F77 kernels])],,enable_power6=no)
if test "$enable_power6" = "yes"; then
  AC_DEFINE(GMX_POWER6,,[Enable IBM Power6-specific F77 kernels])
fi

AC_ARG_ENABLE(bluegene,
              [AC_HELP_STRING([--enable-bluegene],
                              [Build assembly intrinsics kernels on bluegene])],,enable_bluegene=no)


### Optimize for host cpu version
AC_ARG_ENABLE(cpu-optimization,     
              [AC_HELP_STRING([--disable-cpu-optimization],
	                      [No detection or tuning flags for cpu version])],, enable_cpu_optimization=yes)


### Do software 1/sqrt(x)                                                                                                            
AC_ARG_ENABLE(software-invsqrt,
              [AC_HELP_STRING([--disable-software-invsqrt],
                              [No software 1/sqrt (disabled on sgi,ibm,ia64)])],,
[case "${host_cpu}-${host_os}" in
  mips*-irix* | rs6000*-aix* | powerpc*-aix | powerpc*-none | ia64*-*) enable_software_invsqrt=no ;;
  *) enable_software_invsqrt=yes ;;
esac])
if test "$enable_software_invsqrt" = "yes"; then
  AC_DEFINE(GMX_SOFTWARE_INVSQRT,,[Use the GROMACS software 1/sqrt(x)])
  PKG_CFLAGS="$PKG_CFLAGS -DGMX_SOFTWARE_INVSQRT"
fi
AM_CONDITIONAL([GMX_SOFTWARE_INVSQRT],[test "$enable_software_invsqrt" = "yes"])

AC_ARG_ENABLE(fahcore,
              [AC_HELP_STRING([--enable-fahcore],
                              [Create a library with mdrun functionality])],, enable_fahcore=no)



# The four options below used to improve code scheduling, but the newest CPUs can usually rearrange
# instructions dynamically at runtime, so they might not make much difference now:


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


AC_ARG_ENABLE(threads,
              [AC_HELP_STRING([--disable-threads],    
                              [Disable parallel runs using threads])],,enable_threads=yes)


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
fi


AC_ARG_WITH(qmmm_gamess,
              [AC_HELP_STRING([--without-qmmm-gamess],
                              [Use modified Gamess-UK for QM-MM (see website)])],,with_qmmm_gamess=no)
if test "$with_qmmm_gamess" = "yes"; then
  AC_DEFINE(GMX_QMMM_GAMESS,,[Use (modified) Gamess-UK for QM-MM calculations])
fi


AC_ARG_WITH(qmmm_mopac,
              [AC_HELP_STRING([--without-qmmm-mopac],
                              [Use modified Mopac 7 for QM-MM (see website)])],,with_qmmm_mopac=no)
if test "$with_qmmm_mopac" = "yes"; then
  AC_DEFINE(GMX_QMMM_MOPAC,,[Use (modified) Mopac 7 for QM-MM calculations])
fi


AC_ARG_WITH(qmmm_orca,
              [AC_HELP_STRING([--without-qmmm-orca],
                              [Use ORCA for QM-MM])],,with_qmmm_orca=no)
if test "$with_qmmm_orca"="yes"; then
  AC_DEFINE(GMX_QMMM_ORCA,,[Use ORCA for QM-MM calculations])
fi


AC_ARG_WITH(dlopen,
            [AC_HELP_STRING([--without-dlopen],
                            [do not compile with dlopen, needed to read VMD]
                            [supported file formats])],,with_dlopen=yes)



############################################################
# Add some debug info: Who is building, and on what machine?
############################################################
AC_DEFINE_UNQUOTED([BUILD_TIME],["`date`"],[Date and time for build])
AC_DEFINE_UNQUOTED([BUILD_USER],["`whoami`@`hostname`"],[User doing build])
AC_DEFINE_UNQUOTED([BUILD_MACHINE],["`uname -srm`"],[Hardware and OS version for build host])


############################################################
# Checks for requirements for version information generation
############################################################
generate_version=no
if test -d $srcdir/.git ; then
    AC_CHECK_PROG(have_git, git, yes, no)
    if test "$have_git" = "yes" ; then
        AC_MSG_CHECKING(for a recent enough git)
	if git --git-dir=$srcdir/.git rev-list -n1 --pretty=format:%cd HEAD | grep "Mon\|Tue\|Wed\|Thu\|Fri\|Sat\|Sun" >/dev/null 2>&1; then
	    AC_MSG_RESULT(yes)
            generate_version=yes
            AC_DEFINE(USE_VERSION_H,,[Use generated version info])
        else
            AC_MSG_RESULT(no)
        fi
    fi
    if test "$generate_version" = "no" ; then
        AC_MSG_WARN([No compatible git version found, won't be able to generate proper development version information.])
    fi
fi
AM_CONDITIONAL([USE_VERSION_H],[test "$generate_version" = "yes"])


############################################################
# Checks for programs
############################################################

# normally we use the default compiler, but on alpha/linux we try the compaq ones first
case "${host_cpu}-${host_os}" in
  alpha*-linux*)
    cc_names="ccc cc gcc"
    f77_names="fort f77 gfortran g77"
    ;;
  *) 
    cc_names="cc icc xlc gcc"
    f77_names="f77 ifort xlf xlf77 cf77 fl32 gfortran g77 fort77 f90 xlf90 pgf77 cf77 fort fort77 pgf90"
    ;;
esac
# the (non-mpi) c compiler, which is also used for programs executed during build stage
AC_PROG_CC( $cc_names )
# Check for environment variable CC_FOR_BUILD before overwriting it
if test -z "$CC_FOR_BUILD"; then
  CC_FOR_BUILD=$CC
fi

if test "$enable_gcc41_check" = "yes"; then
    if $CC --version 2>&1 | grep 'gcc.* 4\.1' > /dev/null 2>&1; then
        AC_MSG_ERROR([Uh-oh. Your compilers appears to be GCC version 4.1, which unfortunately produces]
[buggy code at high optimimzation levels. It would be a good idea to update or]
[use a different compiler. If you are ABSOLUTELY sure what you are doing, you]
[can override this check with --disable-gcc41-check.])
    fi
fi



AC_SUBST(CC_FOR_BUILD) 

if test "$enable_fortran" = "yes" -o "$enable_power6" = "yes"; then 
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
  enable_threads="no"
  AC_CHECK_PROGS(MPICC, mpxlc mpicc mpcc hcc, $CC)
# now change the normal cc to the MPI one - see the comment above.
  CC=$MPICC
  CXX=$MPICC
  AC_MSG_CHECKING([whether the MPI cc command works]) # be paranoid
  AC_TRY_LINK([#include <mpi.h>],[int argc; char **argv; MPI_Init(&argc,&argv);],
AC_MSG_RESULT([yes]),AC_MSG_ERROR([Cannot compile and link MPI code with $CC]))
  AC_DEFINE(GMX_MPI,,[Make a parallel version of GROMACS using MPI])
  AC_DEFINE(GMX_LIB_MPI,,[Use a real MPI library])
  PKG_CFLAGS="$PKG_CFLAGS -DGMX_LIB_MPI"

  if test "$enable_mpi_in_place" = "yes"; then
    # not all MPIs support MPI_IN_PLACE:
    AC_MSG_CHECKING([for MPI_IN_PLACE in collective operations])
    AC_TRY_COMPILE([
#include <mpi.h>
void test(void)
{
    void *buf;
    MPI_Allreduce(MPI_IN_PLACE, buf, 10, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}
    ],,[
        AC_MSG_RESULT([yes])
        AC_DEFINE(MPI_IN_PLACE_EXISTS,[],[MPI_IN_PLACE available for collective operations])
    ],[
        AC_MSG_RESULT([no])
    ]) 
  fi

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
  PKG_CFLAGS="$PKG_CFLAGS -DGMX_CRAY_XT3"
  AC_DEFINE(GMX_NO_SYSTEM,[],[Ignore calls to system(3)])
  AC_DEFINE(GMX_NO_NICE,[],[Ignore calls to nice(3)])
  AC_DEFINE(GMX_BROKEN_CALLOC,[],[Don't use calloc(3)])
  extra_LDFLAGS="$extra_LDFLAGS -lgmalloc -lacml"
],[
  AC_MSG_RESULT([no])
])
# end of "$enable_mpi" = "yes"
fi


AH_TEMPLATE([F77_OR_C_FUNC],
            [Set to F77_FUNC(name,NAME) if Fortran used, otherwise 'name' for C.])
AH_TEMPLATE([F77_OR_C_FUNC_],
            [Set to F77_FUNC_(name,NAME) if Fortran used, otherwise 'name' for C.])



if test "$enable_fortran" = "yes" -o "$enable_power6" = "yes"; then
  AC_F77_LIBRARY_LDFLAGS
  AC_F77_WRAPPERS
  AC_DEFINE([F77_OR_C_FUNC(name,NAME)],[F77_FUNC(name,NAME)])
  AC_DEFINE([F77_OR_C_FUNC_(name,NAME)],[F77_FUNC_(name,NAME)])
else
  AC_DEFINE([F77_OR_C_FUNC(name,NAME)],[name])
  AC_DEFINE([F77_OR_C_FUNC_(name,NAME)],[name])
fi



AC_PROG_CPP
# Needed by libtool scripts?
AC_PROG_EGREP
AC_PROG_LN_S

if test "$enable_cpu_optimization" = "yes"; then
  ACX_DETECT_GMXCPU
else
  gmxcpu=""
fi
ACX_COMPILER_MAXOPT


# We've postponed the threads check because we need to check with the right
# compilers/flags
SAVED_LIBS="$LIBS"
if test "$enable_threads" = "yes"; then 
  if test "$enable_mpi" = "yes"; then
    AC_MSG_ERROR([Cannot use both threads and MPI - yet.])
  fi
  if test "$with_fft" = "fftw2"; then
    AC_MSG_ERROR([fftw2 can't be used with threads. Use fftw3 or mkl.])
  fi 
  if test "$enable_fortran" = "yes"; then
    AC_MSG_ERROR([FORTRAN is incompatible with threads and only provides a speed-up on certain IBM compilers. Use --disable-threads if you really want to use FORTRAN kernels.])
  fi
  AC_CHECK_HEADERS(unistd.h)
  AC_CHECK_HEADERS(sys/time.h)
  AC_CHECK_HEADERS(sched.h)
  AC_CHECK_FUNCS(sysconf)

  # for now we just assume pthreads
  ACX_PTHREAD(,AC_MSG_ERROR([Cannot find POSIX threads - required for thread support using configure!]))
  LIBS="$PTHREAD_LIBS $LIBS"
  CFLAGS="$CFLAGS $PTHREAD_CFLAGS -I$srcdir/include"
  CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS -I$srcdir/include"
  CC="$PTHREAD_CC "
  AC_DEFINE(THREAD_PTHREADS,,[Use pthreads for multithreading])

  # profiling
  AC_ARG_ENABLE(tmpi-profiling,
              [AC_HELP_STRING([--enable-tmpi-profiling],
              [Enables profiling of thread_mpi calls and wait times])],,
              enable_tmpi_profiling=no)
  if test "x$enable_tmpi_profiling" = "xyes"; then
    AC_DEFINE(TMPI_PROFILE,,[Enable profiling for call and wait times])
  fi

  # check the atomics
  AC_MSG_CHECKING(if atomic operations are supported)
  CFLAGS_RET="$CFLAGS"
  CXXFLAGS_RET="$CXXFLAGS"
  CFLAGS="$CFLAGS -I$srcdir/include -DTMPI_CHECK_ATOMICS"
  CXXFLAGS="$CXXFLAGS -I$srcdir/include -DTMPI_CHECK_ATOMICS"
  # this sets a variable 'tmpi_atomics' to 'yes' if there are atomics. 
  AC_COMPILE_IFELSE([#include "thread_mpi/atomic.h" ],
                      [AC_MSG_RESULT(yes); tmpi_atomics=yes],
                      [AC_MSG_WARN(
[No atomics.

This compiler+CPU combination doesn't have working atomics.
Thread support will be unbearably slow: disable threads.

Atomics should work on all but the most obscure CPU+compiler combinations;
if your system is not obscure -- like, for example, x86 with gcc --  please
contact the developers.
]); tmpi_atomics=no])


  # pthread_setaffinity
  AC_MSG_CHECKING(if pthread_setaffinity_np is available)
  AC_COMPILE_IFELSE([#define _GNU_SOURCE
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
int main(void)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(0, &set);
    pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
    return 0;
}],
  [AC_MSG_RESULT(yes); pthread_setaffinity=yes],
  [AC_MSG_RESULT(no); pthread_setaffinity=no] )

  if test "x$pthread_setaffinity" = "xyes"; then
    AC_DEFINE(HAVE_PTHREAD_SETAFFINITY,1,[pthread_setaffinity_np is available])
  fi

  CFLAGS="$CFLAGS_RET"
  CXXFLAGS="$CXXFLAGS_RET"

  AC_DEFINE(GMX_THREADS,,[Use threads for parallelization])
  PKG_CFLAGS="$PKG_CFLAGS -DGMX_THREADS"
  AC_DEFINE(GMX_MPI,,[Make a parallel version of GROMACS using MPI])
  AM_CONDITIONAL(THREAD_PARALLEL,true)
else
  AM_CONDITIONAL(THREAD_PARALLEL,false)
fi
LIBS="$SAVED_LIBS"



#############
# Check integer endian
#############
AC_C_BIGENDIAN([AC_DEFINE([GMX_INTEGER_BIG_ENDIAN],,[Integer byte order is big endian.])],[],[
AC_MSG_ERROR([Cannot determine endian in compiled output])])
# Little endian implied if GMX_INTEGER_BIG_ENDIAN not defined in source.

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
# Check that SIGUSR1 is defined
#############
AC_MSG_CHECKING([for SIGUSR1])
AC_TRY_COMPILE([#include <signal.h>],[int i=SIGUSR1;],[
AC_MSG_RESULT([yes])
AC_DEFINE(HAVE_SIGUSR1,,[SIGUSR1 supported and defined])
],AC_MSG_RESULT([no]))

#############
# Check if pipes are supported
#############
AC_MSG_CHECKING([for pipes])
AC_TRY_COMPILE([#include <stdio.h>],[
FILE *pipe=popen("/tmp/pipe","r");
],[
AC_MSG_RESULT([yes])
AC_DEFINE(HAVE_PIPES,,[Pipes (popen/pclose) supported])
],AC_MSG_RESULT([no]))



#############
# Check floating-point format and endian
#############
ACX_FLOAT_FORMAT


AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_LIBTOOL_WIN32_DLL
AC_LIBTOOL_DLOPEN
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

#save current LIBS
SAVED_LIBS="$LIBS"

#set empty defaults
PKG_FFT=""
FFT_LIBS=""
case "$with_fft" in
fftw2)
  if test "$enable_float" = "yes"; then
    sizeof_real=4
  else
    sizeof_real=8
  fi
  ACX_CHECK_FFTW2(fftw,$sizeof_real)
  ACX_CHECK_FFTW2(rfftw,$sizeof_real)
  FFT_LIBS="-lfftw -lrfftw"

  case ${ac_fftw_savedprefix} in
   d) AC_DEFINE(FFTW2_NAME_DFFTW,,[Use the d prefix on fftw2 includes]) ;;
   s) AC_DEFINE(FFTW2_NAME_SFFTW,,[Use the s prefix on fftw2 includes]) ;;
   *) AC_DEFINE(FFTW2_NAME_FFTW,,[Dont use any prefix on fftw2 includes]) ;;
  esac

  AC_DEFINE(GMX_FFT_FFTW2,,[Use FFTW2 FFT library])

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
  PKG_FFT="fftw3f"
  FFT_LIBS="-lfftw3f"
else
  AC_CHECK_LIB([fftw3],main,,AC_MSG_ERROR([Cannot find fftw3 library]))
  PKG_FFT="fftw3"
  FFT_LIBS="-lfftw3"
fi
  AC_DEFINE(GMX_FFT_FFTW3,,[Use FFTW3 FFT library])
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
#  AC_CHECK_LIB([mkl],DftiComputeForward,,
#    AC_MSG_ERROR([Cannot find Intel Math Kernel Library >= 6.0]))
  # first check for older MKLs
  AC_CHECK_LIB([mkl],DftiComputeForward,,[have_mkl="no"])
  if test "$have_mkl" = "no"; then
    # we assume it's a new MKL
    LIBS="$LIBS -liomp5 -lmkl_sequential -lmkl_core"
    LMKL=" -liomp5 -lmkl_sequential -lmkl_core"
    AC_CHECK_LIB([mkl_intel_lp64],DftiComputeForward,,[have_mkl10_64="no"])
    FFT_LIBS="$LMKL -lmkl_intel_lp64"
    if test "$have_mkl10_64" = "no"; then
      AC_CHECK_LIB([mkl_intel],DftiComputeForward,,[have_mkl10_32="no"])
      FFT_LIBS="$LMKL -lmkl_intel"
      if test "$have_mkl10_32" = "no"; then
        AC_MSG_ERROR([Cannot find Intel Math Kernel Library >= 6.0])
      fi
    fi
  fi
  AC_DEFINE(GMX_FFT_MKL,,[Use Intel MKL FFT library])
  AC_SUBST(LMKL)
  ;;

fftpack | no)
  AC_MSG_WARN([The built-in FFTPACK routines are slow.
Are you sure you don't want to use FFTW? It is free and much faster...])
  AC_DEFINE(GMX_FFT_FFTPACK,,[Use Built-in FFTPACK FFT library])
  ;;
*)
  AC_MSG_ERROR([Unknown FFT library])
  ;;
esac
#restore LIBS from before fft check
LIBS="$SAVED_LIBS"

#only fftw3 has support for pkg-config
AC_SUBST(PKG_FFT)
if test -z "${PKG_FFT}"; then
  AC_SUBST(PKG_FFT_LIBS,"$FFT_LIBS")
else
  AC_SUBST(PKG_FFT_LIBS,"")
fi
#flag for Makefile
AC_SUBST(FFT_LIBS)

if test "$enable_fftw_measure" != "yes"; then
   AC_DEFINE(GMX_DISABLE_FFTW_MEASURE,,[Do not optimize FFTW setups (not needed with SSE FFT kernels)])
fi

######
# check for xtc headers
case "${host_os}" in
cygwin* | mingw*)   # use our own xdr stuff for cygwin/mingw
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
    PKG_CFLAGS="$PKG_CFLAGS -DGMX_INTERNAL_XDR"
    AC_MSG_WARN([Couldn't find XDR headers and/or libraries - using our own])
  fi
  ;;
esac


########################################################################
# Checks for header files and types
########################################################################

AC_CHECK_HEADERS([limits.h unistd.h],,AC_MSG_ERROR([Include headers not found]))
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(direct.h)
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(sys/time.h)
AC_CHECK_HEADERS(dirent.h,,[AC_MSG_WARN(
[
********************************************************************
* WARNING! We could not find the header file dirent.h              *
*                                                                  *
* Without this header, we will not be able to scan directories for *
* files, which completely breaks pdb2gmx. You will still be able   *
* to run simulations, but you cannot prepare input data.           *
* dirent.h is part of the POSIX 2001-1 standard, and should be     *
* present on all UNIX-like systems, CygWin, and MinGW.             *
********************************************************************])])
## NOTE: On native windows we do not use configure.ac
AC_CHECK_HEADERS(regex.h,,[AC_MSG_WARN(
[No regex.h found.
Regular expression matching in selections will be disabled.
Selections still work, but the default groups can differ from
those obtained with regular expressions enabled.
Consider installing the regular expression library from GNU libc.])])



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
AC_CHECK_FUNCS(gettimeofday)
AC_CHECK_FUNCS(cbrt)
AC_CHECK_FUNCS(isnan)
AC_CHECK_FUNCS(_isnan)
AC_CHECK_FUNCS(isfinite)
AC_CHECK_FUNCS(_isfinite)
AC_CHECK_FUNCS(fsync)
AC_CHECK_FUNCS(fileno)

# check for bool (must come late, since CFLAGS/includes might affect it)
AC_CHECK_TYPES([bool])
AC_SEARCH_LIBS(posix_memalign, c, AC_DEFINE(HAVE_POSIX_MEMALIGN, 1, Define to 1 if you have the posix_memalign function))
AC_SEARCH_LIBS(memalign, c, AC_DEFINE(HAVE_MEMALIGN, 1, Define to 1 if you have the memalign function))
AC_CHECK_FUNCS(_aligned_malloc)
#AC_SEARCH_LIBS(_aligned_malloc, c, AC_DEFINE(HAVE__ALIGNED_MALLOC, 1, Define to 1 if you have the _aligned_malloc function))






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
  if test $ac_cv_sizeof_voidp = 8; then
     enable_ia32_sse=no;
  else
     enable_x86_64_sse=no;
  fi
  ;;
  *) ;;
esac




#####
# Look for X
AC_PATH_XTRA
if test "$no_x" != "yes"; then
  AM_CONDITIONAL(USE_X11,true)
else
  AM_CONDITIONAL(USE_X11,false)
fi

# GNU Scientific Library or not
AC_ARG_WITH(gsl,
 [AC_HELP_STRING([--with-gsl],[Link to the GNU scientific library, ]
                 [enables extra functions in analysis])],,with_gsl=no)

######
PKG_GSL=""
GSL_LIBS=""
if test "$with_gsl" = "yes"; then
  SAVED_LIBS="$LIBS"
  AC_CHECK_HEADERS([gsl/gsl_version.h],[
# header found, check for libraries 
  AC_CHECK_LIB(gslcblas,main,[
# gslcblas library found, check main lib
  AC_CHECK_LIB(gsl,main,[
# everything need for gsl found
  GSL_LIBS="-lgsl -lgslcblas"
#gsl has pkg-config support
  PKG_GSL="gsl"
  AC_DEFINE(HAVE_LIBGSL,1,[Do we have libgsl])
  ],[:],[-lgslcblas])])])
  LIBS="$SAVED_LIBS"
#
fi
AC_SUBST(PKG_GSL)
AC_SUBST(GSL_LIBS)

AC_ARG_WITH(xml,
 [AC_HELP_STRING([--without-xml],[do not link to the xml2 library,] 
                 [disallows the use of certain file formats])],,with_xml=yes)
 
#########
# Check for libxml2, but it is optional, so dont stop
# if it isnt there.
PKG_XML=""
XML_LIBS=""
if test "$with_xml" = "yes"; then
   SAVED_LIBS="$LIBS"
   ac_save_CPPFLAGS="$CPPFLAGS"
   CPPFLAGS="$ac_save_CPPFLAGS -I/usr/include/libxml2"
   AC_CHECK_HEADERS([libxml/parser.h],AC_CHECK_LIB(xml2,main))
   if test "$ac_cv_lib_xml2_main" = "yes"; then
     INCLUDES="$INCLUDES -I/usr/include/libxml2"
     PKG_XML="libxml-2.0"
     XML_LIBS="-lxml2"
     if test "$enable_all_static" = "yes"; then
       LDFLAGS_RET=$LDFLAGS
       LDFLAGS="-static $LDFLAGS"
       LIBS_RET="$LIBS"
       LIBS="$LIBS -lz"
       AC_MSG_CHECKING(for broken libxml2)
       AC_TRY_RUN([
#include <libxml/parser.h>
main() {
  xmlInitParser();
  return 0;
}],
       AC_MSG_RESULT([no]),[AC_MSG_RESULT([yes])
       AC_MSG_ERROR([Your libxml2 is broken for static linking, recompile it without threads or build gromacs without libxml])],[-lz])
       LDFLAGS=$LDFLAGS_RET
       LIBS="$LIBS_RET"
     fi
   fi
# restore CPPFLAGS
   CPPFLAGS="$ac_save_CPPFLAGS"
   LIBS="$SAVED_LIBS"
fi
AC_SUBST(PKG_XML)
AC_SUBST(XML_LIBS)

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






#####
# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_INLINE
AC_C_RESTRICT

# AC_TYPE_UID_T is broken on Mac OS X, so use generic tests instead
AC_CHECK_TYPES([uid_t],[],[AC_DEFINE(uid_t, int, [Define to `int' if <sys/types.h> doesn't define.])],[[#include <sys/types.h>]])
AC_CHECK_TYPES([gid_t],[],[AC_DEFINE(gid_t, int, [Define to `int' if <sys/types.h> doesn't define.])],[[#include <sys/types.h>]])

AC_CHECK_TYPES([sig_atomic_t],[],[AC_DEFINE(sig_atomic_t, int, [Define to `int' (and hope it's atomic) if <signal.h> doesn't define it.])],[[#include <signal.h>]])

    


# Test stuff for ia32 and x86_64 assembly, and add objects/libraries.
if test "$enable_ia32_sse" = "yes" -o "$enable_x86_64_sse" = "yes"; then
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
    AC_DEFINE([GMX_PPC_ALTIVEC],,[Use PowerPC Altivec inner loops])],[
      AC_MSG_RESULT([no])
      AC_MSG_WARN([no altivec support found - On recent PowerPC CPUs]
[you can get a huge performance gain by using an altivec-enabled compiler.]
[On OS X, you need GCC>=3.3, or the IBM compilers. If you are using standard]
[GCC, you need at least version 4.0])
      enable_ppc_altivec=no])
  fi  
fi

if test "$enable_ia64_asm" = "yes"; then
  AC_DEFINE(GMX_IA64_ASM,,[Use ia64 assembly tuned for Itanium2])
  AC_DEFINE(DISABLE_WATER_NLIST,,[Turn off all water neighborlist optimization])
  AC_DEFINE(DISABLE_WATERWATER_NLIST,,[Turn off water-water neighborlist optimization only])
fi

if test "$enable_bluegene" = "yes"; then
  AC_DEFINE(GMX_BLUEGENE,,[Use assembly intrinsics kernels for BlueGene])
  AC_DEFINE_UNQUOTED(GMX_POWERPC_INVSQRT,,[Use the PowerPC hardware 1/sqrt(x)])
  PKG_CFLAGS="$PKG_CFLAGS -DGMX_POWERPC_INVSQRT"
fi

if test "$enable_fortran" = "yes"; then
  AC_DEFINE(GMX_FORTRAN,,[Enable Fortran support])
fi

##
# Add the appropriate assembly loops 
#
AM_CONDITIONAL([GMX_IA32_SSE],[test "$enable_ia32_sse" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_IA32_SSE2],[test "$enable_ia32_sse" = "yes" -a "$enable_float" = "no"])
AM_CONDITIONAL([GMX_X86_64_SSE],[test "$enable_x86_64_sse" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_X86_64_SSE2],[test "$enable_x86_64_sse" = "yes" -a "$enable_float" = "no"])
AM_CONDITIONAL([GMX_FORTRAN],[test "$enable_fortran" = "yes"])
AM_CONDITIONAL([GMX_PPC_ALTIVEC],[test "$enable_ppc_altivec" = "yes" -a "$enable_float" = "yes"])
AM_CONDITIONAL([GMX_IA64_ASM],[test "$enable_ia64_asm" = "yes"])
AM_CONDITIONAL([GMX_BLUEGENE],[test "$enable_bluegene" = "yes"])
AM_CONDITIONAL([GMX_POWER6],[test "$enable_power6" = "yes"])


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


############
# Check dlopen
###########
if test "$with_dlopen" = "yes"; then
      AC_MSG_CHECKING(if dlopen supported)
      CFLAGS_RET=$CFLAGS
      #LDFLAGS="$lt_cv_dlopen_libs $LDFLAGS" #can't make the macro, which is getting lt_cv_dlopen_libs, to work
      LDFLAGS_RET=$LDFLAGS
      LIBS_RET=$LIBS
      if test "$enable_all_static" = "yes"; then  #make sure we test also whether it works static
      	 LDFLAGS="$LDFLAGS -static"
      fi
      CFLAGS="-I$srcdir/include -DGMX_DLOPEN $CFLAGS"
      LIBS="$DLOPEN_LIBS $LIBS"
      AC_TRY_LINK([#include "$srcdir/src/gmxlib/vmddlopen.c"],,[
        AC_MSG_RESULT([yes])
	AC_DEFINE(GMX_DLOPEN,,[Compile with dlopen])
	AC_SUBST(DLOPEN_LIBS)
      ],[
	AC_MSG_RESULT([no])
	AC_SUBST(DLOPEN_LIBS,"")
      ])
      CFLAGS=$CFLAGS_RET
      LIBS=$LIBS_RET
      LDFLAGS=$LDFLAGS_RET
fi




########################################################################
# Final output stage
########################################################################
AS="$CC"
AC_SUBST(ASFLAGS)
AC_SUBST(AS)
AC_SUBST(INCLUDES)          # should be automatic, but doesnt seem to be?

# To tell libraries apart, we use four versions of the name. For instance, libgmx:
# libgmx		single prec, no mpi.
# libgmx_d		double prec, no mpi.
# libgmx_mpi		single prec, mpi.
# libgmx_mpi_d		double prec, mpi.
# However, the non-suffixed names are linked to the _mpi and/or _d suffixed
# ones upon installation if that is the only version available.

AC_ARG_VAR(
  [LIBSUFFIX],
  [To tell libraries apart, a suffix can be used. This is normally done
   automatically during ./configure (_d|_mpi|_mpi_d). But with this
   variable you can overwrite it
  ]
)
if test -z "${LIBSUFFIX}" ; then
  if test "$enable_mpi" = "yes"; then
    LIBSUFFIX="_mpi"
  fi
  if test "$enable_float" = "no"; then
    LIBSUFFIX="${LIBSUFFIX}_d"
  fi
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

if test "$enable_fahcore" = "yes"; then  
  CPPFLAGS="$CPPFLAGS -I\$(top_srcdir)/../corewrap -include swindirect.h -DFULLINDIRECT -DUSE_FAH_XDR"  
  if test "$enable_mpi" = "yes"; then 
    CFLAGS="$CFLAGS -DMPI"  
  fi
fi
AM_CONDITIONAL([GMX_FAHCORE],[test "$enable_fahcore" = "yes"])


if test "$enable_all_static" = "yes"; then
  LDFLAGS="$LDFLAGS -all-static"
fi

# we have asm source, so this is necessary for automake 1.6
CCAS=$CC
CCASFLAGS=$CFLAGS
 
AC_SUBST(CCAS)
AC_SUBST(CCASFLAGS)
AC_SUBST(PKG_CFLAGS)

AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ src/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/libgmx${LIBSUFFIX}.pc:src/gmxlib/libgmx.pc.in ])
AC_CONFIG_FILES([ src/gmxlib/gmx_blas/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/gmx_lapack/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/selection/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/statistics/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/trajana/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/thread_mpi/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_c/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia32_sse/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia32_sse2/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_x86_64_sse/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_x86_64_sse2/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ppc_altivec/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia64_single/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_ia64_double/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_bluegene/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_power6/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_f77_single/Makefile ])
AC_CONFIG_FILES([ src/gmxlib/nonbonded/nb_kernel_f77_double/Makefile ])
AC_CONFIG_FILES([ include/Makefile ])
AC_CONFIG_FILES([ include/types/Makefile ])
AC_CONFIG_FILES([ include/thread_mpi/Makefile ])
AC_CONFIG_FILES([ include/thread_mpi/atomic/Makefile ])
AC_CONFIG_FILES([ src/mdlib/Makefile ])
AC_CONFIG_FILES([ src/mdlib/libmd${LIBSUFFIX}.pc:src/mdlib/libmd.pc.in ])
AC_CONFIG_FILES([ src/kernel/Makefile ])
AC_CONFIG_FILES([ src/kernel/libgmxpreprocess${LIBSUFFIX}.pc:src/kernel/libgmxpreprocess.pc.in ])
AC_CONFIG_FILES([ src/tools/Makefile ])
AC_CONFIG_FILES([ src/tools/libgmxana${LIBSUFFIX}.pc:src/tools/libgmxana.pc.in ])
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
AC_CONFIG_FILES([ share/top/gmx.ff/Makefile ])
AC_CONFIG_FILES([ share/top/gmx2.ff/Makefile ])
AC_CONFIG_FILES([ share/top/encadv.ff/Makefile ])
AC_CONFIG_FILES([ share/top/encads.ff/Makefile ])
AC_CONFIG_FILES([ share/top/oplsaa.ff/Makefile ])
AC_CONFIG_FILES([ share/top/gromos43a1.ff/Makefile ])
AC_CONFIG_FILES([ share/top/gromos43a2.ff/Makefile ])
AC_CONFIG_FILES([ share/top/gromos45a3.ff/Makefile ])
AC_CONFIG_FILES([ share/top/gromos53a5.ff/Makefile ])
AC_CONFIG_FILES([ share/top/gromos53a6.ff/Makefile ])
AC_CONFIG_FILES([ share/top/charmm27.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amber94.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amber96.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amberGS.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amber99.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amber99sb.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amber99sb-ildn.ff/Makefile ])
AC_CONFIG_FILES([ share/top/amber03.ff/Makefile ])
AC_CONFIG_FILES([ share/html/Makefile ])
AC_CONFIG_FILES([ share/html/images/Makefile ])
AC_CONFIG_FILES([ share/html/online/Makefile ])
AC_CONFIG_FILES([ man/Makefile man/man1/Makefile man/man7/Makefile ])
AC_CONFIG_FILES([ Doxyfile ])
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
  echo ""
  echo "  WARNING:"
  echo "  There are known problems with some MPI implementations:"
  echo "             OpenMPI version < 1.4.1" 
  echo "             MVAPICH2 version <= 1.4.1"
fi

if test "$enable_shared" = "no"; then
  echo ""
  echo "* On most platforms you can save 10X space with dynamic libraries, although"
  echo "  the binaries might be less portable. Why not try --enable-shared ?"
fi  


if test "$tmpi_atomics" = "no"; then
  echo ""
  echo "* This compiler+CPU combination doesn't have working atomic operations."
  echo "  Thread support will be unbearably slow: the default number of threads"
  echo "  is set to 1."
  echo "  Atomics should work on all but the most obscure CPU+compiler"
  echo "  combinations; if your system is not obscure (like, for example, x86"
  echo "  with gcc) please contact the developers."
  echo ""
fi


