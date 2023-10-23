#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(SZ, 2.1, sdi1@anl.gov)
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects])
AC_CONFIG_SRCDIR([sz/src/sz.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
: ${CFLAGS=-O3 -std=c99 -Wall}
AC_PROG_CC
AC_PROG_FC
AC_PROG_RANLIB
AM_PROG_AR
AC_PROG_LIBTOOL
#m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AM_MAINTAINER_MODE([disable])

AC_DEFUN([AC_F90_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(Fortran 90)dnl
_AC_PROG_F90_V
AC_CACHE_CHECK([for Fortran 90 libraries], ac_cv_flibs,
[if test "x$FLIBS" != "x"; then
  ac_cv_flibs="$FLIBS" # Let the user override the test.
else

_AC_PROG_F90_V_OUTPUT

ac_cv_flibs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_f90_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
              ac_cv_flibs="$ac_cv_flibs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_flibs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_flibs="$ac_cv_flibs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
             "" | -*);;
             *)
		ac_arg="$ac_arg$[2]"
		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_flibs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_flibs="$ac_cv_flibs $ac_j"])
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
                             ac_cv_flibs="$ac_cv_flibs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_f90_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_flibs)
      ;;
esac
fi # test "x$FLIBS" = "x"
])
FLIBS="$ac_cv_flibs"
AC_SUBST(FLIBS)
AC_LANG_POP(Fortran 90)dnl
])# AC_F90_LIBRARY_LDFLAGS

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h sys/time.h unistd.h])

#
AC_ARG_ENABLE([fortran],
              AS_HELP_STRING([--enable-fortran], [Enable fortran compilation]))
AM_CONDITIONAL([FORTRAN], [test "x$enable_fortran" = "xyes"])


##--------------------------------------------------------------------
##
## GSL  (GNU Scientific Library)
##
##  also need to find out a way to replace -lgslcblas with -lcblas
##  or if ATLAS has been installed, with -lcblas -latlas
## --  Also note the -DHAVE_GSL is a bit of a kludge, since i cannot
## --  figure out how to do this properly through config.h
##--------------------------------------------------------------------


AC_ARG_ENABLE(gsl, [  --enable-gsl            use GSL], ok=$enableval, ok=no)
AM_CONDITIONAL([GSL], [test "x$enable_gsl" = "xyes"])

if test "$ok" = "yes"; then
  AM_PATH_GSL(0.9.2, ok=yes, ok=no)
  if test "$ok" = "yes"; then
    GSL_CFLAGS="$GSL_CFLAGS -DHAVE_GSL"
    AC_DEFINE(HAVE_LIBGSL)
  else
    AC_MSG_WARN([could not find GSL library])
  fi
else
  AC_MSG_WARN([GSL disabled])
fi

AC_SUBST(GSL_CFLAGS)

# Check GSL library

if test "$ok" = "yes"; then
AC_ARG_WITH([gsl-prefix], [AS_HELP_STRING([--with-gsl-prefix=DIR], 
            [installation prefix for gsl])], 
            [with_gsl_include_path=$with_gsl_prefix/include; 
             with_gsl_lib_path=$with_gsl_prefix/lib])

AC_ARG_WITH([gsl-include-path], [AS_HELP_STRING([--with-gsl-include-path=DIR],
            [path to the gsl directory, defaults to /usr/include])], [])

if test "x$with_gsl_include_path" = x ; then
    with_gsl_include_path="/usr/include"
fi

AC_CHECK_FILE([${with_gsl_include_path}/gsl/gsl_cblas.h], [],
              AC_MSG_ERROR([Invalid path to gsl/gsl_cblas.h]))
GSL_HDR="-I$with_gsl_include_path"

AC_ARG_WITH([gsl-lib-path], [AS_HELP_STRING([--with-gsl-lib-path],
            [path to the  gsl libraries, defaults to /usr/lib])], [])
if test "x$with_gsl_lib_path" = x ; then
    with_gsl_lib_path="/usr/lib"
fi
GSL_LIB="-L$with_gsl_lib_path"
LIBS="$GSL_LIB"
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [],
              AC_MSG_ERROR([Function cblas_dgemm in gslcblas library not found]))
AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [],
              AC_MSG_ERROR([Function gsl_blas_dgemm in gsl library not found]))

if test "x$BUILD_STATIC" = xyes ; then
    GSL_STATIC="${with_gsl_lib_path}/libgsl.a ${with_gsl_lib_path}/libgslcblas.a"
    GSL_LIB=
else
    GSL_STATIC="-lgsl -lgslcblas"
fi

AC_SUBST(GSL_HDR)
AC_SUBST(GSL_LIB)
AC_SUBST(GSL_STATIC)
fi

## 
## OPENMP
##

AC_ARG_ENABLE(openmp, [  --enable-openmp            use OPENMP], ok=$enableval, ok=no)
AM_CONDITIONAL([OPENMP], [test "x$enable_openmp" = "xyes"])

if test "$ok" = "yes"; then
        OPENMP_FLAGS="-DHAVE_OPENMP"
fi
AC_SUBST(OPENMP_FLAGS)

##
## TIME-based compression
##

AC_ARG_ENABLE(timecmpr, [  --enable-timecmpr            use TIMECMPR], ok=$enableval, ok=no)
AM_CONDITIONAL([TIMECMPR], [test "x$enable_timecmpr" = "xyes"])

if test "$ok" = "yes"; then
        TIMECMPR_FLAGS="-DHAVE_TIMECMPR"
fi

AC_SUBST(TIMECMPR_FLAGS)


##
## PASTRI
##

AC_ARG_ENABLE(pastri, [  --enable-pastri            use PASTRI], ok=$enableval, ok=no)
AM_CONDITIONAL([PASTRI], [test "x$enable_pastri" = "xyes"])

if test "$ok" = "yes"; then
        PASTRI_FLAGS="-DHAVE_PASTRI"
fi

AC_SUBST(PASTRI_FLAGS)

##
## RANDOMACCESS
##

AC_ARG_ENABLE(randomaccess, [  --enable-randomaccess           use RANDOMACCESS], ok=$enableval, ok=no)
AM_CONDITIONAL([RANDOMACCESS], [test "x$enable_randomaccess" = "xyes"])

if test "$ok" = "yes"; then
        RANDOMACCESS_FLAGS="-DHAVE_RANDOMACCESS"
fi

AC_SUBST(RANDOMACCESS_FLAGS)

##
## WRITESTATS_FLAGS
##

AC_ARG_ENABLE(writestats, [  --enable-writestats            use WRITESTATS], ok=$enableval, ok=no)
AM_CONDITIONAL([WRITESTATS], [test "x$enable_writestats" = "xyes"])

if test "$ok" = "yes"; then
        WRITESTATS_FLAGS="-DHAVE_WRITESTATS"
fi

AC_SUBST(WRITESTATS_FLAGS)



# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memchr memset strchr strtol])

AC_CONFIG_FILES([Makefile
                 sz/Makefile
                 zlib/Makefile
		 zstd/Makefile
		 example/Makefile])
AC_OUTPUT
