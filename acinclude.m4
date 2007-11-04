# ACX_CHECK_FFTW2()
# ----------------
# This macro checks for fftw-2.x header files and libraries,
# including the possible prefixing with s or d to determine precision.
# Arg 1 is the fftw header/library name to check for, without
# prefix or anything else (e.g. rfftw_mpi for real MPI transforms)
# Arg 2 is the size of the real variable used.
AC_DEFUN([ACX_CHECK_FFTW2],
[
if test -z "$ac_fftw_firstname"; then

sizeof_real=$2
if test $sizeof_real = 8; then
  prec="double"
  fftwcheckprefix=d
else
  prec="single"
  fftwcheckprefix=s
fi

xfftwname=${fftwcheckprefix}$1

ok="no"
# check header doesn't work, since we must use mpicc to get includes,
# we cant trust cpp.
AC_MSG_CHECKING([for $xfftwname.h])
AC_TRY_COMPILE([#include <$xfftwname.h>],,
[
fftwname=$xfftwname
AC_MSG_RESULT(yes)
],
AC_MSG_RESULT(no))

# fftwname was set if we found a header

if test -n "$fftwname"; then
# we cannot run the code since an MPI program might not be allowed
# on a login node of a supercomputer
AC_TRY_COMPILE([#include <$fftwname.h>],
[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $sizeof_real)]; ],
[
ok=yes
usedprefix=$fftwcheckprefix
],[ok=no])
fi

if test "$ok" != "yes"; then
  AC_MSG_CHECKING([for $1.h])
  AC_TRY_COMPILE([#include <$1.h>],,AC_MSG_RESULT(yes),
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Cannot find any $prec precision $xfftwname.h or $1.h]
[Do you have $prec precision FFTW-2.x installed? If you are using packages,]
[note that you also need fftw-devel to compile GROMACS. You can find the ]
[software at www.fftw.org, and detailed instructions at www.gromacs.org.]
[If you compiled FFTW-2.x yourself:                                    ]
[Note that the default FFTW-2.x setup is double precision. Change the FFTW]
[configuration to single with --enable-float. If you want MPI support,]
[use --enable-mpi. It is a good idea to install both single & double.]
[If your sysadm doesn't want to install it you can do it to a location]
[in your home directory and provide the correct paths in the CPPFLAGS]
[and LDFLAGS environment variables before running configure.]
[That is also necessary to do if your compiler doesn't search]
[/usr/local/include and /usr/local/lib by default.]
[You can find information at www.gromacs.org, or in the INSTALL file.])
])
AC_TRY_COMPILE([#include <$1.h>],
[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $sizeof_real)];],
[
usedprefix=""
fftwname=$1
],
[
AC_MSG_ERROR([Cannot find any $prec precision $xfftwname.h or $1.h]
[Do you have $prec precision FFTW-2.x installed? If you are using packages,]
[note that you also need fftw-devel to compile GROMACS. You can find the ]
[software at www.fftw.org, and detailed instructions at www.gromacs.org.]
[If you compiled FFTW-2.x yourself:                                   ]
[Note that the default FFTW-2.x setup is double precision. Change the FFTW]
[configuration to single with --enable-float. If you want MPI support,]
[use --enable-mpi. It is a good idea to install both single & double.]
[If your sysadm doesn't want to install it you can do it to a location]
[in your home directory and provide the correct paths in the CPPFLAGS]
[and LDFLAGS environment variables before running configure.]
[That is also necessary to do if your compiler doesn't search]
[/usr/local/include and /usr/local/lib by default.]
[You can find information at www.gromacs.org, or in the INSTALL file.])])
fi

AC_CHECK_LIB($fftwname,main,,
AC_MSG_ERROR([Can't find a library to match the $fftwname header]))
ac_fftw_savedprefix=$usedprefix
ac_fftw_firstname=$fftwname

else

fftwname=${ac_fftw_savedprefix}$1
AC_MSG_CHECKING([for $fftwname.h])
AC_TRY_COMPILE(
[#include <$fftwname.h>],,
[AC_MSG_RESULT(yes)
LIBS="-l$fftwname $LIBS"
AC_TRY_LINK_FUNC([main],,,
AC_MSG_ERROR([Can't find a library to match the $fftwname header]))],
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Cant find $fftwname.h header. Make sure all your
fftw prefixes match - we already use $ac_fftw_firstname.h])
])

fi

])






dnl Check for floating-point format and double precision word order.
dnl We dont require IEEE, but there are optimizations we can only do with it.
dnl Just as for integers, the bytes in a word can be small of big endian.
dnl There is already a standard autoconf macro (AC_C_BIGENDIAN) that you 
dnl should use to check this for integers - I have never heard of a machine
dnl where it is not the same for integer and fp variables, but we still check
dnl it separately for fp variables here to be sure.
dnl
dnl However, in double precision there are also two ways to arrange the words
dnl forming a double (8-byte=2-word) variable.
dnl Normally this order is the same as the endian, but there are 
dnl exceptions (e.g. ARM)
dnl We detect it by compiling a small test program and grepping into it.
dnl
AC_DEFUN([ACX_FLOAT_FORMAT],
[AC_CACHE_CHECK(floating-point format, acx_float_format,
[cat >conftest.$ac_ext <<EOF
[/* Check that a double is 8 bytes - die if it isnt */
extern char xyz [sizeof(double) == 8 ? 1 : -1];
double abc [] = {
  /* "GROMACSX" in ascii    */
  (double)  3.80279098314984902657e+35 , 
  /* "GROMACSX" in ebcdic   */
  (double) -1.37384666579378297437e+38 , 
  /* "D__float" (vax)       */
  (double)  3.53802595280598432000e+18 , 
  /* "IBMHEXFP" s390/ascii  */
  (double)  1.77977764695171661377e+10 , 
  /* "IBMHEXFP" s390/ebcdic */
  (double) -5.22995989424860458374e+10 };
]
EOF
if AC_TRY_EVAL(ac_compile); then
# dont match first and last letter because of rounding errors.
# next: big-endian - string is GROMACSX 
  if   grep 'ROMACS' conftest.o >/dev/null 2>&1; then
    acx_float_format='IEEE754 (big-endian byte and word order)'
# next: big-endian byte order, but little-endian word order - ACSXGROM
  elif grep 'CSXGRO' conftest.o >/dev/null 2>&1; then
    acx_float_format='IEEE754 (big-endian byte, little-endian word order)'
# next: little-endian - XSCAMORG
  elif grep 'SCAMOR' conftest.o >/dev/null 2>&1; then
    acx_float_format='IEEE754 (little-endian byte and word order)'
# next: little-endian byte order, but big-endian word order - MORGXSCA
  elif grep 'ORGXSC' conftest.o >/dev/null 2>&1; then
    acx_float_format='IEEE754 (big-endian byte, little-endian word order)'
  elif grep '__floa' conftest.o >/dev/null 2>&1; then
    acx_float_format='VAX D-float'
  elif grep 'BMHEXF' conftest.o >/dev/null 2>&1; then
    acx_float_format='IBM 370 hex'
  else
    AC_MSG_WARN([Unknown floating-point format])
  fi
else
  AC_MSG_ERROR(compile failed)
fi
rm -rf conftest*])
case $acx_float_format in
    'IEEE754 (big-endian byte and word order)' )
       format=IEEE754
       byteorder=big
       wordorder=big            
       ;;
    'IEEE754 (little-endian byte and word order)' )
       format=IEEE754
       byteorder=little
       wordorder=little
       ;;
    'IEEE754 (big-endian byte, little-endian word order)' )
       format=IEEE754
       byteorder=big
       wordorder=little
       ;;
    'IEEE754 (litte-endian byte, big-endian word order)' )
       format=IEEE754
       byteorder=little
       wordorder=big            
       ;;
    'VAX D-float' )
       AC_DEFINE(FLOAT_FORMAT_VAX,,[VAX floating-point format if set])
       ;;
    'IBM 370 hex' )
       AC_DEFINE(FLOAT_FORMAT_IBM_HEX,,[IBM HEX floating-point format if set (s390?)])
       ;;   
     * )
       format=Unknown   
       ;;
esac
if test "$format" = "IEEE754"; then
       AC_DEFINE(FLOAT_FORMAT_IEEE754,,[IEEE754 floating-point format. Memory layout is defined by
macros IEEE754_BIG_ENDIAN_BYTE_ORDER and IEEE754_BIG_ENDIAN_WORD_ORDER.])
fi
if test "$byteorder" = "big"; then
  AC_DEFINE(IEEE754_BIG_ENDIAN_BYTE_ORDER,,[Bytes in IEEE fp word are in big-endian order if set,
 little-endian if not. Only relevant when FLOAT_FORMAT_IEEE754 is defined.])
fi
if test "$wordorder" = "big"; then
  AC_DEFINE(IEEE754_BIG_ENDIAN_WORD_ORDER,,[The two words in a double precision variable are in b
ig-endian order if set, little-endian if not. Do NOT assume this is the same as the byte order! 
Only relevant when FLOAT_FORMAT_IEEE754 is defined.])
fi
])






dnl
dnl
dnl AC_FIND_MOTIF : find OSF/Motif or LessTif, and provide variables
dnl     to easily use them in a Makefile.
dnl
dnl Adapted from a macro by Andreas Zeller.
dnl
dnl The variables provided are :
dnl     link_motif              (e.g. -L/usr/lesstif/lib -lXm -lXt)
dnl     include_motif           (e.g. -I/usr/lesstif/lib)
dnl     motif_libraries         (e.g. /usr/lesstif/lib)
dnl     motif_includes          (e.g. /usr/lesstif/include)
dnl
dnl The link_motif and include_motif variables should be fit to put on
dnl your application's link line in your Makefile.
dnl
AC_DEFUN([AC_FIND_MOTIF],
[
AC_REQUIRE([AC_PATH_XTRA])

motif_includes=
motif_libraries=

dnl AC_ARG_WITH(motif,
dnl [  --without-motif         do not use Motif widgets])
dnl Treat --without-motif like
dnl --without-motif-includes --without-motif-libraries.

if test "$no_x" = "yes"
then
  motif_includes=no
  motif_libraries=no
fi

AC_ARG_WITH(motif-includes,
[  --with-motif-includes=DIR     Motif include files are in DIR],
motif_includes="$withval")

AC_ARG_WITH(motif-libraries,
[  --with-motif-libraries=DIR    Motif libraries are in DIR],
motif_libraries="$withval")


AC_MSG_CHECKING(for Motif)

#
#
# Search the include files.
#
if test "$motif_includes" = ""; then
AC_CACHE_VAL(ac_cv_motif_includes,
[
ac_motif_save_LIBS="$LIBS"
ac_motif_save_INCLUDES="$INCLUDES"
ac_motif_save_CPPFLAGS="$CPPFLAGS"
ac_motif_save_LDFLAGS="$LDFLAGS"
#
LIBS="$X_PRE_LIBS -lXm -lXt -lX11 $X_EXTRA_LIBS $LIBS"
INCLUDES="$X_CFLAGS $INCLUDES"
CPPFLAGS="$X_CFLAGS $CPPFLAGS"
LDFLAGS="$X_LIBS $LDFLAGS"
#
ac_cv_motif_includes="no"
AC_TRY_COMPILE([#include <Xm/Xm.h>],[int a;],
[
# Xm/Xm.h is in the standard search path.
ac_cv_motif_includes=
],
[
# Xm/Xm.h is not in the standard search path.
# Locate it and put its directory in `motif_includes'
#
# /usr/include/Motif* are used on HP-UX (Motif).
# /usr/include/X11* are used on HP-UX (X and Athena).
# /usr/dt is used on Solaris (Motif).
# /usr/openwin is used on Solaris (X and Athena).
# Other directories are just guesses.
for dir in "$x_includes" "${prefix}/include" /usr/include /usr/local/include \
           /usr/include/Motif2.0 /usr/include/Motif1.2 /usr/include/Motif1.1 \
           /usr/include/X11R6 /usr/include/X11R5 /usr/include/X11R4 \
           /usr/dt/include /usr/openwin/include \
           /usr/dt/*/include /opt/*/include /usr/include/Motif* \
           "${prefix}"/*/include /usr/*/include /usr/local/*/include \
           "${prefix}"/include/* /usr/include/* /usr/local/include/*; do
if test -f "$dir/Xm/Xm.h"; then
ac_cv_motif_includes="$dir"
break
fi
done
])
#
LIBS="$ac_motif_save_LIBS"
INCLUDES="$ac_motif_save_INCLUDES"
CPPFLAGS="$ac_motif_save_CPPFLAGS"
LDFLAGS="$ac_motif_save_LDFLAGS"
])
motif_includes="$ac_cv_motif_includes"
fi
#
#
# Now for the libraries.
#
if test "$motif_libraries" = ""; then
AC_CACHE_VAL(ac_cv_motif_libraries,
[
ac_motif_save_LIBS="$LIBS"
ac_motif_save_INCLUDES="$INCLUDES"
ac_motif_save_CPPFLAGS="$CPPFLAGS"
ac_motif_save_LDFLAGS="$LDFLAGS"
#
LIBS="$X_PRE_LIBS -lXm -lXt -lX11 $X_EXTRA_LIBS $LIBS"
INCLUDES="$X_CFLAGS $INCLUDES"
CPPFLAGS="$X_CFLAGS $CPPFLAGS"
LDFLAGS="$X_LIBS $LDFLAGS"
#
ac_cv_motif_libraries="no"
AC_TRY_LINK([#include <Xm/Xm.h>],[XtToolkitInitialize();],
[
# libXm.a is in the standard search path.
ac_cv_motif_libraries=
],
[
# libXm.a is not in the standard search path.
# Locate it and put its directory in `motif_libraries'
#
# /usr/lib/Motif* are used on HP-UX (Motif).
# /usr/lib/X11* are used on HP-UX (X and Athena).
# /usr/dt is used on Solaris (Motif).
# /usr/lesstif is used on Linux (Lesstif).
# /usr/openwin is used on Solaris (X and Athena).
# Other directories are just guesses.
for dir in "$x_libraries" "${prefix}/lib" /usr/lib /usr/local/lib \
           /usr/lib/Motif2.0 /usr/lib/Motif1.2 /usr/lib/Motif1.1 \
           /usr/lib/X11R6 /usr/lib/X11R5 /usr/lib/X11R4 /usr/lib/X11 \
           /usr/dt/lib /usr/openwin/lib \
           /usr/dt/*/lib /opt/*/lib /usr/lib/Motif* \
           /usr/lesstif*/lib /usr/lib/Lesstif* \
           "${prefix}"/*/lib /usr/*/lib /usr/local/*/lib \
           "${prefix}"/lib/* /usr/lib/* /usr/local/lib/*; do
if test -d "$dir" && test "`ls $dir/libXm.* 2> /dev/null`" != ""; then
ac_cv_motif_libraries="$dir"
break
fi
done
])
#
LIBS="$ac_motif_save_LIBS"
INCLUDES="$ac_motif_save_INCLUDES"
CPPFLAGS="$ac_motif_save_CPPFLAGS"
LDFLAGS="$ac_motif_save_LDFLAGS"
])
#
motif_libraries="$ac_cv_motif_libraries"
fi
#
# Provide an easier way to link
#
if test "$motif_includes" = "no" -o "$motif_libraries" = "no"; then
        with_motif="no"
else
        with_motif="yes"
fi

if test "$with_motif" != "no"; then
        if test "$motif_libraries" = ""; then
                link_motif="-lXm -lXt"
                MOTIF_LIBS="-lXm -lXt"
        else
                link_motif="-L$motif_libraries -lXm -lXt"
                MOTIF_LIBS="-L$motif_libraries -lXm -lXt"
        fi
        if test "$motif_includes" != ""; then
                include_motif="-I$motif_includes"
                MOTIF_INCLUDES="-I$motif_includes"
        fi
	LIBS="$LIBS $MOTIF_LIBS"
	INCLUDES="$INCLUDES $MOTIF_INCLUDES"
        AC_DEFINE(HAVE_MOTIF,,[Use motif/lesstif libraries])
else
        with_motif="no"
fi
#
#
#
#
motif_libraries_result="$motif_libraries"
motif_includes_result="$motif_includes"
test "$motif_libraries_result" = "" && motif_libraries_result="in default path"
test "$motif_includes_result" = "" && motif_includes_result="in default path"
test "$motif_libraries_result" = "no" && motif_libraries_result="(none)"
test "$motif_includes_result" = "no" && motif_includes_result="(none)"
AC_MSG_RESULT([libraries $motif_libraries_result, headers $motif_includes_result])
	
# seems as if Xm depends on -lXext and/or -lXp on old redhat and OS X. 
ac_motif_save_LIBS="$LIBS"
ac_motif_save_INCLUDES="$INCLUDES"
ac_motif_save_CPPFLAGS="$CPPFLAGS"
ac_motif_save_LDFLAGS="$LDFLAGS"
CPPFLAGS="$CPPFLAGS $X_CFLAGS"
INCLUDE="$INCLUDE $X_CFLAGS"
LDFLAGS="$X_LIBS $LDFLAGS"
# first try both - they are crossdependent! urk...
LIBS="$X_PRE_LIBS -lX11 $X_EXTRA_LIBS $ac_motif_save_LIBS -lXext -lXp"
AC_MSG_CHECKING(for libXext and libXp)
AC_TRY_LINK([#include <Xm/Xm.h>],[XtToolkitInitialize();],
  [AC_MSG_RESULT(yes)
   X_PRE_LIBS="$X_PRE_LIBS -lXext -lXp"],[
   AC_MSG_RESULT(no)
   # both libs didnt work, try libXext separately
   LIBS="$X_PRE_LIBS -lX11 $X_EXTRA_LIBS $ac_motif_save_LIBS -lXext"
   AC_MSG_CHECKING(for only libXext)
   AC_TRY_LINK([#include <Xm/Xm.h>],[XtToolkitInitialize();],
  [AC_MSG_RESULT(yes)
  X_PRE_LIBS="$X_PRE_LIBS -lXext"],[AC_MSG_RESULT(no)])
  ])
LIBS=$ac_motif_save_LIBS
INCLUDES="$ac_motif_save_INCLUDES"
CPPFLAGS=$ac_motif_save_CPPFLAGS
LDFLAGS="$ac_motif_save_LDFLAGS"
])dnl


dnl macro modified from the fftw distribution (www.fftw.org)
AC_DEFUN([ACX_CHECK_CC_FLAGS],
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether $CC accepts $1, ac_$2,
[echo 'void f(){}' > conftest.c
res=`$CC $1 -c conftest.c 2>&1`
#
# The stupid intel compiler echos the filename on stderr...
# 
if test -z "$res" -o "$res" = "conftest.c:"; then
	ac_$2=yes
else
	ac_$2=no
fi
rm -rf conftest*
])
if test "$ac_$2" = yes; then
	:
	$3
else
	:
	$4
fi
])

dnl macro modified from the fftw distribution (www.fftw.org)
AC_DEFUN([ACX_CHECK_F77_FLAGS],
[
AC_REQUIRE([AC_PROG_F77])
AC_CACHE_CHECK(whether $F77 accepts $1, ac_$2,
[cat > conftest.f << EOF
      subroutine f
      return 
      end
EOF
if test -z "`$F77 $1 -c conftest.f `"; then
	ac_$2=yes
else
	ac_$2=no
fi
rm -rf conftest*
])
if test "$ac_$2" = yes; then
	:
	$3
else
	:
	$4
fi
])


# ACX_DETECT_GMXCPU
# ---------------------------
# Macro to extend the exact CPU for some hosts
AC_DEFUN([ACX_DETECT_GMXCPU],
[
AC_REQUIRE([AC_CANONICAL_HOST])

#
# Determine the exact cpu type on some common systems where it is 
# not visible from the host triplet.
# (on e.g. intel and dec/tru64 the host type is enough)

gmxcpu="";

case "${host_cpu}-${host_os}" in

*-aix*)
  # some versions of config.status says these systems are PowerPC even
  # when they have Power3 CPUs (they used to be recognized as rs6000), 
  # so we need to work around that.
  # 
  # we need to fool the combination of m4, sh and awk - thus the seemingly unnecessary n
  if test -f /usr/sbin/lsdev && test -f /usr/sbin/lsattr; then
    IBM_CPU_ID=`/usr/sbin/lsdev -C -c processor -S available | head -1 | awk '{ n=1; print $n }'`
    if /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER5 >/dev/null 2>&1; then
      gmxcpu=power5
    elif /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER4 >/dev/null 2>&1; then
      gmxcpu=power4
    elif /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER3 >/dev/null 2>&1; then
      gmxcpu=power3
    elif /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER2 >/dev/null 2>&1; then
      gmxcpu=power2
    fi
  fi
  if test -z "${gmxcpu}" && test -f /usr/sbin/lscfg; then
    if /usr/sbin/lscfg -vp | grep PowerPC | grep 604 >/dev/null 2>&1; then
      gmxcpu=ppc604
    elif /usr/sbin/lscfg -vp | grep PowerPC | grep 603 >/dev/null 2>&1; then
      gmxcpu=ppc603
    elif /usr/sbin/lscfg -vp | grep PowerPC | grep rs64a >/dev/null 2>&1; then
      gmxcpu=rs64a
    elif /usr/sbin/lscfg -vp | grep PowerPC | grep rs64b >/dev/null 2>&1; then
      gmxcpu=rs64b
    elif /usr/sbin/lscfg -vp | grep PowerPC | grep rs64c >/dev/null 2>&1; then
      gmxcpu=rs64c
    elif /usr/sbin/lscfg -vp | grep POWER2 >/dev/null 2>&1; then
      gmxcpu=power2
    elif /usr/sbin/lscfg -vp | grep POWER3 >/dev/null 2>&1; then
      gmxcpu=power3
    elif /usr/sbin/lscfg -vp | grep POWER4 >/dev/null 2>&1; then
      gmxcpu=power4
    fi
  fi
  ;;

mips*-irix*)
  if /sbin/hinv | grep CPU | grep R12000 >/dev/null 2>&1; then
    gmxcpu=r12000
  elif /sbin/hinv | grep CPU | grep R10000 >/dev/null 2>&1; then
    gmxcpu=r10000
  elif /sbin/hinv | grep CPU | grep R8000 >/dev/null 2>&1; then
    gmxcpu=r8000
  elif /sbin/hinv | grep CPU | grep R5000 >/dev/null 2>&1; then
    gmxcpu=r5000
  else
    gmxcpu=""
  fi
  ;;

sparc*-solaris*)
  if /usr/sbin/prtconf | grep UltraSPARC-III >/dev/null 2>&1; then
    gmxcpu=ultrasparc3
  elif /usr/sbin/prtconf | grep UltraSPARC-IIi >/dev/null 2>&1; then
    gmxcpu=ultrasparc2i
  elif /usr/sbin/prtconf | grep UltraSPARC-II >/dev/null 2>&1; then
    gmxcpu=ultrasparc2
  elif /usr/sbin/prtconf | grep UltraSPARC >/dev/null 2>&1; then
    gmxcpu=ultrasparc
  else
    gmxcpu=""
  fi
  ;;
*)
  gmxcpu=""
  ;;

esac
])



###############################################################
# Macro modified from the fftw distribution (www.fftw.org)
# to determine optimization flags.
# Note that we have modified config.guess and config.sub
# to provide extended information on the detailed type of CPU.
# In general we assume you have recent versions of the compilers
# that support the highest optimization we know of. If not, you 
# can always override these flags, but it's better to upgrade :-)
###############################################################
AC_DEFUN([ACX_COMPILER_MAXOPT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_F77])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine "good" native compiler flags if none specified on command
# line. To avoid repeating the entire procedure for fortran flags, we first
# determine our suggested choices for both C and fortran, and then possibly
# override them with user choices.

cc_vendor="unknown"

case "${host_cpu}-${host_os}" in

  *-solaris2*) 
    case "${gmxcpu}" in
      ultrasparc3*)
        xCFLAGS="-fast -xO5 -xtarget=ultra3 -fsimple=2 -fnonstd -dalign"
        xFFLAGS=$xCFLAGS
        ;;
      ultrasparc2i*)
        xCFLAGS="-fast -xO5 -xtarget=ultra2i -fsimple=2 -fnonstd -dalign"
        xFFLAGS=$xCFLAGS
        ;;
      ultrasparc2*)
        xCFLAGS="-fast -xO5 -xtarget=ultra2 -fsimple=2 -fnonstd -dalign"
        xFFLAGS=$xCFLAGS
        ;;
      ultrasparc*)
        xCFLAGS="-fast -xO5 -xtarget=ultra -fsimple=2 -fnonstd -dalign"
        xFFLAGS=$xCFLAGS
        ;;
      *)
        xCFLAGS="-native -fast -xO5 -fsimple=2 -fnonstd -dalign"
        xFFLAGS=$xCFLAGS
        ;;
    esac
    ;;

  *-hpux*)  
    xCFLAGS="-Ae +O3 +Oall"
    xFFLAGS=$xCFLAGS
    # If you haven't noticed, we don't like hp very much...
    # but perhaps that will change if they make something nice out of ia64.
    ;;

  ia64*-*)
    # The GNU compilers are checked outside this case statement.
    # Check for Intel Compilers. The SGI one was killed before
    # it went final, so I cant imagine anyone is using it...

    # Apparently, -O2 is better than -O3 for villin at least,
    # but I have not yet had time to test all the other benchmarks
    # on both optimization levels. Might need further tweaking.

    # The Intel compilers are _really_ chatty when it comes to
    # warnings, and also echo a lot of incomprehensible internal
    # stuff (not gromacs-related) when we are using ia64 assembly.
    # For this reason we disable warnings...

   if $CC -V 2>&1 | grep 'Intel' > /dev/null 2>&1; then
     xCFLAGS="-O3 -w"
     xASFLAGS=$xCFLAGS
     ac_cv_prog_gcc="no"	
   fi  
   if $F77 -V 2>&1 | grep 'Intel' > /dev/null 2>&1; then
     xFFLAGS="-O3 -w90 -w95 -w"
     ac_cv_prog_g77="no"
   fi  
   # PORTME 2. Check for intel compilers when we get our hands on one!
   ;;	
  *-aix*)
    # dont use inter-procedure analysis for the innerloops - they take
    # forever to compile with it, and it doesnt help at all.

    # use 8 segments (max 2Gb) instead of 1 (max 256Meg) by default.
    xLDFLAGS="$xLDFLAGS -bmaxdata:0x80000000"
    case "${gmxcpu}" in
      power5*)
        xCFLAGS="-O3 -qarch=pwr5 -qtune=pwr5 -qmaxmem=16384"
        xFFLAGS="-O3 -Q -qarch=pwr5 -qtune=pwr5 -qmaxmem=16384 -qhot -qnoipa"
        ;;
      power4*)
	xCFLAGS="-O3 -qarch=pwr4 -qtune=pwr4 -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr4 -qtune=pwr4 -qmaxmem=16384 -qhot -qnoipa"
	;;
      power3*)
	xCFLAGS="-O3 -qarch=pwr3 -qtune=pwr3 -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr3 -qtune=pwr3 -qmaxmem=16384 -qhot -qnoipa"
	;;
      power2*)
	xCFLAGS="-O3 -qarch=pwr2 -qtune=pwr2 -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr2 -qtune=pwr2 -qmaxmem=16384 -qhot -qnoipa"
	;;
      power)
	xCFLAGS="-O3 -qarch=pwr -qtune=pwr -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr -qtune=pwr -qmaxmem=16384 -qhot -qnoipa"
	;;
      ppc604)
	xCFLAGS="-O3 -qarch=604 -qtune=604 -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=604 -qtune=604 -qmaxmem=16384 -qhot"
	;;
      ppc603)
	xCFLAGS="-O3 -qarch=603 -qtune=603 -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=603 -qtune=603 -qmaxmem=16384 -qhot"
	;;
      rs64a)
	xCFLAGS="-O3 -qarch=rs64a -qtune=rs64a -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=rs64a -qtune=rs64a -qmaxmem=16384 -qhot"
	;;
      rs64b)
	xCFLAGS="-O3 -qarch=rs64b -qtune=rs64b -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=rs64b -qtune=rs64b -qmaxmem=16384 -qhot"
	;;
      rs64c)
	xCFLAGS="-O3 -qarch=rs64c -qtune=rs64c -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=rs64c -qtune=rs64c -qmaxmem=16384 -qhot"
	;;
      *)
	xCFLAGS="-O3 -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qmaxmem=16384 -qhot"
	;;
    esac
    ;;

  powerpc*-darwin* | powerpc*-linux* )
    # Check for IBM compilers on OS X     
    if $CC 2>&1 | grep 'IBM' > /dev/null 2>&1; then
       xCFLAGS="-O4 -Q=500 -qaltivec -qnoipa"
    fi
    if $F77 -V 2>&1 | grep 'IBM' > /dev/null 2>&1; then
      xFFLAGS="-O4 -Q=500 -qnoipa"
    fi
    ;;

  mips*-irix*)
    xCFLAGS="-O3 -OPT:IEEE_arithmetic=3 -OPT:rsqrt=ON -SWP:loop_overhead -INLINE:=ON -LNO:opt=1 -LNO:ou_further=3 -OPT:Olimit=0:roundoff=3:alias=typed -woff 1174 -D__INLINE_INTRINSICS"
    xFFLAGS="-O3 -OPT:IEEE_arithmetic=3 -OPT:rsqrt=ON -SWP:loop_overhead -INLINE:=ON -LNO:opt=1 -LNO:ou_further=3 -OPT:Olimit=0:roundoff=3:alias=typed -OPT:cray_ivdep=TRUE"
    
    if $CC -version | grep "Version 7.1" > /dev/null 2>&1; then
      xCFLAGS="$xCFLAGS -GCM:aggressive_speculation -GCM:array_speculation" 
      xFFLAGS="$xFFLAGS -GCM:aggressive_speculation -GCM:array_speculation" 
    fi

    if $CC -version | grep "Version 7.3" > /dev/null 2>&1; then
      xCFLAGS="$xCFLAGS -SWP:heur=fdms,nhms,fdnms" 
      xFFLAGS="$xFFLAGS -SWP:heur=fdms,nhms,fdnms" 
    fi
    xLDFLAGS="-woff 84"

    # I have removed -n32 from the flags since it causes too many problems.
    # New SGIs should use the right objects automatically, and it's not
    # worth the hassle for 5-10 year old machines...  

    case "${gmxcpu}" in
      r12000*)
	xCFLAGS="$IRIXOBJFLAG -r12000 -mips4 $xCFLAGS"
	xFFLAGS="$IRIXOBJFLAG -r12000 -mips4 $xFFLAGS"
	xLDFLAGS="$IRIXOBJFLAG -r12000 -mips4 $xLDFLAGS"
	;;
      r10000*)
	xCFLAGS="$IRIXOBJFLAG -r10000 -mips4 $xCFLAGS"
	xFFLAGS="$IRIXOBJFLAG -r10000 -mips4 $xFFLAGS"
	xLDFLAGS="$IRIXOBJFLAG -r10000 -mips4 $xLDFLAGS"
	;;
      r8000*)
	xCFLAGS="$IRIXOBJFLAG -r8000 -mips4 $xCFLAGS"
	xFFLAGS="$IRIXOBJFLAG -r8000 -mips4 $xFFLAGS"
	xLDFLAGS="$IRIXOBJFLAG -r8000 -mips4 $xLDFLAGS"
	;;
      r5000*)
	xCFLAGS="$IRIXOBJFLAG -r5000 -mips4 $xCFLAGS"
	xFFLAGS="$IRIXOBJFLAG -r5000 -mips4 $xFFLAGS"
	xLDFLAGS="$IRIXOBJFLAG -r5000 -mips4 $xLDFLAGS"
	;;
      *)		
	xCFLAGS="$IRIXOBJFLAG $xCFLAGS"
	xFFLAGS="$IRIXOBJFLAG $xFFLAGS"
	xLDFLAGS="$IRIXOBJFLAG $xLDFLAGS"
	;;
    esac
    ;;

  alpha*-osf*) 
     # NB: -arch implies -tune according to the cc manual.
     # We dont use -ifo since it conflicts with dependency
     # generation on old versions of the compiler.
    case "${host_cpu}" in
      alphaev*)
        # extract the processor from cpu type (e.g. alphaev56 -> ev56)
        evtype=`echo ${host_cpu} | sed 's/alpha//'`
        xCFLAGS="-std1 -fast -O4 -no_ifo -arch $evtype -unroll 2 -fp_reorder"
        xFFLAGS="$xCFLAGS -assume noaccuracy_sensitive"
	xASFLAGS="-O4 -no_ifo -arch $evtype"
        xLDFLAGS="-O4"
        ;;
      *)
	xCFLAGS="-std1 -fast -O4 -no_ifo -arch host -unroll 2 -fp_reorder"
	xFFLAGS="$xCFLAGS -assume noaccuracy_sensitive"
	xASFLAGS="-O4 -no_ifo -arch host"
	xLDFLAGS="-O4"
	;;
    esac
    ;;

  alpha*-linux*)
    case "${host_cpu}" in
      alphaev*)
	# extract the processor from cpu type (e.g. alphaev56 -> ev56)
	evtype=`echo ${host_cpu} | sed 's/alpha//'`
	tmpCFLAGS="-std1 -fast -O4 -no_ifo -arch $evtype -unroll 2 -fp_reorder"
	tmpFFLAGS="$tmpCFLAGS -assume noaccuracy_sensitive"
	tmpASFLAGS="-O4 -no_ifo -arch $evtype"
	tmpLDFLAGS="-O4"
	;;
      *)
	tmpCFLAGS="-std1 -fast -O4 -no_ifo -arch host -unroll 2 -fp_reorder"
	tmpFFLAGS="$tmpCFLAGS -assume noaccuracy_sensitive"
	tmpASFLAGS="-O4 -no_ifo -arch host"
	tmpLDFLAGS="-O4"
	;;
    esac
	# Compaq sometimes uses -version and sometimes -V
	# Not 100% sure if ccc always has -V and F77 -version, so 
	# we check both alternatives to be sure.
    if (($CC -V 2>&1 | grep ompaq > /dev/null) || 
	($CC -version 2>&1 | grep ompaq > /dev/null)); then
      xCFLAGS="$tmpCFLAGS"
      xASFLAGS="$tmpASFLAGS"
      cc_vendor="Compaq"
    fi
    if test "$enable_fortran" = "yes"; then
      if (($F77 -V 2>&1 | grep ompaq > /dev/null) || 
	  ($F77 -version 2>&1 | grep ompaq > /dev/null)); then
        xFFLAGS="$tmpFFLAGS"
      fi
    fi
    ;;

  *-*)
    # most of these systems (e.g. linux, FreeBSD) use gcc which is treated
    # further down, but check for some specific compilers.
    # Portland group compilers:
    if $CC -V 2>  /dev/null | grep ortland > /dev/null 2>&1; then
      case "${host_cpu}" in
	i586)
	  pgiopt="-tp p5" 
          ;;
	i686)
	  pgiopt="-tp p6" 
 	  ;;
      esac
      xCFLAGS="$pgiopt -fast -pc 32"
      xASFLAGS="$xCFLAGS"
    fi
    if test "$enable_fortran" = "yes"; then
      if $F77 -version 2>  /dev/null | grep Portland > /dev/null 2>&1; then
	xFFLAGS="$xCFLAGS"
      fi	
    fi

    # Intel compilers
    # The Intel compilers are _really_ chatty when it comes to
    # warnings, and also echo a lot of incomprehensible internal
    # stuff (not gromacs-related) when we are using assembly.
    # For this reason we disable warnings...

    if $CC -V 2>&1 | grep 'Intel' > /dev/null 2>&1; then
      ac_cv_prog_gcc="no"	
      case "${host_cpu}" in
        x86_64)
          xCFLAGS="-O3 -tpp7 -axW -ip -w"
          ;;
	i686)
	  xCFLAGS="-O3 -tpp6 -axK -ip -w" 
 	  ;;
	ia64)
	  xCFLAGS="-O3 -ip -w" 
 	  ;;
      esac
      xASFLAGS="$xCFLAGS"
      # search in /usr/local/include too, just as gcc does. (handy for fftw)
      CPPFLAGS="$CPPFLAGS -I/usr/local/include"
    fi
    if test "$enable_fortran" = "yes"; then
      if $F77 -V 2>&1 | grep 'Intel' > /dev/null 2>&1; then
        ac_cv_prog_g77="no"
	xFFLAGS="$xCFLAGS -w90 -w95"
      fi	
    fi
	
    ;;
esac	
# Phew, end of all those operating systems and processors!			

# use default flags for gcc/g77 on all systems
if test $ac_cv_prog_gcc = yes; then
  ACX_CHECK_CC_FLAGS(-O3,o3,xCFLAGS="$xCFLAGS -O3")
  xCFLAGS="$xCFLAGS -fomit-frame-pointer -finline-functions -Wall -Wno-unused"
  # For alpha axp assembly we need the preprocessor to tell elf from ecoff.
  # The compaq ccc compiler only knows .s files, and always runs them
  # through cpp. We support this by telling gcc to preprocess .s files.
  case "${host_cpu}" in
    alphaev*)
      xASFLAGS="$xCFLAGS -x assembler-with-cpp"
      ;;
    *)
      ;;
  esac
  # -malign-double for x86 systems
  ACX_CHECK_CC_FLAGS(-malign-double,align_double,xCFLAGS="$xCFLAGS -malign-double")
fi
  
if test $enable_fortran = yes; then
  if test $ac_cv_prog_g77 = yes; then
    xFFLAGS="-O3 -ffast-math -fomit-frame-pointer -finline-functions -funroll-all-loops -Wall -Wno-unused"
    # -malign-double for f77 on x86 systems - haven't checked that this works yet.
    #ACX_CHECK_F77_FLAGS(-malign-double,align_double,xFFLAGS="$xFFLAGS -malign-double")
  fi
fi

CPU_FLAGS=""

if test "$GCC" = "yes"; then
  # try to guess correct CPU flags, at least for powerpc linux
  case "${host_cpu}" in
    # i586/i686 cpu flags don't improve speed, thus no need to use them.
    # don't check f77 separately - we assume f77 and gcc are similar	  
    powerpc*)
        # don't use the separate apple cpp on OS X
#        ACX_CHECK_CC_FLAGS(-no-cpp-precomp,no_cpp_precomp,xCFLAGS="$xCFLAGS -no-cpp-precomp")
        if test "$enable_ppc_altivec" = "yes"; then
            # Apple (darwin) uses a hacked version of gcc with special flags 
            case "${host_os}" in
            darwin*)       	            	
                ACX_CHECK_CC_FLAGS(-faltivec,faltivec,xCFLAGS="$xCFLAGS -faltivec")
                ;;
            *)
                # Need to update CPPFLAGS too, since we later call 
                # AC_CHECK_HEADER for altivec.h, and then autoconf complains
                # if it cannot process it with the preprocessor.
                ACX_CHECK_CC_FLAGS(-maltivec,maltivec,xCFLAGS="$xCFLAGS -maltivec" CPPFLAGS="$CPPFLAGS -maltivec")
                ACX_CHECK_CC_FLAGS(-mabi=altivec,mabialtivec,xCFLAGS="$xCFLAGS -mabi=altivec" CPPFLAGS="$CPPFLAGS -mabi=altivec")
                ;;
            esac 
        fi
        # -funroll-all-loops exposes a bug in altivec-enabled gcc-2.95.3
        # on powerpc, so we only enable it on other platforms or gcc3.    
        # The gcc 2.95 instruction scheduler also destroys our handcoded altivec,
        # so disable instruction scheduling on 2.95
        if $CC --version 2>&1 | grep '2.95' > /dev/null 2>&1; then
	  echo "*****************************************************************************"
          echo "* IMPORTANT INFO: You are using gcc-2.95.x on PowerPC. This compiler works, *"
          echo "* but you will get better performance with gcc-3.3 or later. If you are     *"
          echo "* running OS X, download the latest devtools from http://developer.apple.com*"
	  echo "*****************************************************************************"
          ACX_CHECK_CC_FLAGS(-fno-schedule-insns,fno_schedule_insns,xCFLAGS="$xCFLAGS -fno-schedule-insns")
        fi
	ACX_CHECK_CC_FLAGS(-mcpu=7450,m_cpu_7450,CPU_FLAGS="-mcpu=7450")
	ACX_CHECK_CC_FLAGS(-mtune=970,m_tune_970,CPU_FLAGS="$CPU_FLAGS -mtune=970")
	if test -z "$CPU_FLAGS"; then
  	  ACX_CHECK_CC_FLAGS(-mcpu=powerpc,m_cpu_powerpc,CPU_FLAGS="-mcpu=powerpc")
        fi	
      ;;
   *)
        ACX_CHECK_CC_FLAGS(-funroll-all-loops,funroll_all_loops,xCFLAGS="$xCFLAGS -funroll-all-loops")
      ;;
   esac
fi

if test -n "$CPU_FLAGS"; then
  xCFLAGS="$xCFLAGS $CPU_FLAGS"
  xFFLAGS="$xFFLAGS $CPU_FLAGS"
  xASFLAGS="$xASFLAGS $CPU_FLAGS"
fi

# Now check if the user provided anything special for C or fortran...
# Not nice to have checked everything then, but otherwise we would have
# to use entirely separate checks for C and fortran flags, doubling the code.
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS="$xCFLAGS"
  # Use the extra link optimization flags on e.g. irix only when
  # we are using our own C compiler flags
  LDFLAGS="$LDFLAGS $xLDFLAGS"
  
  if test -z "$CFLAGS"; then
    echo "********************************************************************"
    echo "* Note: We have not optimized the C compiler flags on your target  *"
    echo "* yet, but the default CFLAGS=-O3 should be OK in most cases.      *"
    echo "* You can override this by setting the CFLAGS environment variable.*"
    echo "*******************************************************************"
    CFLAGS="-O3"
  fi
  ACX_CHECK_CC_FLAGS(${CFLAGS}, guessed_cflags, , [
    echo "*******************************************************************"
    echo "* Sorry, these optimization settings don't seem to work for       *"
    echo "* your C compiler. Use make CFLAGS=..., or edit the top Makefile. *"
    echo "*******************************************************************"
    CFLAGS=""
  ])
else
  echo "******************************************"
  echo "* Using CFLAGS from environment variable *"
  echo "******************************************"
fi

if test "$enable_fortran" = "yes"; then	
  if test "$ac_test_FFLAGS" != "set"; then
    FFLAGS="$xFFLAGS"
    if test -z "$FFLAGS"; then
    echo "********************************************************************"
    echo "* Note: We have not optimized the Fortran compiler flags on your   *"
    echo "* target, but the default FFLAGS=-O3 should be OK in most cases.   *"
    echo "* You can override this by setting the CFLAGS environment variable.*"
    echo "********************************************************************"
      FFLAGS="-O3"
    fi
    ACX_CHECK_F77_FLAGS(${FFLAGS}, guessed_fflags, , [
      echo "*******************************************************************"
      echo "* Sorry, these optimization settings don't seem to work for       *"
      echo "* your f77 compiler. Use make FFLAGS=.., or edit the top Makefile.*"
      echo "*******************************************************************"
      FFLAGS=""
    ])
  else
    echo "******************************************"
    echo "* Using FFLAGS from environment variable *"
    echo "******************************************"
  fi
fi
# Be silent for assembly flags, they are usually not important anyway
if test "${ASFLAGS+set}" != set; then
  if test "${xASFLAGS+set}" != set; then
    xASFLAGS="$CFLAGS"
  fi
  ASFLAGS="$xASFLAGS"
fi

])


dnl @synopsis ACX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro figures out how to build C programs using POSIX
dnl threads. It sets the PTHREAD_LIBS output variable to the threads
dnl library and linker flags, and the PTHREAD_CFLAGS output variable
dnl to any special C compiler flags that are needed. (The user can also
dnl force certain compiler flags/libs to be tested by setting these
dnl environment variables.)
dnl
dnl Also sets PTHREAD_CC to any special C compiler that is needed for
dnl multi-threaded programs (defaults to the value of CC otherwise).
dnl (This is necessary on AIX to use the special cc_r compiler alias.)
dnl
dnl If you are only building threads programs, you may wish to
dnl use these variables in your default LIBS, CFLAGS, and CC:
dnl
dnl LIBS="$PTHREAD_LIBS $LIBS"
dnl CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
dnl CC="$PTHREAD_CC"
dnl
dnl In addition, if the PTHREAD_CREATE_JOINABLE thread-attribute
dnl constant has a nonstandard name, defines PTHREAD_CREATE_JOINABLE
dnl to that name (e.g. PTHREAD_CREATE_UNDETACHED on AIX).
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a threads
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found. If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_PTHREAD.
dnl
dnl Please let the authors know if this macro fails on any platform,
dnl or if you have any other suggestions or comments. This macro was
dnl based on work by SGJ on autoconf scripts for FFTW (www.fftw.org)
dnl (with help from M. Frigo), as well as ac_pthread and hb_pthread
dnl macros posted by AFC to the autoconf macro repository. We are also
dnl grateful for the helpful feedback of numerous users.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu> and Alejandro Forero Cuervo <bachue@bachue.com>

AC_DEFUN([ACX_PTHREAD], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG_C
acx_pthread_ok=no

# We used to check for pthread.h first, but this fails if pthread.h
# requires special compiler flags (e.g. on True64 or Sequent).
# It gets checked for in the link test anyway.

# First of all, check if the user has set any of the PTHREAD_LIBS,
# etcetera environment variables, and if threads linking works using
# them:
if test x"$PTHREAD_LIBS$PTHREAD_CFLAGS" != x; then
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        AC_MSG_CHECKING([for pthread_join in LIBS=$PTHREAD_LIBS with CFLAGS=$PTHREAD_CFLAGS])
        AC_TRY_LINK_FUNC(pthread_join, acx_pthread_ok=yes)
        AC_MSG_RESULT($acx_pthread_ok)
        if test x"$acx_pthread_ok" = xno; then
                PTHREAD_LIBS=""
                PTHREAD_CFLAGS=""
        fi
        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"
fi

# We must check for the threads library under a number of different
# names; the ordering is very important because some systems
# (e.g. DEC) have both -lpthread and -lpthreads, where one of the
# libraries is broken (non-POSIX).

# Create a list of thread flags to try. Items starting with a "-" are
# C compiler flags, and other items are library names, except for "none"
# which indicates that we try without any flags at all.

acx_pthread_flags="pthreads none -Kthread -kthread lthread -pthread -pthreads -mthreads pthread --thread-safe -mt"

# The ordering *is* (sometimes) important. Some notes on the
# individual items follow:

# pthreads: AIX (must check this before -lpthread)
# none: in case threads are in libc; should be tried before -Kthread and
# other compiler flags to prevent continual compiler warnings
# -Kthread: Sequent (threads in libc, but -Kthread needed for pthread.h)
# -kthread: FreeBSD kernel threads (preferred to -pthread since SMP-able)
# lthread: LinuxThreads port on FreeBSD (also preferred to -pthread)
# -pthread: Linux/gcc (kernel threads), BSD/gcc (userland threads)
# -pthreads: Solaris/gcc
# -mthreads: Mingw32/gcc, Lynx/gcc
# -mt: Sun Workshop C (may only link SunOS threads [-lthread], but it
# doesn't hurt to check since this sometimes defines pthreads too;
# also defines -D_REENTRANT)
# pthread: Linux, etcetera
# --thread-safe: KAI C++

case "${host_cpu}-${host_os}" in
        *solaris*)

        # On Solaris (at least, for some versions), libc contains stubbed
        # (non-functional) versions of the pthreads routines, so link-based
        # tests will erroneously succeed. (We need to link with -pthread or
        # -lpthread.) (The stubs are missing pthread_cleanup_push, or rather
        # a function called by this macro, so we could check for that, but
        # who knows whether they'll stub that too in a future libc.) So,
        # we'll just look for -pthreads and -lpthread first:

        acx_pthread_flags="-pthread -pthreads pthread -mt $acx_pthread_flags"
        ;;
esac

if test x"$acx_pthread_ok" = xno; then
for flag in $acx_pthread_flags; do

        case $flag in
                none)
                AC_MSG_CHECKING([whether pthreads work without any flags])
                ;;

                -*)
                AC_MSG_CHECKING([whether pthreads work with $flag])
                PTHREAD_CFLAGS="$flag"
                ;;

                *)
                AC_MSG_CHECKING([for the pthreads library -l$flag])
                PTHREAD_LIBS="-l$flag"
                ;;
        esac

        save_LIBS="$LIBS"
        save_CFLAGS="$CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Check for various functions. We must include pthread.h,
        # since some functions may be macros. (On the Sequent, we
        # need a special flag -Kthread to make this header compile.)
        # We check for pthread_join because it is in -lpthread on IRIX
        # while pthread_create is in libc. We check for pthread_attr_init
        # due to DEC craziness with -lpthreads. We check for
        # pthread_cleanup_push because it is one of the few pthread
        # functions on Solaris that doesn't have a non-functional libc stub.
        # We try pthread_create on general principles.
        AC_TRY_LINK([#include <pthread.h>],
                    [pthread_t th; pthread_join(th, 0);
                     pthread_attr_init(0); pthread_cleanup_push(0, 0);
                     pthread_create(0,0,0,0); pthread_cleanup_pop(0); ],
                    [acx_pthread_ok=yes])

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        AC_MSG_RESULT($acx_pthread_ok)
        if test "x$acx_pthread_ok" = xyes; then
                break;
        fi

        PTHREAD_LIBS=""
        PTHREAD_CFLAGS=""
done
fi

# Various other checks:
if test "x$acx_pthread_ok" = xyes; then
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Detect AIX lossage: threads are created detached by default
        # and the JOINABLE attribute has a nonstandard name (UNDETACHED).
        AC_MSG_CHECKING([for joinable pthread attribute])
        AC_TRY_LINK([#include <pthread.h>],
                    [int attr=PTHREAD_CREATE_JOINABLE;],
                    ok=PTHREAD_CREATE_JOINABLE, ok=unknown)
        if test x"$ok" = xunknown; then
                AC_TRY_LINK([#include <pthread.h>],
                            [int attr=PTHREAD_CREATE_UNDETACHED;],
                            ok=PTHREAD_CREATE_UNDETACHED, ok=unknown)
        fi
        if test x"$ok" != xPTHREAD_CREATE_JOINABLE; then
                AC_DEFINE(PTHREAD_CREATE_JOINABLE, $ok,
                          [Define to the necessary symbol if this constant
                           uses a non-standard name on your system.])
        fi
        AC_MSG_RESULT(${ok})
        if test x"$ok" = xunknown; then
                AC_MSG_WARN([we do not know how to create joinable pthreads])
        fi

        AC_MSG_CHECKING([if more special flags are required for pthreads])
        flag=no
        case "${host_cpu}-${host_os}" in
                *-aix* | *-freebsd*) flag="-D_THREAD_SAFE";;
                *solaris* | alpha*-osf*) flag="-D_REENTRANT";;
        esac
        AC_MSG_RESULT(${flag})
        if test "x$flag" != xno; then
                PTHREAD_CFLAGS="$flag $PTHREAD_CFLAGS"
        fi

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        # More AIX lossage: must compile with cc_r
        AC_CHECK_PROG(PTHREAD_CC, cc_r, cc_r, ${CC})
else
        PTHREAD_CC="$CC"
fi

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_CC)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pthread_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.]),[$1])
        :
else
        acx_pthread_ok=no
        $2
fi
AC_LANG_RESTORE
])dnl ACX_PTHREAD 




