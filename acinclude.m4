 
# ACX_CHECK_FFTW()
# ----------------
# This macro checks for fftw header files and libraries,
# including the possible prefixing with s or d to determine precision.
# Arg 1 is the fftw header/library name to check for, without 
# prefix or anything else (e.g. rfftw_mpi for real MPI transforms)
# Arg 2 is the floating point size we want (real=4, double=8)
AC_DEFUN(ACX_CHECK_FFTW,
[
if test -z "$ac_fftw_firstname"; then

if test "$2" = "8"; then
  prec="double"
  fftwcheckprefix=d
else
  prec="single"
  fftwcheckprefix=s
fi

usedprefix=""
ok="no"
# check header doesn't work, since we must use mpicc to get includes, not just /lib/cpp
AC_MSG_CHECKING([for $1.h])
AC_TRY_COMPILE([#include <$1.h>],,
[
fftwname=$1 
AC_MSG_RESULT(yes)
],
AC_MSG_RESULT(no))


if test -n "$fftwname"; then
# we cannot run the code since MPI program might not be allowed outside a charge queue
AC_TRY_COMPILE([#include <$fftwname.h>],
[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $2)]; ],
ok="yes",ok="no")
fi

fftwname=$1

if test "$ok" != "yes"; then
  xfftwname=${fftwcheckprefix}${fftwname}
  AC_MSG_CHECKING([for $xfftwname.h])
  AC_TRY_COMPILE([#include <$xfftwname.h>],,AC_MSG_RESULT(yes),
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Cannot find any $prec precision $fftwname.h or $xfftwname.h]
[Do you have $prec precision FFTW installed? You can find it at www.fftw.org] 
[Note that the default FFTW setup is double precision. You change the]
[FFTW configuration to single with --enable-float and turn on MPI support]
[with --enable-mpi. It is a good idea to install both single & double.] 
[If your sysadm doesn't want to install it you can do it to a location]
[in your home directory and provide Gromacs configure with the correct]
[paths by setting the CPPFLAGS and LDFLAGS environment variables.]
[Check the Gromacs INSTALL file for additional information.])
])
AC_TRY_COMPILE([#include <$xfftwname.h>],[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $2)];],
[
fftwname=$xfftwname 
usedprefix=$fftwcheckprefix
],
AC_MSG_ERROR([Cannot find any $prec precision $fftwname.h or $xfftwname.h]
[Do you have $prec precision FFTW installed? You can find it at www.fftw.org] 
[Note that the default FFTW setup is double precision. You change the]
[FFTW configuration to single with --enable-float and turn on MPI support]
[with --enable-mpi. It is a good idea to install both single & double.] 
[If your sysadm doesn't want to install it you can do it to a location]
[in your home directory and provide Gromacs configure with the correct]
[paths by setting the CPPFLAGS and LDFLAGS environment variables.]
[Check the Gromacs INSTALL file for additional information.]))
fi

AC_CHECK_LIB($fftwname,main,,AC_MSG_ERROR([Can't find a library to match the $fftwname header]))
ac_fftw_savedprefix=$usedprefix
ac_fftw_firstname=$fftwname

else

fftwname=${ac_fftw_savedprefix}$1
AC_MSG_CHECKING([for $fftwname.h])
AC_TRY_COMPILE(
[#include <$fftwname.h>],,AC_MSG_RESULT(yes)
AC_CHECK_LIB($fftwname,main,,AC_MSG_ERROR([Can't find a library to match the $fftwname header])),
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Cant find $fftwname.h header. Make sure all your fftw prefixes match - we already use $ac_fftw_firstname.h])
])

fi

])



dnl macro from the fftw distribution (http://www.fftw.org)
dnl like AC_SUBST, but replace XXX_variable_XXX instead of @variable@
dnl This macro protects VARIABLE from being diverted twice
dnl if this macro is called twice for it.
dnl AC_SUBST(VARIABLE)
define(ACX_SUBST_XXX,
[ifdef([ACX_SUBST_XXX_$1], ,
[define([ACX_SUBST_XXX_$1], )dnl
AC_DIVERT_PUSH(AC_DIVERSION_SED)dnl
s=XXX_$1_XXX=[$]$1=g
AC_DIVERT_POP()dnl
])])


# ACX_F77_NAME_MANGLING
# ---------------------
# This is a macro adopted from the fftw distribution (http://www.fftw.org)
# to determine how we should call fortran functions from C.
AC_DEFUN(ACX_F77_NAME_MANGLING,
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_F77])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_MSG_CHECKING(fortran name mangling)
cat > mangle-func.f <<EOF
      subroutine foobar()
      return
      end
      subroutine foo_bar()
      return
      end
EOF
ac_try='$F77 -c $FFLAGS mangle-func.f 1>&AC_FD_CC'
if AC_TRY_EVAL(ac_try); then
  ac_try=""
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat mangle-func.f >&AC_FD_CC
  rm -f mangle-func*
  AC_MSG_ERROR(failed to compile fortran test program)
fi

ac_f77_mangle_type=unknown
AC_LANG_SAVE
AC_LANG_C
ac_save_LIBS="$LIBS"
LIBS="mangle-func.o $FLIBS $LIBS"
AC_TRY_LINK(,foobar();,
     ac_f77_mangle_type=lowercase,
     AC_TRY_LINK(,foobar_();,
          ac_f77_mangle_type=lowercase-underscore,
          AC_TRY_LINK(,FOOBAR();,
               ac_f77_mangle_type=uppercase,
               AC_TRY_LINK(,FOOBAR_();,
                    ac_f77_mangle_type=uppercase-underscore))))
LIBS="$ac_save_LIBS"
AC_LANG_RESTORE
AC_MSG_RESULT($ac_f77_mangle_type)

mangle_try=unknown
case $ac_f77_mangle_type in
        lowercase)
                AC_DEFINE(F77_NAME_LOWERCASE,,[call f77 with lowercase, no underscore])
                mangle_try=foo_bar_
                ;;
        lowercase-underscore)
                AC_DEFINE(F77_NAME_LOWERCASE_UNDERSCORE,,[call f77 with lowercase and underscore])
                mangle_try=foo_bar__
                ;;
        uppercase)
                AC_DEFINE(F77_NAME_UPPERCASE,,[call f77 with uppercase, no underscore])
                mangle_try=FOO_BAR_
                ;;
        uppercase-underscore)
                AC_DEFINE(F77_NAME_UPPERCASE_UNDERSCORE,,[call f77 with uppercase and underscore])
                mangle_try=FOO_BAR__
                ;;
esac

AC_MSG_CHECKING(whether f77 functions with underscore get an extra underscore)

AC_LANG_SAVE
AC_LANG_C
ac_save_LIBS="$LIBS"
LIBS="mangle-func.o $FLIBS $LIBS"
AC_TRY_LINK(,$mangle_try();,
            [ac_f77_mangle_underscore=yes;
             AC_DEFINE(F77_NAME_EXTRA_UNDERSCORE,,[Append extra underscore to already underscored names])],
            [ac_f77_mangle_underscore=no])
LIBS="$ac_save_LIBS"
AC_LANG_RESTORE
rm -f mangle-func*
AC_MSG_RESULT($ac_f77_mangle_underscore)
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
AC_DEFUN(AC_FIND_MOTIF,
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
  motif_includes=none
  motif_libraries=none
fi

AC_ARG_WITH(motif-includes,
[  --with-motif-includes=DIR    Motif include files are in DIR],
motif_includes="$withval")

AC_ARG_WITH(motif-libraries,
[  --with-motif-libraries=DIR   Motif libraries are in DIR],
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
ac_cv_motif_includes="none"
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
ac_cv_motif_libraries="none"
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
if test "$motif_includes" = "none" -o "$motif_libraries" = "none"; then
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
test "$motif_libraries_result" = "none" && motif_libraries_result="(none)"
test "$motif_includes_result" = "none" && motif_includes_result="(none)"
AC_MSG_RESULT(
  [libraries $motif_libraries_result, headers $motif_includes_result])
])dnl


dnl macro modified from the fftw distribution (www.fftw.org)
AC_DEFUN(ACX_CHECK_CC_FLAGS,
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether $CC accepts $1, ac_$2,
[echo 'void f(){}' > conftest.c
if test -z "`$CC $1 -c conftest.c 2>&1`"; then
	ac_$2=yes
else
	ac_$2=no
fi
rm -f conftest*
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
AC_DEFUN(ACX_CHECK_F77_FLAGS,
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
rm -f conftest*
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
AC_DEFUN(ACX_DETECT_GMXCPU,
[
AC_REQUIRE([AC_CANONICAL_HOST])

#
# Determine the exact cpu type on some common systems where it is 
# not visible from the host triplet.
# (on e.g. intel and dec/tru64 the host type is enough)

case "${host_cpu}-${host_os}" in

rs6000*-aix*)
  # we need to fool the combination of m4, sh and awk - thus the seemingly unnecessary n
  IBM_CPU_ID=`/usr/sbin/lsdev -C -c processor -S available | head -1 | awk '{ n=1; print $n }'`
  if /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER4 >/dev/null 2>&1; then
    gmxcpu=power4
  elif /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER3 >/dev/null 2>&1; then
    gmxcpu=power3
  elif /usr/sbin/lsattr -EHl ${IBM_CPU_ID} | grep POWER2 >/dev/null 2>&1; then
    gmxcpu=power2
  else
    gmxcpu=""
  fi
  ;;

powerpc*-aix*)
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
  else
    gmxcpu=""
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
AC_DEFUN(ACX_COMPILER_MAXOPT,
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_F77])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine "good" native compiler flags if none specified on command
# line. To avoid repeating the entire procedure for fortran flags, we first
# determine our suggested choices for both C and fortran, and then possibly
# override them with user choices.

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

  rs6000*-aix*)
    # dont use inter-procedure analysis for the innerloops - they take
    # forever to compile with it, and it doesnt help at all.
    case "${gmxcpu}" in
      power4*)
	xCFLAGS="-O3 -qarch=pwr4 -qtune=pwr4 -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr4 -qtune=pwr4 -qmaxmem=16384 -qhot -qnoipa"
	;;
      power3*)
	xCFLAGS="-O3 -qarch=pwr3 -qtune=pwr3 -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr3 -qtune=pwr3 -qmaxmem=16384 -qhot -qnoipa"
	;;
      power2*)
	xCFLAGS="-O3 -qarch=pwr2 -qtune=pwr2 -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr2 -qtune=pwr2 -qmaxmem=16384 -qhot -qnoipa"
	;;
      power)
	xCFLAGS="-O3 -qarch=pwr -qtune=pwr -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=pwr -qtune=pwr -qmaxmem=16384 -qhot -qnoipa"
	;;
      *)
	# I don't think people are using anything older than power2, so we tune for
        # pwr, but dont set the arch since it is nice to have common binaries 
        # that run also on powerpc.
	xCFLAGS="-O3 -qlanglvl=ansi -qtune=pwr -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qtune=pwr -qmaxmem=16384 -qhot"
	;;
    esac
    ;;

  powerpc*-aix*)
    case "${gmxcpu}" in
      ppc604)
	xCFLAGS="-O3 -qarch=604 -qtune=604 -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=604 -qtune=604 -qmaxmem=16384 -qhot"
	;;
      ppc603)
	xCFLAGS="-O3 -qarch=603 -qtune=603 -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=603 -qtune=603 -qmaxmem=16384 -qhot"
	;;
      rs64a)
	xCFLAGS="-O3 -qarch=rs64a -qtune=rs64a -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=rs64a -qtune=rs64a -qmaxmem=16384 -qhot"
	;;
      rs64b)
	xCFLAGS="-O3 -qarch=rs64b -qtune=rs64b -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=rs64b -qtune=rs64b -qmaxmem=16384 -qhot"
	;;
      rs64c)
	xCFLAGS="-O3 -qarch=rs64c -qtune=rs64c -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qarch=rs64c -qtune=rs64c -qmaxmem=16384 -qhot"
	;;
      *)
	xCFLAGS="-O3 -qlanglvl=ansi -qmaxmem=16384"
	xFFLAGS="-O3 -Q -qmaxmem=16384 -qhot"
	;;
    esac
    ;;

  mips*-irix*)
    xCFLAGS="-O3 -OPT:IEEE_arithmetic=3 -OPT:rsqrt=ON -SWP:loop_overhead -INLINE:=ON -LNO:opt=1 -LNO:ou_further=3 -OPT:Olimit=0:roundoff=3:alias=typed -fullwarn -woff 1174 -D__INLINE_INTRINSICS"
    xFFLAGS="-O3 -OPT:IEEE_arithmetic=3 -OPT:rsqrt=ON -SWP:loop_overhead -INLINE:=ON -LNO:opt=1 -LNO:ou_further=3 -OPT:Olimit=0:roundoff=3:alias=typed -OPT:cray_ivdep=TRUE"

    if $CC -version | grep "Version 7.1" > /dev/null 2>&1; then
      xCFLAGS="$xCFLAGS -GCM:aggressive_speculation -GCM:array_speculation" 
      xFFLAGS="$xFFLAGS -GCM:aggressive_speculation -GCM:array_speculation" 
    fi

    if $CC -version | grep "Version 7.3" > /dev/null 2>&1; then
      xCFLAGS="$xCFLAGS -SWP:heur=fdms,nhms,fdnms" 
      xFFLAGS="$xFFLAGS -SWP:heur=fdms,nhms,fdnms" 
    fi

    case "${gmxcpu}" in
      r12000*)
	xCFLAGS="-n32 -r12000 -mips4 $xCFLAGS"
	xFFLAGS="-n32 -r12000 -mips4 $xFFLAGS"
	xLDFLAGS="-n32 -r12000 -mips4"
	;;
      r10000*)
	xCFLAGS="-n32 -r10000 -mips4 $xCFLAGS"
	xFFLAGS="-n32 -r10000 -mips4 $xFFLAGS"
	xLDFLAGS="-n32 -r10000 -mips4"
	;;
      r8000*)
	xCFLAGS="-n32 -r8000 -mips4 $xCFLAGS"
	xFFLAGS="-n32 -r8000 -mips4 $xFFLAGS"
	xLDFLAGS="-n32 -r8000 -mips4"
	;;
      r5000*)
	xCFLAGS="-n32 -r5000 -mips4 $xCFLAGS"
	xFFLAGS="-n32 -r5000 -mips4 $xFFLAGS"
	xLDFLAGS="-n32 -r5000 -mips4"
	;;
      *)		
	xCFLAGS="-n32 $xCFLAGS"
	xFFLAGS="-n32 $xFFLAGS"
	xLDFLAGS="-n32"
	;;
    esac
    ;;

  alpha*-osf*)
    case "${host_cpu}" in
      alphaev*)
        # extract the processor from cpu type (e.g. alphaev56 -> ev56)
        evtype=`echo ${host_cpu} | sed 's/alpha//'`
        xCFLAGS="-O5 -arch $evtype -tune $evtype -fast -unroll 2 -fp_reorder"
        xFFLAGS="$xCFLAGS -assume noaccuracy_sensitive"
        xLDFLAGS="-O4"
        ;;
      *)
	xCFLAGS="-O5 -arch host -tune host -fast -unroll 2 -fp_reorder"
	xFFLAGS="$xCFLAGS -assume noaccuracy_sensitive"
	xLDFLAGS="-O4"
	;;
    esac
    ;;

  alpha*-linux*)
    case "${host_cpu}" in
      alphaev*)
	# extract the processor from cpu type (e.g. alphaev56 -> ev56)
	evtype=`echo ${host_cpu} | sed 's/alpha//'`
	tmpCFLAGS="-O5 -arch $evtype -tune $evtype -fast -unroll 2 -fp_reorder"
	tmpFFLAGS="$tmpCFLAGS -assume noaccuracy_sensitive"
	tmpLDFLAGS="-O4"
	;;
      *)
	tmpCFLAGS="-O5 -arch host -tune host -fast -unroll 2 -fp_reorder"
	tmpFFLAGS="$tmpCFLAGS -assume noaccuracy_sensitive"
	tmpLDFLAGS="-O4"
	;;
    esac
    if $CC -V 2>  /dev/null | grep Compaq > /dev/null 2>&1; then
      xCFLAGS="$tmpCFLAGS"
    fi
    if test "$enable_fortran" = "yes"; then
      if $F77 -V 2>  /dev/null | grep Compaq > /dev/null 2>&1; then
        xFFLAGS="$tmpFFLAGS"
      fi
    fi
    ;;

  *-*)
    # most of these systems (e.g. linux, FreeBSD) use gcc which is treated
    # further down, but we can at least check if the Portland compilers are used:
    if $CC -V 2>  /dev/null | grep Portland > /dev/null 2>&1; then
      case "${host_cpu}" in
	i586)
	  pgiopt="-tp p5" 
          ;;
	i686)
	  pgiopt="-tp p6" 
 	  ;;
      esac
      xCFLAGS="$pgiopt -fast -Minfo=loop -pc 32"
    fi
    if test "$enable_fortran" = "yes"; then
      if $F77 -V 2>  /dev/null | grep Portland /dev/null 2>&1; then
	xFFLAGS="$xCFLAGS -Mneginfo=loop"
      fi	
    fi
    ;;
esac	
# Phew, end of all those operating systems and processors!			

# use default flags for gcc/g77 on all systems
if test $ac_cv_prog_gcc = yes; then
  xCFLAGS="-O6 -fomit-frame-pointer -finline-functions -funroll-all-loops -Wall -Wno-unused"
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
  # try to guess correct CPU flags, at least for linux
  case "${host_cpu}" in
    # i586/i686 cpu flags don't improve speed, thus no need to use them.
    # don't check f77 separately - we assume f77 and gcc are similar
	  
    powerpc*)
      cputype=`(grep cpu /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/ //g') 2> /dev/null`
      is60x=`echo $cputype | egrep "^60[0-9]e?$"`
      if test -n "$is60x"; then
	ACX_CHECK_CC_FLAGS(-mcpu=$cputype,m_cpu_60x,CPU_FLAGS=-mcpu=$cputype)
      elif test "$cputype" = 750; then
        ACX_CHECK_CC_FLAGS(-mcpu=750,m_cpu_750,CPU_FLAGS=-mcpu=750)
      fi
      if test -z "$CPU_FLAGS"; then
        ACX_CHECK_CC_FLAGS(-mcpu=powerpc,m_cpu_powerpc,CPU_FLAGS=-mcpu=powerpc)
      fi
      if test -z "$CPU_FLAGS"; then
	ACX_CHECK_CC_FLAGS(-mpowerpc,m_powerpc,CPU_FLAGS=-mpowerpc)
      fi
   esac
fi

if test -n "$CPU_FLAGS"; then
  xCFLAGS="$xCFLAGS $CPU_FLAGS"
  xFFLAGS="$xFFLAGS $CPU_FLAGS"
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
    echo "*******************************************************************"
    echo "* WARNING: No special optimization settings found for the C       *"
    echo "* compiler. Use  make CFLAGS=..., or edit the top level Makefile. *"
    echo "* Reverting to the default setting CFLAGS=-O3. (mail us about it)  *"
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
      echo "*******************************************************************"
      echo "* WARNING: No special optimization settings found for the fortran *"
      echo "* compiler. Use  make FFLAGS=..., or edit the top level Makefile. *"
      echo "* Reverting to the default setting FFLAGS=-O3. (mail us about it) *"
      echo "*******************************************************************"
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
])






