 
# ACX_CHECK_FFTW()
# ----------------
# This macro checks for fftw header files and libraries,
# including the possible prefixing with s or d to determine precision.
# Arg 1 is the fftw header/library name to check for, without 
# prefix or anything else (e.g. rfftw_mpi for real MPI transforms)
# Arg 2 is the value of $enable_float, i.e. yes/no for precision=4/8.
AC_DEFUN(ACX_CHECK_FFTW,
[
if test -z "$ac_fftw_firstname"; then

if test "$2" = "no"; then
  prec="double"
  realsize=8
  fftwcheckprefix=d
else
  prec="single"
  realsize=4
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
# we cannot run the code since an MPI program might not be allowed on a login node of a supercomputer
AC_TRY_COMPILE([#include <$fftwname.h>],
[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $realsize)]; ],
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
[Do you have $prec precision FFTW installed? If you are using packages,]
[note that you also need fftw-devel to compile GROMACS. You can find the ]
[software at www.fftw.org, and detailed instructions at www.gromacs.org.]
[If you compiled FFTW yourself:                                        ]
[Note that the default FFTW setup is double precision. Change the FFTW]
[configuration to single with --enable-float. If you want MPI support,]
[use --enable-mpi. It is a good idea to install both single & double.] 
[If your sysadm doesn't want to install it you can do it to a location]
[in your home directory and provide the correct paths in the CPPFLAGS]
[and LDFLAGS environment variables before running configure.]
[That is also necessary to do if your compiler doesn't search]
[/usr/local/include and /usr/local/lib by default.]
[You can find information at www.gromacs.org, or in the INSTALL file.])
])
AC_TRY_COMPILE([#include <$xfftwname.h>],[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $realsize)];],
[
fftwname=$xfftwname 
usedprefix=$fftwcheckprefix
],
[
AC_MSG_ERROR([Cannot find any $prec precision $fftwname.h or $xfftwname.h]
[Do you have $prec precision FFTW installed? If you are using packages,]
[note that you also need fftw-devel to compile GROMACS. You can find the ]
[software at www.fftw.org, and detailed instructions at www.gromacs.org.]
[If you compiled FFTW yourself:                                       ]
[Note that the default FFTW setup is double precision. Change the FFTW]
[configuration to single with --enable-float. If you want MPI support,]
[use --enable-mpi. It is a good idea to install both single & double.] 
[If your sysadm doesn't want to install it you can do it to a location]
[in your home directory and provide the correct paths in the CPPFLAGS]
[and LDFLAGS environment variables before running configure.]
[That is also necessary to do if your compiler doesn't search]
[/usr/local/include and /usr/local/lib by default.]
[You can find information at www.gromacs.org, or in the INSTALL file.])])
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
    LDFLAGS="$LDFLAGS -woff 84"

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
      cc_vendor="Compaq"
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
      if $F77 -V 2>  /dev/null | grep Portland > /dev/null 2>&1; then
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
  ASFLAGS="$ASFLAGS -x assembler-with-cpp"
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
else
  AM_CONDITIONAL(GNU_CC,false)
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
    echo "* Reverting to the default setting CFLAGS=-O3. (mail us about it!)*"
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
      echo "* Reverting to the default setting FFLAGS=-O3. (mail us about it!)*"
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






## libtool.m4 - Configure libtool for the host system. -*-Shell-script-*-
## Copyright 1996, 1997, 1998, 1999, 2000, 2001
## Free Software Foundation, Inc.
## Originally by Gordon Matzigkeit <gord@gnu.ai.mit.edu>, 1996
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## As a special exception to the GNU General Public License, if you
## distribute this file as part of a program that contains a
## configuration script generated by Autoconf, you may include it under
## the same distribution terms that you use for the rest of that program.

# serial 46 AC_PROG_LIBTOOL
AC_DEFUN([AC_PROG_LIBTOOL],
[AC_REQUIRE([_AC_PROG_LIBTOOL])dnl
dnl If AC_PROG_CXX has already been expanded, run AC_LIBTOOL_CXX
dnl immediately, otherwise, hook it in at the end of AC_PROG_CXX.
  AC_PROVIDE_IFELSE([AC_PROG_CXX],
    [AC_LIBTOOL_CXX],
    [define([AC_PROG_CXX], defn([AC_PROG_CXX])[AC_LIBTOOL_CXX
])])
  AC_PROVIDE_IFELSE([AC_PROG_F77],
    [AC_LIBTOOL_F77],
    [define([AC_PROG_F77], defn([AC_PROG_F77])[AC_LIBTOOL_F77
])])

dnl Quote A][M_PROG_GCJ so that aclocal doesn't bring it in needlessly.
dnl If either AC_PROG_GCJ or A][M_PROG_GCJ have already been expanded, run
dnl AC_LIBTOOL_GCJ immediately, otherwise, hook it in at the end of both.
  AC_PROVIDE_IFELSE([AC_PROG_GCJ],
    [AC_LIBTOOL_GCJ],
    [AC_PROVIDE_IFELSE([A][M_PROG_GCJ],
        [AC_LIBTOOL_GCJ],
	[AC_PROVIDE_IFELSE([LT_AC_PROG_GCJ],
	  [AC_LIBTOOL_GCJ],
	[ifdef([AC_PROG_GCJ],
	       [define([AC_PROG_GCJ], defn([AC_PROG_GCJ])[AC_LIBTOOL_GCJ
])])
	 ifdef([A][M_PROG_GCJ],
	       [define([A][M_PROG_GCJ], defn([A][M_PROG_GCJ])[AC_LIBTOOL_GCJ
])])
	 ifdef([LT_AC_PROG_GCJ],
	       [define([LT_AC_PROG_GCJ], defn([LT_AC_PROG_GCJ])[AC_LIBTOOL_GCJ
])])])])])])

AC_DEFUN([_AC_PROG_LIBTOOL],
[AC_REQUIRE([AC_LIBTOOL_SETUP])dnl
AC_BEFORE([$0],[AC_LIBTOOL_CXX])dnl
AC_BEFORE([$0],[AC_LIBTOOL_GCJ])dnl

# Save cache, so that ltconfig can load it
AC_CACHE_SAVE

# Actually configure libtool.  ac_aux_dir is where install-sh is found.
AR="$AR" LTCC="$CC" CC="$CC" CFLAGS="$CFLAGS" CPPFLAGS="$CPPFLAGS" \
MAGIC_CMD="$MAGIC_CMD" LD="$LD" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB" STRIP="$STRIP" \
AS="$AS" DLLTOOL="$DLLTOOL" OBJDUMP="$OBJDUMP" \
objext="$OBJEXT" exeext="$EXEEXT" reload_flag="$reload_flag" \
deplibs_check_method="$deplibs_check_method" file_magic_cmd="$file_magic_cmd" \
${CONFIG_SHELL-/bin/sh} $ac_aux_dir/ltconfig --no-reexec \
$libtool_flags --no-verify --build="$build" $ac_aux_dir/ltmain.sh $host \
|| AC_MSG_ERROR([libtool configure failed])

# Reload cache, that may have been modified by ltconfig
AC_CACHE_LOAD

# This can be used to rebuild libtool when needed
LIBTOOL_DEPS="$ac_aux_dir/ltconfig $ac_aux_dir/ltmain.sh $ac_aux_dir/ltcf-c.sh"

# Always use our own libtool.
LIBTOOL='$(SHELL) $(top_builddir)/libtool'
AC_SUBST(LIBTOOL)dnl

# Redirect the config.log output again, so that the ltconfig log is not
# clobbered by the next message.
exec 5>>./config.log
])

AC_DEFUN([AC_LIBTOOL_SETUP],
[AC_PREREQ(2.13)dnl
AC_REQUIRE([AC_ENABLE_SHARED])dnl
AC_REQUIRE([AC_ENABLE_STATIC])dnl
AC_REQUIRE([AC_ENABLE_FAST_INSTALL])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_LD])dnl
AC_REQUIRE([AC_PROG_LD_RELOAD_FLAG])dnl
AC_REQUIRE([AC_PROG_NM])dnl
AC_REQUIRE([AC_PROG_LN_S])dnl
AC_REQUIRE([AC_DEPLIBS_CHECK_METHOD])dnl
# Autoconf 2.13's AC_OBJEXT and AC_EXEEXT macros only works for C compilers!
AC_REQUIRE([AC_OBJEXT])dnl
AC_REQUIRE([AC_EXEEXT])dnl
dnl

# Only perform the check for file, if the check method requires it
case $deplibs_check_method in
file_magic*)
  if test "$file_magic_cmd" = '$MAGIC_CMD'; then
    AC_PATH_MAGIC
  fi
  ;;
esac

AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_CHECK_TOOL(STRIP, strip, :)

# Check for any special flags to pass to ltconfig.
libtool_flags="--cache-file=$cache_file"
test "$enable_shared" = no && libtool_flags="$libtool_flags --disable-shared"
test "$enable_static" = no && libtool_flags="$libtool_flags --disable-static"
test "$enable_fast_install" = no && libtool_flags="$libtool_flags --disable-fast-install"
test "$GCC" = yes && libtool_flags="$libtool_flags --with-gcc"
test "$lt_cv_prog_gnu_ld" = yes && libtool_flags="$libtool_flags --with-gnu-ld"
ifdef([AC_PROVIDE_AC_LIBTOOL_DLOPEN],
[libtool_flags="$libtool_flags --enable-dlopen"])
ifdef([AC_PROVIDE_AC_LIBTOOL_WIN32_DLL],
[libtool_flags="$libtool_flags --enable-win32-dll"])
AC_ARG_ENABLE(libtool-lock,
  [  --disable-libtool-lock        avoid locking (might break parallel builds)])
test "x$enable_libtool_lock" = xno && libtool_flags="$libtool_flags --disable-lock"
test x"$silent" = xyes && libtool_flags="$libtool_flags --silent"

AC_ARG_WITH(pic,
  [  --with-pic                    try to use only PIC/non-PIC [default=both]],
     pic_mode="$withval", pic_mode=default)
test x"$pic_mode" = xyes && libtool_flags="$libtool_flags --prefer-pic"
test x"$pic_mode" = xno && libtool_flags="$libtool_flags --prefer-non-pic"

# Some flags need to be propagated to the compiler or linker for good
# libtool support.
case $host in
*-*-irix6*)
  # Find out which ABI we are using.
  echo '[#]line __oline__ "configure"' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.$ac_objext` in
    *32-bit*)
      LD="${LD-ld} -32"
      ;;
    *N32*)
      LD="${LD-ld} -n32"
      ;;
    *64-bit*)
      LD="${LD-ld} -64"
      ;;
    esac
  fi
  rm -rf conftest*
  ;;

*-*-sco3.2v5*)
  # On SCO OpenServer 5, we need -belf to get full-featured binaries.
  SAVE_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -belf"
  AC_CACHE_CHECK([whether the C compiler needs -belf], lt_cv_cc_needs_belf,
    [AC_LANG_SAVE
     AC_LANG_C
     AC_TRY_LINK([],[],[lt_cv_cc_needs_belf=yes],[lt_cv_cc_needs_belf=no])
     AC_LANG_RESTORE])
  if test x"$lt_cv_cc_needs_belf" != x"yes"; then
    # this is probably gcc 2.8.0, egcs 1.0 or newer; no need for -belf
    CFLAGS="$SAVE_CFLAGS"
  fi
  ;;

ifdef([AC_PROVIDE_AC_LIBTOOL_WIN32_DLL],
[*-*-cygwin* | *-*-mingw* | *-*-pw32*)
  AC_CHECK_TOOL(DLLTOOL, dlltool, false)
  AC_CHECK_TOOL(AS, as, false)
  AC_CHECK_TOOL(OBJDUMP, objdump, false)

  # recent cygwin and mingw systems supply a stub DllMain which the user
  # can override, but on older systems we have to supply one
  AC_CACHE_CHECK([if libtool should supply DllMain function], lt_cv_need_dllmain,
    [AC_TRY_LINK([],
      [extern int __attribute__((__stdcall__)) DllMain(void*, int, void*);
      DllMain (0, 0, 0);],
      [lt_cv_need_dllmain=no],[lt_cv_need_dllmain=yes])])

  case $host/$CC in
  *-*-cygwin*/gcc*-mno-cygwin*|*-*-mingw*)
    # old mingw systems require "-dll" to link a DLL, while more recent ones
    # require "-mdll"
    SAVE_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS -mdll"
    AC_CACHE_CHECK([how to link DLLs], lt_cv_cc_dll_switch,
      [AC_TRY_LINK([], [], [lt_cv_cc_dll_switch=-mdll],[lt_cv_cc_dll_switch=-dll])])
    CFLAGS="$SAVE_CFLAGS" ;;
  *-*-cygwin* | *-*-pw32*)
    # cygwin systems need to pass --dll to the linker, and not link
    # crt.o which will require a WinMain@16 definition.
    lt_cv_cc_dll_switch="-Wl,--dll -nostartfiles" ;;
  esac
  ;;
  ])
esac
])

# AC_LIBTOOL_DLOPEN - enable checks for dlopen support
AC_DEFUN([AC_LIBTOOL_DLOPEN], [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])])

# AC_LIBTOOL_WIN32_DLL - declare package support for building win32 dll's
AC_DEFUN([AC_LIBTOOL_WIN32_DLL], [AC_BEFORE([$0], [AC_LIBTOOL_SETUP])])

# AC_ENABLE_SHARED - implement the --enable-shared flag
# Usage: AC_ENABLE_SHARED[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN([AC_ENABLE_SHARED],
[define([AC_ENABLE_SHARED_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(shared,
changequote(<<, >>)dnl
<<  --enable-shared[=PKGS]        build shared libraries [default=>>AC_ENABLE_SHARED_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case $enableval in
yes) enable_shared=yes ;;
no) enable_shared=no ;;
*)
  enable_shared=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_shared=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_shared=AC_ENABLE_SHARED_DEFAULT)dnl
])

# AC_DISABLE_SHARED - set the default shared flag to --disable-shared
AC_DEFUN([AC_DISABLE_SHARED], [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_SHARED(no)])

# AC_ENABLE_STATIC - implement the --enable-static flag
# Usage: AC_ENABLE_STATIC[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN([AC_ENABLE_STATIC],
[define([AC_ENABLE_STATIC_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(static,
changequote(<<, >>)dnl
<<  --enable-static[=PKGS]        build static libraries [default=>>AC_ENABLE_STATIC_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case $enableval in
yes) enable_static=yes ;;
no) enable_static=no ;;
*)
  enable_static=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_static=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_static=AC_ENABLE_STATIC_DEFAULT)dnl
])

# AC_DISABLE_STATIC - set the default static flag to --disable-static
AC_DEFUN([AC_DISABLE_STATIC],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_STATIC(no)])


# AC_ENABLE_FAST_INSTALL - implement the --enable-fast-install flag
# Usage: AC_ENABLE_FAST_INSTALL[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN([AC_ENABLE_FAST_INSTALL],
[define([AC_ENABLE_FAST_INSTALL_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(fast-install,
changequote(<<, >>)dnl
<<  --enable-fast-install[=PKGS]  optimize for fast installation [default=>>AC_ENABLE_FAST_INSTALL_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case $enableval in
yes) enable_fast_install=yes ;;
no) enable_fast_install=no ;;
*)
  enable_fast_install=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_fast_install=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_fast_install=AC_ENABLE_FAST_INSTALL_DEFAULT)dnl
])

# AC_DISABLE_FAST_INSTALL - set the default to --disable-fast-install
AC_DEFUN([AC_DISABLE_FAST_INSTALL],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_FAST_INSTALL(no)])

# AC_LIBTOOL_PICMODE - implement the --with-pic flag
# Usage: AC_LIBTOOL_PICMODE[(MODE)]
#   Where MODE is either `yes' or `no'.  If omitted, it defaults to
#   `both'.
AC_DEFUN([AC_LIBTOOL_PICMODE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
pic_mode=ifelse($#,1,$1,default)])


# AC_PATH_TOOL_PREFIX - find a file program which can recognise shared library
AC_DEFUN([AC_PATH_TOOL_PREFIX],
[AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(lt_cv_path_MAGIC_CMD,
[case $MAGIC_CMD in
  /*)
  lt_cv_path_MAGIC_CMD="$MAGIC_CMD" # Let the user override the test with a path.
  ;;
  ?:/*)
  lt_cv_path_MAGIC_CMD="$MAGIC_CMD" # Let the user override the test with a dos path.
  ;;
  *)
  ac_save_MAGIC_CMD="$MAGIC_CMD"
  IFS="${IFS=   }"; ac_save_ifs="$IFS"; IFS=":"
dnl $ac_dummy forces splitting on constant user-supplied paths.
dnl POSIX.2 word splitting is done only on the output of word expansions,
dnl not every word.  This closes a longstanding sh security hole.
  ac_dummy="ifelse([$2], , $PATH, [$2])"
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$1; then
      lt_cv_path_MAGIC_CMD="$ac_dir/$1"
      if test -n "$file_magic_test_file"; then
	case $deplibs_check_method in
	"file_magic "*)
	  file_magic_regex="`expr \"$deplibs_check_method\" : \"file_magic \(.*\)\"`"
	  MAGIC_CMD="$lt_cv_path_MAGIC_CMD"
	  if eval $file_magic_cmd \$file_magic_test_file 2> /dev/null |
	    egrep "$file_magic_regex" > /dev/null; then
	    :
	  else
	    cat <<EOF 1>&2

*** Warning: the command libtool uses to detect shared libraries,
*** $file_magic_cmd, produces output that libtool cannot recognize.
*** The result is that libtool may fail to recognize shared libraries
*** as such.  This will affect the creation of libtool libraries that
*** depend on shared libraries, but programs linked with such libtool
*** libraries will work regardless of this problem.  Nevertheless, you
*** may want to report the problem to your system manager and/or to
*** bug-libtool@gnu.org

EOF
	  fi ;;
	esac
      fi
      break
    fi
  done
  IFS="$ac_save_ifs"
  MAGIC_CMD="$ac_save_MAGIC_CMD"
  ;;
esac])
MAGIC_CMD="$lt_cv_path_MAGIC_CMD"
if test -n "$MAGIC_CMD"; then
  AC_MSG_RESULT($MAGIC_CMD)
else
  AC_MSG_RESULT(no)
fi
])


# AC_PATH_MAGIC - find a file program which can recognise a shared library
AC_DEFUN([AC_PATH_MAGIC],
[AC_REQUIRE([AC_CHECK_TOOL_PREFIX])dnl
AC_PATH_TOOL_PREFIX(${ac_tool_prefix}file, /usr/bin:$PATH)
if test -z "$lt_cv_path_MAGIC_CMD"; then
  if test -n "$ac_tool_prefix"; then
    AC_PATH_TOOL_PREFIX(file, /usr/bin:$PATH)
  else
    MAGIC_CMD=:
  fi
fi
])


# AC_PROG_LD - find the path to the GNU or non-GNU linker
AC_DEFUN([AC_PROG_LD],
[AC_ARG_WITH(gnu-ld,
[  --with-gnu-ld                 assume the C compiler uses GNU ld [default=no]],
test "$withval" = no || with_gnu_ld=yes, with_gnu_ld=no)
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
ac_prog=ld
if test "$GCC" = yes; then
  # Check if gcc -print-prog-name=ld gives a path.
  AC_MSG_CHECKING([for ld used by GCC])
  case $host in
  *-*-mingw*)
    # gcc leaves a trailing carriage return which upsets mingw
    ac_prog=`($CC -print-prog-name=ld) 2>&5 | tr -d '\015'` ;;
  *)
    ac_prog=`($CC -print-prog-name=ld) 2>&5` ;;
  esac
  case $ac_prog in
    # Accept absolute paths.
    [[\\/]* | [A-Za-z]:[\\/]*)]
      re_direlt=['/[^/][^/]*/\.\./']
      # Canonicalize the path of ld
      ac_prog=`echo $ac_prog| sed 's%\\\\%/%g'`
      while echo $ac_prog | grep "$re_direlt" > /dev/null 2>&1; do
	ac_prog=`echo $ac_prog| sed "s%$re_direlt%/%"`
      done
      test -z "$LD" && LD="$ac_prog"
      ;;
  "")
    # If it fails, then pretend we aren't using GCC.
    ac_prog=ld
    ;;
  *)
    # If it is relative, then search for the first ld in PATH.
    with_gnu_ld=unknown
    ;;
  esac
elif test "$with_gnu_ld" = yes; then
  AC_MSG_CHECKING([for GNU ld])
else
  AC_MSG_CHECKING([for non-GNU ld])
fi
AC_CACHE_VAL(lt_cv_path_LD,
[if test -z "$LD"; then
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}${PATH_SEPARATOR-:}"
  for ac_dir in $PATH; do
    test -z "$ac_dir" && ac_dir=.
    if test -f "$ac_dir/$ac_prog" || test -f "$ac_dir/$ac_prog$ac_exeext"; then
      lt_cv_path_LD="$ac_dir/$ac_prog"
      # Check to see if the program is GNU ld.  I'd rather use --version,
      # but apparently some GNU ld's only accept -v.
      # Break only if it was the GNU/non-GNU ld that we prefer.
      if "$lt_cv_path_LD" -v 2>&1 < /dev/null | egrep '(GNU|with BFD)' > /dev/null; then
	test "$with_gnu_ld" != no && break
      else
	test "$with_gnu_ld" != yes && break
      fi
    fi
  done
  IFS="$ac_save_ifs"
else
  lt_cv_path_LD="$LD" # Let the user override the test with a path.
fi])
LD="$lt_cv_path_LD"
if test -n "$LD"; then
  AC_MSG_RESULT($LD)
else
  AC_MSG_RESULT(no)
fi
test -z "$LD" && AC_MSG_ERROR([no acceptable ld found in \$PATH])
AC_PROG_LD_GNU
])

AC_DEFUN([AC_PROG_LD_GNU],
[AC_CACHE_CHECK([if the linker ($LD) is GNU ld], lt_cv_prog_gnu_ld,
[# I'd rather use --version here, but apparently some GNU ld's only accept -v.
if $LD -v 2>&1 </dev/null | egrep '(GNU|with BFD)' 1>&5; then
  lt_cv_prog_gnu_ld=yes
else
  lt_cv_prog_gnu_ld=no
fi])
with_gnu_ld=$lt_cv_prog_gnu_ld
])

# AC_PROG_LD_RELOAD_FLAG - find reload flag for linker
#   -- PORTME Some linkers may need a different reload flag.
AC_DEFUN([AC_PROG_LD_RELOAD_FLAG],
[AC_CACHE_CHECK([for $LD option to reload object files], lt_cv_ld_reload_flag,
[lt_cv_ld_reload_flag='-r'])
reload_flag=$lt_cv_ld_reload_flag
test -n "$reload_flag" && reload_flag=" $reload_flag"
])

# AC_DEPLIBS_CHECK_METHOD - how to check for library dependencies
#  -- PORTME fill in with the dynamic library characteristics
AC_DEFUN([AC_DEPLIBS_CHECK_METHOD],
[AC_CACHE_CHECK([how to recognise dependant libraries],
lt_cv_deplibs_check_method,
[lt_cv_file_magic_cmd='$MAGIC_CMD'
lt_cv_file_magic_test_file=
lt_cv_deplibs_check_method='unknown'
# Need to set the preceding variable on all platforms that support
# interlibrary dependencies.
# 'none' -- dependencies not supported.
# `unknown' -- same as none, but documents that we really don't know.
# 'pass_all' -- all dependencies passed with no checks.
# 'test_compile' -- check by making test program.
# 'file_magic [regex]' -- check by looking for files in library path
# which responds to the $file_magic_cmd with a given egrep regex.
# If you have `file' or equivalent on your system and you're not sure
# whether `pass_all' will *always* work, you probably want this one.

case $host_os in
aix*)
  lt_cv_deplibs_check_method=pass_all
  ;;

beos*)
  lt_cv_deplibs_check_method=pass_all
  ;;

bsdi4*)
  lt_cv_deplibs_check_method=['file_magic ELF [0-9][0-9]*-bit [ML]SB (shared object|dynamic lib)']
  lt_cv_file_magic_cmd='/usr/bin/file -L'
  lt_cv_file_magic_test_file=/shlib/libc.so
  ;;

cygwin* | mingw* |pw32*)
  lt_cv_deplibs_check_method='file_magic file format pei*-i386(.*architecture: i386)?'
  lt_cv_file_magic_cmd='$OBJDUMP -f'
  ;;

darwin* | rhapsody*)
  lt_cv_deplibs_check_method='file_magic Mach-O dynamically linked shared library'
  lt_cv_file_magic_cmd='/usr/bin/file -L'
  case "$host_os" in
  rhapsody* | darwin1.[012])
    lt_cv_file_magic_test_file='/System/Library/Frameworks/System.framework/System'
    ;;
  *) # Darwin 1.3 on
    lt_cv_file_magic_test_file='/usr/lib/libSystem.dylib'
    ;;
  esac
  ;;

freebsd* )
  if echo __ELF__ | $CC -E - | grep __ELF__ > /dev/null; then
    case $host_cpu in
    i*86 )
      # Not sure whether the presence of OpenBSD here was a mistake.
      # Let's accept both of them until this is cleared up.
      lt_cv_deplibs_check_method=['file_magic (FreeBSD|OpenBSD)/i[3-9]86 (compact )?demand paged shared library']
      lt_cv_file_magic_cmd=/usr/bin/file
      lt_cv_file_magic_test_file=`echo /usr/lib/libc.so.*`
      ;;
    esac
  else
    lt_cv_deplibs_check_method=pass_all
  fi
  ;;

gnu*)
  lt_cv_deplibs_check_method=pass_all
  ;;

hpux10.20*|hpux11*)
  lt_cv_deplibs_check_method=['file_magic (s[0-9][0-9][0-9]|PA-RISC[0-9].[0-9]) shared library']
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=/usr/lib/libc.sl
  ;;

irix5* | irix6*)
  case $host_os in
  irix5*)
    # this will be overridden with pass_all, but let us keep it just in case
    lt_cv_deplibs_check_method="file_magic ELF 32-bit MSB dynamic lib MIPS - version 1"
    ;;
  *)
    case $LD in
    *-32|*"-32 ") libmagic=32-bit;;
    *-n32|*"-n32 ") libmagic=N32;;
    *-64|*"-64 ") libmagic=64-bit;;
    *) libmagic=never-match;;
    esac
    # this will be overridden with pass_all, but let us keep it just in case
    lt_cv_deplibs_check_method=["file_magic ELF ${libmagic} MSB mips-[1234] dynamic lib MIPS - version 1"]
    ;;
  esac
  lt_cv_file_magic_test_file=`echo /lib${libsuff}/libc.so*`
  lt_cv_deplibs_check_method=pass_all
  ;;

# This must be Linux ELF.
linux-gnu*)
  case $host_cpu in
  alpha* | hppa* | i*86 | powerpc* | sparc* | ia64* )
    lt_cv_deplibs_check_method=pass_all ;;
  *)
    # glibc up to 2.1.1 does not perform some relocations on ARM
    lt_cv_deplibs_check_method=['file_magic ELF [0-9][0-9]*-bit [LM]SB (shared object|dynamic lib )'] ;;
  esac
  lt_cv_file_magic_test_file=`echo /lib/libc.so* /lib/libc-*.so`
  ;;

netbsd*)
  if echo __ELF__ | $CC -E - | grep __ELF__ > /dev/null; then
    [lt_cv_deplibs_check_method='match_pattern /lib[^/\.]+\.so\.[0-9]+\.[0-9]+$']
  else
    [lt_cv_deplibs_check_method='match_pattern /lib[^/\.]+\.so$']
  fi
  ;;

newsos6)
  [lt_cv_deplibs_check_method='file_magic ELF [0-9][0-9]*-bit [ML]SB (executable|dynamic lib)']
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=/usr/lib/libnls.so
  ;;

osf3* | osf4* | osf5*)
  # this will be overridden with pass_all, but let us keep it just in case
  lt_cv_deplibs_check_method='file_magic COFF format alpha shared library'
  lt_cv_file_magic_test_file=/shlib/libc.so
  lt_cv_deplibs_check_method=pass_all
  ;;

sco3.2v5*)
  lt_cv_deplibs_check_method=pass_all
  ;;

solaris*)
  lt_cv_deplibs_check_method=pass_all
  lt_cv_file_magic_test_file=/lib/libc.so
  ;;

[sysv5uw[78]* | sysv4*uw2*)]
  lt_cv_deplibs_check_method=pass_all
  ;;

sysv4 | sysv4.2uw2* | sysv4.3* | sysv5*)
  case $host_vendor in
  ncr)
    lt_cv_deplibs_check_method=pass_all
    ;;
  motorola)
    lt_cv_deplibs_check_method=['file_magic ELF [0-9][0-9]*-bit [ML]SB (shared object|dynamic lib) M[0-9][0-9]* Version [0-9]']
    lt_cv_file_magic_test_file=`echo /usr/lib/libc.so*`
    ;;
  esac
  ;;
esac
])
file_magic_cmd=$lt_cv_file_magic_cmd
deplibs_check_method=$lt_cv_deplibs_check_method
])


# AC_PROG_NM - find the path to a BSD-compatible name lister
AC_DEFUN([AC_PROG_NM],
[AC_MSG_CHECKING([for BSD-compatible nm])
AC_CACHE_VAL(lt_cv_path_NM,
[if test -n "$NM"; then
  # Let the user override the test.
  lt_cv_path_NM="$NM"
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}${PATH_SEPARATOR-:}"
  for ac_dir in $PATH /usr/ccs/bin /usr/ucb /bin; do
    test -z "$ac_dir" && ac_dir=.
    tmp_nm=$ac_dir/${ac_tool_prefix}nm
    if test -f $tmp_nm || test -f $tmp_nm$ac_exeext ; then
      # Check to see if the nm accepts a BSD-compat flag.
      # Adding the `sed 1q' prevents false positives on HP-UX, which says:
      #   nm: unknown option "B" ignored
      # Tru64's nm complains that /dev/null is an invalid object file
      if ($tmp_nm -B /dev/null 2>&1 | sed '1q'; exit 0) | egrep '(/dev/null|Invalid file or object type)' >/dev/null; then
	lt_cv_path_NM="$tmp_nm -B"
	break
      elif ($tmp_nm -p /dev/null 2>&1 | sed '1q'; exit 0) | egrep /dev/null >/dev/null; then
	lt_cv_path_NM="$tmp_nm -p"
	break
      else
	lt_cv_path_NM=${lt_cv_path_NM="$tmp_nm"} # keep the first match, but
	continue # so that we can try to find one that supports BSD flags
      fi
    fi
  done
  IFS="$ac_save_ifs"
  test -z "$lt_cv_path_NM" && lt_cv_path_NM=nm
fi])
NM="$lt_cv_path_NM"
AC_MSG_RESULT([$NM])
])

# AC_CHECK_LIBM - check for math library
AC_DEFUN([AC_CHECK_LIBM],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
LIBM=
case $host in
*-*-beos* | *-*-cygwin* | *-*-pw32*)
  # These system don't have libm
  ;;
*-ncr-sysv4.3*)
  AC_CHECK_LIB(mw, _mwvalidcheckl, LIBM="-lmw")
  AC_CHECK_LIB(m, main, LIBM="$LIBM -lm")
  ;;
*)
  AC_CHECK_LIB(m, main, LIBM="-lm")
  ;;
esac
])

# AC_LIBLTDL_CONVENIENCE[(dir)] - sets LIBLTDL to the link flags for
# the libltdl convenience library and INCLTDL to the include flags for
# the libltdl header and adds --enable-ltdl-convenience to the
# configure arguments.  Note that LIBLTDL and INCLTDL are not
# AC_SUBSTed, nor is AC_CONFIG_SUBDIRS called.  If DIR is not
# provided, it is assumed to be `libltdl'.  LIBLTDL will be prefixed
# with '${top_builddir}/' and INCLTDL will be prefixed with
# '${top_srcdir}/' (note the single quotes!).  If your package is not
# flat and you're not using automake, define top_builddir and
# top_srcdir appropriately in the Makefiles.
AC_DEFUN([AC_LIBLTDL_CONVENIENCE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  case $enable_ltdl_convenience in
  no) AC_MSG_ERROR([this package needs a convenience libltdl]) ;;
  "") enable_ltdl_convenience=yes
      ac_configure_args="$ac_configure_args --enable-ltdl-convenience" ;;
  esac
  LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdlc.la
  INCLTDL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
])

# AC_LIBLTDL_INSTALLABLE[(dir)] - sets LIBLTDL to the link flags for
# the libltdl installable library and INCLTDL to the include flags for
# the libltdl header and adds --enable-ltdl-install to the configure
# arguments.  Note that LIBLTDL and INCLTDL are not AC_SUBSTed, nor is
# AC_CONFIG_SUBDIRS called.  If DIR is not provided and an installed
# libltdl is not found, it is assumed to be `libltdl'.  LIBLTDL will
# be prefixed with '${top_builddir}/' and INCLTDL will be prefixed
# with '${top_srcdir}/' (note the single quotes!).  If your package is
# not flat and you're not using automake, define top_builddir and
# top_srcdir appropriately in the Makefiles.
# In the future, this macro may have to be called after AC_PROG_LIBTOOL.
AC_DEFUN([AC_LIBLTDL_INSTALLABLE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  AC_CHECK_LIB(ltdl, main,
  [test x"$enable_ltdl_install" != xyes && enable_ltdl_install=no],
  [if test x"$enable_ltdl_install" = xno; then
     AC_MSG_WARN([libltdl not installed, but installation disabled])
   else
     enable_ltdl_install=yes
   fi
  ])
  if test x"$enable_ltdl_install" = x"yes"; then
    ac_configure_args="$ac_configure_args --enable-ltdl-install"
    LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdl.la
    INCLTDL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
  else
    ac_configure_args="$ac_configure_args --enable-ltdl-install=no"
    LIBLTDL="-lltdl"
    INCLTDL=
  fi
])

# If this macro is not defined by Autoconf, define it here.
ifdef([AC_PROVIDE_IFELSE],
      [],
      [define([AC_PROVIDE_IFELSE],
              [ifdef([AC_PROVIDE_$1],
                     [$2], [$3])])])

# AC_LIBTOOL_F77 - enable support for fortran libraries
AC_DEFUN([AC_LIBTOOL_F77], [AC_REQUIRE([_AC_LIBTOOL_F77])])

AC_DEFUN([_AC_LIBTOOL_F77],
[
if test "$ac_cv_prog_f77_works" = "yes"; then
AC_REQUIRE([AC_PROG_F77])
LIBTOOL_DEPS=$LIBTOOL_DEPS" $ac_aux_dir/ltcf-f77.sh"
lt_save_CC="$CC"
lt_save_CFLAGS="$CFLAGS"
dnl Make sure LTCC is set to the C compiler, i.e. set LTCC before CC
dnl is set to the fortran compiler.
AR="$AR" LTCC="$CC" CC="$F77" F77="$F77" CFLAGS="$FFLAGS" CPPFLAGS="" \
MAGIC_CMD="$MAGIC_CMD" LD="$LD" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB" STRIP="$STRIP" \
AS="$AS" DLLTOOL="$DLLTOOL" OBJDUMP="$OBJDUMP" \
objext="$OBJEXT" exeext="$EXEEXT" reload_flag="$reload_flag" \
deplibs_check_method="$deplibs_check_method" \
file_magic_cmd="$file_magic_cmd" \
${CONFIG_SHELL-/bin/sh} $ac_aux_dir/ltconfig -o libtool $libtool_flags \
--build="$build" --add-tag=F77 $ac_aux_dir/ltcf-f77.sh $host \
|| AC_MSG_ERROR([libtool tag configuration failed])
CC="$lt_save_CC"
CFLAGS="$lt_save_CFLAGS"

# Redirect the config.log output again, so that the ltconfig log is not
# clobbered by the next message.
exec 5>>./config.log
fi
])

# AC_LIBTOOL_CXX - enable support for C++ libraries
AC_DEFUN([AC_LIBTOOL_CXX], [AC_REQUIRE([_AC_LIBTOOL_CXX])])

AC_DEFUN([_AC_LIBTOOL_CXX],
[AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AC_PROG_CXXCPP])
LIBTOOL_DEPS=$LIBTOOL_DEPS" $ac_aux_dir/ltcf-cxx.sh"
lt_save_CC="$CC"
lt_save_CFLAGS="$CFLAGS"
dnl Make sure LTCC is set to the C compiler, i.e. set LTCC before CC
dnl is set to the C++ compiler.
AR="$AR" LTCC="$CC" CC="$CXX" CXX="$CXX" CFLAGS="$CXXFLAGS" CPPFLAGS="$CPPFLAGS" \
MAGIC_CMD="$MAGIC_CMD" LD="$LD" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB" STRIP="$STRIP" \
AS="$AS" DLLTOOL="$DLLTOOL" OBJDUMP="$OBJDUMP" \
objext="$OBJEXT" exeext="$EXEEXT" reload_flag="$reload_flag" \
deplibs_check_method="$deplibs_check_method" \
file_magic_cmd="$file_magic_cmd" \
${CONFIG_SHELL-/bin/sh} $ac_aux_dir/ltconfig -o libtool $libtool_flags \
--build="$build" --add-tag=CXX $ac_aux_dir/ltcf-cxx.sh $host \
|| AC_MSG_ERROR([libtool tag configuration failed])
CC="$lt_save_CC"
CFLAGS="$lt_save_CFLAGS"

# Redirect the config.log output again, so that the ltconfig log is not
# clobbered by the next message.
exec 5>>./config.log
])

# AC_LIBTOOL_GCJ - enable support for GCJ libraries
AC_DEFUN([AC_LIBTOOL_GCJ],[AC_REQUIRE([_AC_LIBTOOL_GCJ])])

AC_DEFUN([_AC_LIBTOOL_GCJ],
[AC_REQUIRE([AC_PROG_LIBTOOL])
AC_PROVIDE_IFELSE([AC_PROG_GCJ],[],
  [AC_PROVIDE_IFELSE([A][M_PROG_GCJ],[],
    [AC_PROVIDE_IFELSE([LT_AC_PROG_GCJ],[],
      [ifdef([AC_PROG_GCJ],[AC_REQUIRE([AC_PROG_GCJ])],
         [ifdef([A][M_PROG_GCJ],[AC_REQUIRE([A][M_PROG_GCJ])],
           [AC_REQUIRE([A][C_PROG_GCJ_OR_A][M_PROG_GCJ])])])])])])
LIBTOOL_DEPS=$LIBTOOL_DEPS" $ac_aux_dir/ltcf-gcj.sh"
lt_save_CC="$CC"
lt_save_CFLAGS="$CFLAGS"
dnl Make sure LTCC is set to the C compiler, i.e. set LTCC before CC
dnl is set to the C++ compiler.
AR="$AR" LTCC="$CC" CC="$GCJ" CFLAGS="$GCJFLAGS" CPPFLAGS="$CPPFLAGS" \
MAGIC_CMD="$MAGIC_CMD" LD="$LD" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB" STRIP="$STRIP" \
AS="$AS" DLLTOOL="$DLLTOOL" OBJDUMP="$OBJDUMP" \
objext="$OBJEXT" exeext="$EXEEXT" reload_flag="$reload_flag" \
deplibs_check_method="$deplibs_check_method" \
file_magic_cmd="$file_magic_cmd" \
${CONFIG_SHELL-/bin/sh} $ac_aux_dir/ltconfig -o libtool $libtool_flags \
--build="$build" --add-tag=GCJ $ac_aux_dir/ltcf-gcj.sh $host \
|| AC_MSG_ERROR([libtool tag configuration failed])
CC="$lt_save_CC"
CFLAGS="$lt_save_CFLAGS"

# Redirect the config.log output again, so that the ltconfig log is not
# clobbered by the next message.
exec 5>>./config.log
])

dnl old names
AC_DEFUN([AM_PROG_LIBTOOL],   [AC_PROG_LIBTOOL])
AC_DEFUN([AM_ENABLE_SHARED],  [AC_ENABLE_SHARED($@)])
AC_DEFUN([AM_ENABLE_STATIC],  [AC_ENABLE_STATIC($@)])
AC_DEFUN([AM_DISABLE_SHARED], [AC_DISABLE_SHARED($@)])
AC_DEFUN([AM_DISABLE_STATIC], [AC_DISABLE_STATIC($@)])
AC_DEFUN([AM_PROG_LD],        [AC_PROG_LD])
AC_DEFUN([AM_PROG_NM],        [AC_PROG_NM])

dnl This is just to silence aclocal about the macro not being used
ifelse([AC_DISABLE_FAST_INSTALL])dnl
ifelse([AC_DISABLE_SHARED])dnl

AC_DEFUN([LT_AC_PROG_GCJ],
[AC_CHECK_TOOL(GCJ, gcj, no)
  test "x${GCJFLAGS+set}" = xset || GCJFLAGS="-g -O2"
  AC_SUBST(GCJFLAGS)
])
