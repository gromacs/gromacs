#
# RPM specification file to make normal gromacs packages (without mpi)
# If you have mpi installed, you can create an mpi mdrun executable
# and libs with the gromacs-mpi spec file!

#
# Main package - only dynamic libs, and no header files
#
Summary: Molecular dynamics package (non-parallel)
Name: gromacs
Version: 3.0
Release: 1
Copyright: GPL
Group: Applications/Science
Prefix: /usr/local/gromacs
Requires: fftw >= 2.1.3 
Source: ftp://ftp.gromacs.org/pub/gromacs/source/gromacs-%{version}.tar.gz
URL: http://www.gromacs.org
Packager: Erik Lindahl <lindahl@gromacs.org>
%description
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high 
performance means it is used also in several other field 
like polymer chemistry and solid state physics. This
version has the dynamic libs and executables; to hack new
utility programs you also need the headers and static
libs in gromacs-dev. Linux kernel 2.4 or later is STRONGLY
recommended on Pentium III and later processors since
GROMACS then can use assembly loops with SSE instructions.

#
# The header files and static libraries go into gromacs-devel...
#
%package devel
Summary: Header files and static libraries for GROMACS
Group: Applications/Science
Prefix: %{prefix}
Requires: fftw-devel >=2.1.3, gromacs = %{version}-%{release}
%description devel
This package contains header files, static libraries,
and a program example for the GROMACS molecular
dynamics software. You need it if you want to write your
own analysis programs. A word of warning, though; 
this package is still somewhat untidy, and might put a 
lot of files in your standard include dir if you change
the default prefix!

%prep
%setup

%build
# We can actually install in /usr/local/bin, /usr/local/lib, etc,  by using the options
# --prefix=/usr/local, --exec-prefix=/usr/local  and --datadir=/usr/local/share/gromacs, 
# but since the development package puts a lot of include files in {prefix}/include we
# default to /usr/local/gromacs. Will try to fix that in gromacs 4.0 :-)

./configure
make 

%install
make DESTDIR=${RPM_BUILD_ROOT} install
# Move mdrun to mdrun_nompi - we make a link at post-install script time!
(cd ${RPM_BUILD_ROOT}%{prefix}/%{_target}/bin && mv mdrun mdrun_nompi)

%clean
rm -rf ${RPM_BUILD_ROOT}

%post
#
# Add our (final) library directory to /etc/ld.so.conf if it is not already there
#
if test -z `grep ${RPM_INSTALL_PREFIX}/%{_target}/lib /etc/ld.so.conf`; then
     cat >> /etc/ld.so.conf < ${RPM_INSTALL_PREFIX}/%{_target}/lib
fi
# Make a link from mdrun_nompi to mdrun - if it doesn't already exist!
(cd ${RPM_INSTALL_PREFIX}/%{_target}/bin && test ! -e mdrun && ln -s mdrun_nompi mdrun)

# run ldconfig to update the runtime linker database with the new libraries
# (make sure /sbin is in the $PATH)
PATH="/sbin:$PATH" ldconfig

%postun
# after uninstall, run ldconfig to remove the libs from the linker database
PATH="/sbin:$PATH" ldconfig
# and remove the link from mdrun_nompi to mdrun
(cd ${RPM_INSTALL_PREFIX}/%{_target}/bin && rm -f mdrun)

%files
%defattr(-,root,root)
%doc /usr/local/gromacs/README
/usr/local/gromacs/%{_target}/bin/*
/usr/local/gromacs/share/top/*
/usr/local/gromacs/share/tutor/*
%docdir /usr/local/gromacs/share/html
/usr/local/share/gromacs/html/
%docdir /usr/local/gromacs/man
/usr/local/gromacs/man/*
/usr/local/gromacs/%{_target}/lib/libgmx.so.1.0.0
/usr/local/gromacs/%{_target}/lib/libgmx.so.1
/usr/local/gromacs/%{_target}/lib/libmd.so.1.0.0
/usr/local/gromacs/%{_target}/lib/libmd.so.1
%files devel
%defattr(-,root,root)
/usr/local/gromacs/share/template/*
/usr/local/gromacs/%{_target}/lib/libgmx.so
/usr/local/gromacs/%{_target}/lib/libgmx.a
/usr/local/gromacs/%{_target}/lib/libgmx.la
/usr/local/gromacs/%{_target}/lib/libmd.so
/usr/local/gromacs/%{_target}/lib/libmd.a
/usr/local/gromacs/%{_target}/lib/libmd.la
/usr/local/gromacs/include/*
