#
# RPM specification file to make normal gromacs packages (without mpi)
# If you have mpi installed, you can create an mpi mdrun executable
# and libs with the gromacs-mpi spec file!

#
# Main package - only dynamic libs, and no header files
#
Summary: Molecular dynamics package (non-parallel version)
Name: gromacs
Version: 3.1.beta_20020215
Release: 1
Copyright: GPL
Group: Applications/Science
Prefix: /usr/local
Buildroot: %{_topdir}/buildroot
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
You can also perform parallel simulations if you install
gromacs-lammpi.

#
# The header files and static libraries go into gromacs-devel...
#
%package devel
Summary: Header files and static libraries for GROMACS
Group: Applications/Science
Prefix: %{prefix}
Requires: fftw-devel >= 2.1.3, gromacs = %{version}-%{release}
%description devel
This package contains header files, static libraries,
and a program example for the GROMACS molecular
dynamics software. You need it if you want to write your
own analysis programs. 

%prep
%setup

%build
# Use the standard /usr/local setup on linux, even if that's slightly
# different from the normal gromacs directory standard. Don't use
# the automatic gromacs architecture exec-prefix.
# Since 'gromacs' isnt present in the prefix it will be added to datadir
# and includedir.
# (This way the package won't interfere with a manual gromacs installation)
# dont use motif since it is not standard on linux.
./configure --enable-shared --prefix=%{prefix} --exec-prefix=%{prefix} --without-motif-libraries
make 

%install
make DESTDIR=${RPM_BUILD_ROOT} install

%clean
rm -rf ${RPM_BUILD_ROOT}

%post
#
# Add our (final) library directory to /etc/ld.so.conf if it is not already there
#
if test -z `grep ${RPM_INSTALL_PREFIX}/lib  /etc/ld.so.conf`; then
     cat >> /etc/ld.so.conf < ${RPM_INSTALL_PREFIX}/lib
fi
# run ldconfig to update the runtime linker database with the new libraries
# (make sure /sbin is in the $PATH)
PATH="/sbin:$PATH" ldconfig

%postun
# after uninstall, run ldconfig to remove the libs from the linker database
PATH="/sbin:$PATH" ldconfig

%files
%defattr(-,root,root)
%{prefix}/bin/*
%{prefix}/share/gromacs/top/*
%{prefix}/share/gromacs/tutor/*
%docdir %{prefix}/share/gromacs/html
%{prefix}/share/gromacs/html/
%{prefix}/man/*
%{prefix}/lib/libgmx.so.2.0.0
%{prefix}/lib/libgmx.so.2
%{prefix}/lib/libmd.so.2.0.0
%{prefix}/lib/libmd.so.2
%files devel
%defattr(-,root,root)
%{prefix}/share/gromacs/template/*
%{prefix}/lib/libgmx.so
%{prefix}/lib/libgmx.a
%{prefix}/lib/libgmx.la
%{prefix}/lib/libmd.so
%{prefix}/lib/libmd.a
%{prefix}/lib/libmd.la
%{prefix}/include/gromacs/*




