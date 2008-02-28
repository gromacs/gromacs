#
# RPM specification file to make normal gromacs packages (without mpi)
# If you have mpi installed, you can create an mpi mdrun executable
# and libs with the gromacs-mpi spec file!


#
# Main package - only dynamic libs, and no header files
#
Summary: Molecular dynamics package (non-parallel version)
Name: gromacs
Version: 3.3.3
Release: 1
License: GPL
Group: Applications/Science
Buildroot: %{_topdir}/buildroot
Requires: fftw3 >= 3.0.1
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
libs in gromacs-dev. 

#
# The header files and static libraries go into gromacs-devel...
#
%package devel
Summary: Header files and static libraries for GROMACS
Group: Applications/Science
Requires: fftw3-devel >= 3.0.1, gromacs = %{version}-%{release}
%description devel
This package contains header files, static libraries,
and a program example for the GROMACS molecular
dynamics software. You need it if you want to write your
own analysis programs. 

%prep
%setup

%build
%configure 	--enable-shared \
	   	--without-motif-libraries 
make %{?_smp_mflags}

%install
make DESTDIR=${RPM_BUILD_ROOT} install

%clean
rm -rf ${RPM_BUILD_ROOT}

%post
# add libraries to system database
# (make sure /sbin is in the $PATH)
PATH="/sbin:$PATH" ldconfig

%postun
# after uninstall, run ldconfig to remove the libs from the linker database
PATH="/sbin:$PATH" ldconfig

%files
%defattr(-,root,root)
%{_bindir}/*
%{_libdir}/*.so.*
%{_mandir}/*
%{_datadir}/gromacs/html/*
%{_datadir}/gromacs/top/*
%{_datadir}/gromacs/tutor/*
%{_datadir}/gromacs/template/*

%files devel
%defattr(-,root,root)
%exclude %{_libdir}/*.la
%{_includedir}/gromacs/*
%{_libdir}/*.a
%{_libdir}/*.so

