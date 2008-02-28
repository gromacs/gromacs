#
# RPM specification file to make gromacs-mpi mdrun executable,
# and libraries. For the rest, use the gromacs (non-mpi) spec file.

#
# Main package - only dynamic libs, and no header files
#
Summary: Molecular dynamics package (parallel)
Name: gromacs-mpi
Version: 3.3.3
Release: 1
License: GPL
Group: Applications/Science
Buildroot: %{_topdir}/buildroot
Requires: fftw3 >= 3.0.1 , lam, gromacs = %{version}-%{release}
Source: ftp://ftp.gromacs.org/pub/gromacs/source/gromacs-%{version}.tar.gz
URL: http://www.gromacs.org
Packager: Erik Lindahl <lindahl@gromacs.org>
%description
This is the MPI support files and mdrun executable of GROMACS, 
a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. All auxiliary programs
and data files are located in the gromacs (non-mpi) package.

#
# The header files and static libraries go into gromacs-mpi-devel...
#
%package devel
Summary: Header files and static libs for parallel GROMACS
Group: Applications/Science
Requires: fftw3 >= 3.0.1, fftw3-devel >= 3.0.1, lam, gromacs = %{version}-%{release}, gromacs-devel = %{version}-%{release}, gromacs-mpi = %{version}-%{release}
%description devel
This package contains the static libraries for
the parallel GROMACS development. You will only need
it if you are hacking parallel mdrun stuff, and then 
you probably want the full source anyway...

%prep
%setup -n gromacs-%{version}

%build
# Call it mdrun_mpi
%configure     --enable-shared \
               --without-motif-libraries \
	       --enable-mpi \
	       --program-suffix=_mpi
make %{?_smp_mflags} mdrun

%install
make DESTDIR=${RPM_BUILD_ROOT} install-mdrun

%clean
rm -rf ${RPM_BUILD_ROOT}

%post

%postun


%files
%defattr(-,root,root)
%{_bindir}/mdrun_mpi
%{_libdir}/*_mpi.so.*


%files devel
%defattr(-,root,root)
%exclude %{_libdir}/*.la
%{_libdir}/*_mpi.a
%{_libdir}/*_mpi.so













