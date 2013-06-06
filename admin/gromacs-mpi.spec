#
# RPM specification file to make gromacs-mpi mdrun executable,
# and libraries. For the rest, use the gromacs (non-mpi) spec file.

#
# Main package - only dynamic libs, and no header files
#
Summary: Molecular dynamics package (parallel)
Name: gromacs-mpi
Version: 4.5
Release: 1
Copyright: LGPLv2.1
Group: Applications/Science
Prefix: /usr/local
Buildroot: %{_topdir}/buildroot
Requires: fftw-lammpi >= 2.1.3 , lam, gromacs = %{version}-%{release}
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
Prefix: %{prefix}
Requires: fftw-lammpi >= 2.1.3, fftw-lammpi-devel >= 2.1.3, lam, gromacs = %{version}-%{release}, gromacs-devel = %{version}-%{release}, gromacs-mpi = %{version}-%{release}
%description devel
This package contains the static libraries for
the parallel GROMACS development. You will only need
it if you are hacking parallel mdrun stuff, and then 
you probably want the full source anyway...

%prep
%setup -n gromacs-%{version}

%build
# Call it mdrun_mpi
./configure --enable-shared --enable-mpi --prefix=%{prefix} --exec-prefix=%{prefix} --program-suffix=_mpi --without-motif-libraries
make mdrun

%install
make DESTDIR=${RPM_BUILD_ROOT} install-mdrun

%clean
rm -rf ${RPM_BUILD_ROOT}

%post

%postun


%files 
%defattr(-,root,root)
%{prefix}/bin/mdrun_mpi
%{prefix}/lib/libgmx_mpi.so.3.0.0
%{prefix}/lib/libgmx_mpi.so.3
%{prefix}/lib/libmd_mpi.so.3.0.0
%{prefix}/lib/libmd_mpi.so.3
%{prefix}/lib/libgmxana_mpi.so.3.0.0
%{prefix}/lib/libgmxana_mpi.so.3
%files devel
%defattr(-,root,root)
%{prefix}/lib/libgmx_mpi.so
%{prefix}/lib/libgmx_mpi.a
%{prefix}/lib/libgmx_mpi.la
%{prefix}/lib/libmd_mpi.so
%{prefix}/lib/libmd_mpi.a
%{prefix}/lib/libmd_mpi.la
%{prefix}/lib/libgmxana_mpi.so
%{prefix}/lib/libgmxana_mpi.a
%{prefix}/lib/libgmxana_mpi.la













