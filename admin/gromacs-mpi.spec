#
# RPM specification file to make gromacs-mpi mdrun executable,
# and libraries. For the rest, use the gromacs (non-mpi) spec file.

#
# Main package - only dynamic libs, and no header files
#
Summary: Molecular dynamics package (parallel)
Name: gromacs-mpi
Version: 3.0
Release: 1
Copyright: GPL
Group: Applications/Science
Requires: fftw-lammpi >= 2.1.3 , lam = 6.5.2, gromacs = %{version}-%{release}
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
Requires: fftw-lammpi >= 2.1.3, fftw-lammpi-devel >= 2.1.3, lam = 6.5.2, gromacs = %{version}-%{release}, gromacs-devel = %{version}-%{release}, gromacs-mpi = %{version}-%{release}
%description devel
This package contains the static libraries for
the parallel GROMACS development. You will only need
it if you are hacking parallel mdrun stuff, and then 
you probably want the full source anyway...

%prep
%setup

%build
# Call it mdrun_mpi
./configure --enable-shared --program_suffix="_mpi" 
make mdrun

%install
make DESTDIR=${RPM_BUILD_ROOT} install-mdrun

%clean
rm -rf ${RPM_BUILD_ROOT}

%post
# /etc/ld.so.conf should have been updated by the normal gromacs package.
# Overwrite the mdrun link - it should point to mdrun_mpi iso mdrun_nompi now!
(cd ${RPM_INSTALL_PREFIX}/%{_host}/bin && ln -sf mdrun_mpi mdrun)


%postun
# If we removed the gromacs-mpi package, while the non-mpi version is still present,
# we should restore the mdrun link:
(cd ${RPM_INSTALL_PREFIX}/%{_host}/bin && test ! -e mdrun && ln -s mdrun_nompi mdrun)

%files 
%defattr(-,root,root)
/usr/local/gromacs/%{_host}/bin/mdrun_mpi
/usr/local/gromacs/%{_host}/lib/libgmx_mpi.so.1.0.0
/usr/local/gromacs/%{_host}/lib/libgmx_mpi.so.1
/usr/local/gromacs/%{_host}/lib/libmd_mpi.so.1.0.0
/usr/local/gromacs/%{_host}/lib/libmd_mpi.so.1
%files devel
%defattr(-,root,root)
/usr/local/gromacs/%{_host}/lib/libgmx_mpi.so
/usr/local/gromacs/%{_host}/lib/libgmx_mpi.a
/usr/local/gromacs/%{_host}/lib/libgmx_mpi.la
/usr/local/gromacs/%{_host}/lib/libmd_mpi.so
/usr/local/gromacs/%{_host}/lib/libmd_mpi.a
/usr/local/gromacs/%{_host}/lib/libmd_mpi.la