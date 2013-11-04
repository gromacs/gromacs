Name:		gromacs
Version:	4.5
Release:	7%{?dist}
Summary:	GROMACS binaries
Group:		Applications/Engineering
License:	LGPLv2.1
URL:		http://www.gromacs.org
Source0:	ftp://ftp.gromacs.org/pub/gromacs/gromacs-%{version}.tar.gz
Source1:	ftp://ftp.gromacs.org/pub/manual/manual-4.0.pdf
Source2:	gromacs-template-makefile-single
Source3:	gromacs-template-makefile-double
Source4:	gromacs-template-makefile-mpi-single
Source5:	gromacs-template-makefile-mpi-double
Source6:	gromacs-README.fedora

# Add shebangs to scripts
Patch0:		gromacs-GMXRC.patch

BuildRoot:	%{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
Requires:	gromacs-common  = %{version}-%{release}

BuildRequires:	fftw-devel
BuildRequires:	gsl-devel
BuildRequires:	openmpi-devel

%if 0%{?rhel} == 4
BuildRequires:	blas
BuildRequires:	lapack
BuildRequires:	xorg-x11-devel
%else
BuildRequires:	blas-devel
BuildRequires:	lapack-devel
BuildRequires:	libX11-devel
%endif

# Check for mpi-selector or environment-modules

%define selector 0
%define modules 0

%if 0%{?fedora} > 9
%define modules 1
%endif

%if 0%{?rhel} == 4
%define selector 1
%endif

%if 0%{?rhel} == 5
%define selector 1
%endif

%if %modules == 1
BuildRequires:	environment-modules
%endif

%if %selector == 1
BuildRequires:	mpi-selector
%endif


%description
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides single and double precision binaries.
The documentation is in the package gromacs-common.

N.B. All binaries have names starting with g_, for example
mdrun has been renamed to g_mdrun.

%package libs
Summary:	GROMACS libraries
Group:		Applications/Engineering
Requires:	gromacs-common = %{version}-%{release}
# Need to have this so that yum doesn't pull atlas instead
Requires:	blas
Requires:	lapack

%description libs
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides runtime libraries needed for the
single and double precision binaries.


%package mpi
Summary:	GROMACS MPI binaries
Group:		Applications/Engineering
Requires:	gromacs-common = %{version}-%{release}
# Need to have this so that yum doesn't pull atlas instead
Requires:	blas
Requires:	lapack

%description mpi
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides MPI single precision and double
precision binaries.


%package common
Summary:	GROMACS shared data and documentation
Group:		Applications/Engineering

%description common
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package includes architecture independent data and
documentation.


%package devel
Summary:	GROMACS header files and development libraries
Group:		Applications/Engineering
Requires:	gromacs-common = %{version}-%{release}
Requires:	gromacs-libs = %{version}-%{release}

%description devel
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package contains header files, development libraries,
and a program example for the GROMACS molecular
dynamics software. You need it if you want to write your
own analysis programs.


%package mpi-devel
Summary:	GROMACS MPI development libraries
Group:		Applications/Engineering
Requires:	gromacs-mpi-libs = %{version}-%{release}
Requires:	gromacs-devel =  %{version}-%{release}
# Need to have this so that yum doesn't install LAM instead
Requires:	openmpi

%description mpi-devel
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package contains development libraries for GROMACS MPI.
You need it if you want to write your own analysis programs.


%package mpi-libs
Summary:	GROMACS libraries
Group:		Applications/Engineering
Requires:	gromacs-common = %{version}-%{release}
# Need to have this so that yum doesn't install LAM instead
Requires:	openmpi
# Need to have this so that yum doesn't pull atlas instead
Requires:	blas
Requires:	lapack

%description mpi-libs
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides runtime libraries needed for the
MPI single and double precision binaries.


%package bash
Summary:	GROMACS bash completion
Group:		Applications/Engineering
Requires:	bash-completion

%description bash
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides the needed 
bash completion for GROMACS


%package zsh
Summary:	GROMACS zsh support
Group:		Applications/Engineering
Requires:	zsh

%description zsh
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides scripts needed to run GROMACS with
zsh, also it provides zsh completion.


%package csh
Summary:	GROMACS csh support
Group:		Applications/Engineering
Requires:	csh

%description csh
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides scripts needed to run GROMACS with
csh and a completion script.

%package tutor
Summary:	GROMACS tutorial files
Group:		Applications/Engineering
Requires:	gromacs-common = %{version}-%{release}

%description tutor
GROMACS is a versatile and extremely well optimized package
to perform molecular dynamics computer simulations and
subsequent trajectory analysis. It is developed for
biomolecules like proteins, but the extremely high
performance means it is used also in several other field
like polymer chemistry and solid state physics.

This package provides tutorials for the use of GROMACS.

%prep
%setup -q
%patch0 -p1

# Fix incorrect permission
chmod a-x src/tools/gmx_xpm2ps.c



%build
# Assembly kernels haven't got .note.GNU-stack sections
# because of incompatibilies with Microsoft Assembler.
# Add noexecstack to compiler flags

export CFLAGS="%optflags -Wa,--noexecstack -fPIC"
export LIBS="-lblas -llapack"

# Single precision
mkdir single
cd single
ln -s ../configure .
%configure --enable-shared \
	--disable-static --enable-float \
	--with-external-blas --with-external-lapack \
	--with-gsl --with-x
sed -i 's|^hardcode_libdir_flag_spec=.*|hardcode_libdir_flag_spec=""|g' libtool
sed -i 's|^runpath_var=LD_RUN_PATH|runpath_var=DIE_RPATH_DIE|g' libtool

make %{?_smp_mflags}
cd ..

# Double precision
mkdir double
cd double
ln -s ../configure .
%configure --disable-rpath --enable-shared \
	--disable-static --disable-float \
	--with-external-blas --with-external-lapack \
	--with-gsl --with-x \
	--program-suffix=_d
sed -i 's|^hardcode_libdir_flag_spec=.*|hardcode_libdir_flag_spec=""|g' libtool
sed -i 's|^runpath_var=LD_RUN_PATH|runpath_var=DIE_RPATH_DIE|g' libtool

make %{?_smp_mflags}
cd ..

# Load MPI enviroment

%if %modules == 1
. /etc/profile.d/modules.sh
module load %{_libdir}/openmpi/*/openmpi.module
%endif

%if %selector == 1
# Set MPI environment
mpi-selector --set `mpi-selector --list | grep openmpi`
source /etc/profile.d/mpi-selector.sh
%endif


# MPI, single precision

mkdir mpi-single
cd mpi-single
ln -s ../configure .
%configure --enable-shared \
	--disable-static --enable-float \
	--with-external-blas --with-external-lapack \
	--with-gsl --with-x --enable-mpi \
	--program-suffix=_mpi
sed -i 's|^hardcode_libdir_flag_spec=.*|hardcode_libdir_flag_spec=""|g' libtool
sed -i 's|^runpath_var=LD_RUN_PATH|runpath_var=DIE_RPATH_DIE|g' libtool

make %{?_smp_mflags} mdrun
#make %{?_smp_mflags}
cd ..

# MPI, double precision
mkdir mpi-double
cd mpi-double
ln -s ../configure .
%configure --enable-shared \
	--disable-static --disable-float \
	--with-external-blas --with-external-lapack \
	--with-gsl --with-x --enable-mpi \
	--program-suffix=_mpi_d
sed -i 's|^hardcode_libdir_flag_spec=.*|hardcode_libdir_flag_spec=""|g' libtool
sed -i 's|^runpath_var=LD_RUN_PATH|runpath_var=DIE_RPATH_DIE|g' libtool

make %{?_smp_mflags} mdrun
#make %{?_smp_mflags}
cd ..


%install
rm -rf %{buildroot}

# Single precision
cd single
make DESTDIR=%{buildroot} INSTALL="install -p" install
cd ..

# Double precision
cd double
make DESTDIR=%{buildroot} INSTALL="install -p" install
cd ..



# MPI, single precision
cd mpi-single
make DESTDIR=%{buildroot} INSTALL="install -p" install-mdrun
cd ..

# MPI, double precision
cd mpi-double
make DESTDIR=%{buildroot} INSTALL="install -p" install-mdrun
cd ..

# Install manual & packager's note
install -cpm 644 %{SOURCE1} .
install -cpm 644 %{SOURCE6} README.fedora

# Remove broken makefiles generated by build process
rm -rf %{buildroot}%{_datadir}/%{name}/template/Makefil*
# Install template makefiles
install -cpm 644 %{SOURCE2} %{buildroot}%{_datadir}/%{name}/template/Makefile.single
install -cpm 644 %{SOURCE3} %{buildroot}%{_datadir}/%{name}/template/Makefile.double
install -cpm 644 %{SOURCE4} %{buildroot}%{_datadir}/%{name}/template/Makefile.mpi.single
install -cpm 644 %{SOURCE5} %{buildroot}%{_datadir}/%{name}/template/Makefile.mpi.double


# Fix GMXRC file permissions
chmod a+x %{buildroot}%{_bindir}/GMXRC %{buildroot}%{_bindir}/GMXRC.*

# Rename binaries and man pages to prevent clashes
# (This is done here so that we don't need to mess with machine generated makefiles.
for bin in anadock do_dssp editconf eneconv genbox genconf genion genrestr gmxcheck gmxdump grompp highway luck make_edi make_ndx mdrun mk_angndx ngmx pdb2gmx protonate sigeps tpbconv trjcat trjconv trjorder wheel x2top xpm2ps xrama ; do 
mv %{buildroot}%{_bindir}/${bin} %{buildroot}%{_bindir}/g_${bin}
mv %{buildroot}%{_bindir}/${bin}_d %{buildroot}%{_bindir}/g_${bin}_d
done

for bin in demux.pl xplor2gmx.pl; do
mv %{buildroot}%{_bindir}/$bin %{buildroot}%{_bindir}/g_${bin}
done

# MPI-enabled binaries (list will continue when the makefile has
# the possibility to compile all mpi-enabled files
for mpibin in mdrun; do
mv %{buildroot}%{_bindir}/${mpibin}_mpi %{buildroot}%{_bindir}/g_${mpibin}_mpi
mv %{buildroot}%{_bindir}/${mpibin}_mpi_d %{buildroot}%{_bindir}/g_${mpibin}_mpi_d
done

# Man pages
for bin in anadock do_dssp editconf eneconv genbox genconf genion genrestr gmxcheck gmxdump grompp highway make_edi make_ndx mdrun mk_angndx ngmx pdb2gmx protonate sigeps tpbconv trjcat trjconv trjorder wheel x2top xpm2ps xrama ; do 
mv %{buildroot}%{_mandir}/man1/${bin}.1 %{buildroot}%{_mandir}/man1/g_${bin}.1
mv %{buildroot}%{_mandir}/man1/${bin}_d.1 %{buildroot}%{_mandir}/man1/g_${bin}_d.1
done

# Move completion files around
chmod a-x %{buildroot}%{_bindir}/completion.*
# Zsh
mkdir -p %{buildroot}%{_datadir}/zsh/site-functions
mv %{buildroot}%{_bindir}/completion.zsh %{buildroot}%{_datadir}/zsh/site-functions/gromacs
# Bash
mkdir -p %{buildroot}%{_sysconfdir}/bash_completion.d
mv %{buildroot}%{_bindir}/completion.bash %{buildroot}/etc/bash_completion.d/gromacs
# Tcsh
mv %{buildroot}%{_bindir}/completion.csh . 

# Remove .la files
rm -rf %{buildroot}/%{_libdir}/*.la

# Post install for libs

%post libs -p /sbin/ldconfig

%postun libs -p /sbin/ldconfig

%post mpi-libs -p /sbin/ldconfig

%postun mpi-libs -p /sbin/ldconfig


%clean
rm -rf %{buildroot}




# Files section

%files
%defattr(-,root,root,-)
%{_bindir}/*
%exclude %{_bindir}/g_mdrun_mpi
%exclude %{_bindir}/g_mdrun_mpi_d
%exclude %{_bindir}/GMXRC*

%files libs
%defattr(-,root,root,-)
%{_libdir}/libgmx.so.*
%{_libdir}/libgmx_d.so.*
%{_libdir}/libgmxana.so.*
%{_libdir}/libgmxana_d.so.*
%{_libdir}/libmd.so.*
%{_libdir}/libmd_d.so.*

%files mpi
%defattr(-,root,root,-)
%{_bindir}/g_mdrun_mpi
%{_bindir}/g_mdrun_mpi_d


%files mpi-libs
%defattr(-,root,root,-)
%{_libdir}/libgmx_mpi.so.*
%{_libdir}/libgmx_mpi_d.so.*
%{_libdir}/libmd_mpi.so.*
%{_libdir}/libmd_mpi_d.so.*



%files common
%defattr(-,root,root,-)
%doc AUTHORS COPYING README manual-4.0.pdf README.fedora
%{_bindir}/GMXRC
%{_bindir}/GMXRC.bash
%{_mandir}/man1/*
%{_datadir}/%{name}
%exclude %{_datadir}/%{name}/template
%exclude %{_datadir}/%{name}/tutor

%files devel
%defattr(-,root,root,-)
%{_includedir}/%{name}
%{_libdir}/libgmx.so
%{_libdir}/libgmx_d.so
%{_libdir}/libgmxana.so
%{_libdir}/libgmxana_d.so
%{_libdir}/libmd.so
%{_libdir}/libmd_d.so
%{_datadir}/%{name}/template
%exclude %{_datadir}/%{name}/template/Makefile.mpi.*

%files mpi-devel
%defattr(-,root,root,-)
%{_libdir}/libgmx_mpi.so
%{_libdir}/libgmx_mpi_d.so
%{_libdir}/libmd_mpi.so
%{_libdir}/libmd_mpi_d.so
%{_datadir}/%{name}/template/Makefile.mpi.*


%files zsh
%defattr(-,root,root,-)
%{_datadir}/zsh/site-functions/gromacs
%{_bindir}/GMXRC.zsh

%files bash
%defattr(-,root,root,-)
%config(noreplace) %{_sysconfdir}/bash_completion.d/gromacs


%files csh
%defattr(-,root,root,-)
%doc completion.csh
%{_bindir}/GMXRC.csh

%files tutor
%defattr(-,root,root,-)
%{_datadir}/%{name}/tutor


%changelog
* Wed Jan 14 2009 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0.2-7
- Update manual to latest version.
- Removed Requires: blas and lapack.

* Mon Nov 10 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0.2-6
- Update to 4.0.2.

* Sun Nov 09 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0.1-5
- Add Requires: blas too.

* Sun Nov 09 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0.1-4
- Update to 4.0.1.
- Add Requires: lapack and openmpi to prevent yum from pulling atlas and lam
instead.

* Wed Oct 15 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-3
- Rename also man pages.

* Mon Oct 13 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-2
- Added noreplace to bash completion file.
- Changed double precision mpi binary suffix to _mpi_d.

* Sun Oct 12 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-1
- Update to Gromacs 4.0.
- Remove module system and patch file names to begin with g_.

* Wed Oct 08 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.15.rc3
- Changed location of binaries.
- Removed conflict of module file, as the program is binary compatible with older versions.

* Wed Oct 08 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.14.rc3
- The gromacs module is loaded automatically and it conflicts with gromacs3.

* Tue Oct 07 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.13.rc3
- Renamed module files from %%{name}-%%{version} to %%{name}.

* Mon Oct 06 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.12.rc3
- Fix BR to get GROMACS to build in mock for epel-4.

* Sat Oct 04 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.11.rc3
- Fix to get GROMACS to build in mock for epel-5.

* Sat Oct 04 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.10.rc3
- Implement module system & remove binary renaming.
- No need for autoreconf anymore.
- Update to rc3.

* Sat Oct 04 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.9.rc2
- Fall back to autoreconf due to binary renaming.

* Fri Oct 03 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.8.rc2
- Modified install commands to preserve timestamps.

* Fri Oct 03 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.7.rc2
- Even more review fixes.
- Binaries renamed:
	highway	->	g_highway
	luck	->	g_luck
	sigeps	->	g_sigeps
	wheel	->	g_wheel

* Thu Oct 02 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.6.rc2
- Final review fixes.

* Wed Oct 01 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.5.rc2
- Strip down requires by branching tutor to its own package.

* Tue Sep 30 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.4.rc2
- Extensive package review fixes.
- Unclear licenses on some files, filed upstream bug 217.
  http://bugzilla.gromacs.org/show_bug.cgi?id=217

* Mon Sep 29 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.3.rc2
- Move .so files to -devel package.
- Remove .la files.

* Mon Sep 29 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.2.rc2
- Implement out-of-tree-builds.
- Add --noexecstack to CFLAGS.
- Remove execstack procedure and prelink from buildreqs.
- Filed upstream bug 215 to add .note.GNU-stack .
- Fix incorrect file permission on src/tools/gmx_xpm2ps.c .

* Mon Sep 29 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.1.rc2
- Alphabetized buildrequires.
- Changed gromacs-share to gromacs-common.

* Fri Sep 26 2008 Jussi Lehtola <jussi.lehtola@iki.fi> - 4.0-0.0.rc2
- Initial build.
