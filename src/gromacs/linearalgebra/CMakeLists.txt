#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

add_library(linearalgebra OBJECT)
target_sources(linearalgebra PRIVATE
        binary_information.cpp
        eigensolver.cpp
        gmx_arpack.cpp
        matrix.cpp
        sparsematrix.cpp
)

# Note that the use of MKL for FFT determines that MKL will be used
# also for linear algebra, for simplicity.
set_source_files_properties(binary_information.cpp PROPERTIES COMPILE_FLAGS
    "-DGMX_DESCRIBE_BLAS=\"${GMX_DESCRIBE_BLAS}\" -DGMX_DESCRIBE_LAPACK=\"${GMX_DESCRIBE_LAPACK}\" -DGMX_FFT_MKL=${GMX_FFT_MKL}")

if (NOT GMX_EXTERNAL_BLAS)
    set(BLAS_SOURCES
        gmx_blas/dasum.cpp
        gmx_blas/daxpy.cpp
        gmx_blas/dcopy.cpp
        gmx_blas/ddot.cpp
        gmx_blas/dgemm.cpp
        gmx_blas/dgemv.cpp
        gmx_blas/dger.cpp
        gmx_blas/dnrm2.cpp
        gmx_blas/drot.cpp
        gmx_blas/dscal.cpp
        gmx_blas/dswap.cpp
        gmx_blas/dsymv.cpp
        gmx_blas/dsyr2.cpp
        gmx_blas/dsyr2k.cpp
        gmx_blas/dtrmm.cpp
        gmx_blas/dtrmv.cpp
        gmx_blas/dtrsm.cpp
        gmx_blas/idamax.cpp
        gmx_blas/isamax.cpp
        gmx_blas/sasum.cpp
        gmx_blas/saxpy.cpp
        gmx_blas/scopy.cpp
        gmx_blas/sdot.cpp
        gmx_blas/sgemm.cpp
        gmx_blas/sgemv.cpp
        gmx_blas/sger.cpp
        gmx_blas/snrm2.cpp
        gmx_blas/srot.cpp
        gmx_blas/sscal.cpp
        gmx_blas/sswap.cpp
        gmx_blas/ssymv.cpp
        gmx_blas/ssyr2.cpp
        gmx_blas/ssyr2k.cpp
        gmx_blas/strmm.cpp
        gmx_blas/strmv.cpp
        gmx_blas/strsm.cpp
        )
    target_sources(linearalgebra PRIVATE ${BLAS_SOURCES})
endif()

if (NOT GMX_EXTERNAL_LAPACK)
    set(LAPACK_SOURCES
        gmx_lapack/dbdsdc.cpp
        gmx_lapack/dbdsqr.cpp
        gmx_lapack/dgebd2.cpp
        gmx_lapack/dgebrd.cpp
        gmx_lapack/dgelq2.cpp
        gmx_lapack/dgelqf.cpp
        gmx_lapack/dgeqr2.cpp
        gmx_lapack/dgeqrf.cpp
        gmx_lapack/dgesdd.cpp
        gmx_lapack/dgetf2.cpp
        gmx_lapack/dgetrf.cpp
        gmx_lapack/dgetri.cpp
        gmx_lapack/dgetrs.cpp
        gmx_lapack/dlabrd.cpp
        gmx_lapack/dlacpy.cpp
        gmx_lapack/dlae2.cpp
        gmx_lapack/dlaebz.cpp
        gmx_lapack/dlaed6.cpp
        gmx_lapack/dlaev2.cpp
        gmx_lapack/dlagtf.cpp
        gmx_lapack/dlagts.cpp
        gmx_lapack/dlamrg.cpp
        gmx_lapack/dlange.cpp
        gmx_lapack/dlanst.cpp
        gmx_lapack/dlansy.cpp
        gmx_lapack/dlapy2.cpp
        gmx_lapack/dlar1vx.cpp
        gmx_lapack/dlarf.cpp
        gmx_lapack/dlarfb.cpp
        gmx_lapack/dlarfg.cpp
        gmx_lapack/dlarft.cpp
        gmx_lapack/dlarnv.cpp
        gmx_lapack/dlarrbx.cpp
        gmx_lapack/dlarrex.cpp
        gmx_lapack/dlarrfx.cpp
        gmx_lapack/dlarrvx.cpp
        gmx_lapack/dlartg.cpp
        gmx_lapack/dlaruv.cpp
        gmx_lapack/dlas2.cpp
        gmx_lapack/dlascl.cpp
        gmx_lapack/dlasd0.cpp
        gmx_lapack/dlasd1.cpp
        gmx_lapack/dlasd2.cpp
        gmx_lapack/dlasd3.cpp
        gmx_lapack/dlasd4.cpp
        gmx_lapack/dlasd5.cpp
        gmx_lapack/dlasd6.cpp
        gmx_lapack/dlasd7.cpp
        gmx_lapack/dlasd8.cpp
        gmx_lapack/dlasda.cpp
        gmx_lapack/dlasdq.cpp
        gmx_lapack/dlasdt.cpp
        gmx_lapack/dlaset.cpp
        gmx_lapack/dlasq1.cpp
        gmx_lapack/dlasq2.cpp
        gmx_lapack/dlasq3.cpp
        gmx_lapack/dlasq4.cpp
        gmx_lapack/dlasq5.cpp
        gmx_lapack/dlasq6.cpp
        gmx_lapack/dlasr.cpp
        gmx_lapack/dlasrt.cpp
        gmx_lapack/dlasrt2.cpp
        gmx_lapack/dlassq.cpp
        gmx_lapack/dlasv2.cpp
        gmx_lapack/dlaswp.cpp
        gmx_lapack/dlatrd.cpp
        gmx_lapack/dorg2r.cpp
        gmx_lapack/dorgbr.cpp
        gmx_lapack/dorgl2.cpp
        gmx_lapack/dorglq.cpp
        gmx_lapack/dorgqr.cpp
        gmx_lapack/dorm2l.cpp
        gmx_lapack/dorm2r.cpp
        gmx_lapack/dormbr.cpp
        gmx_lapack/dorml2.cpp
        gmx_lapack/dormlq.cpp
        gmx_lapack/dormql.cpp
        gmx_lapack/dormqr.cpp
        gmx_lapack/dormtr.cpp
        gmx_lapack/dstebz.cpp
        gmx_lapack/dstegr.cpp
        gmx_lapack/dstein.cpp
        gmx_lapack/dsteqr.cpp
        gmx_lapack/dsterf.cpp
        gmx_lapack/dstevr.cpp
        gmx_lapack/dsyevr.cpp
        gmx_lapack/dsytd2.cpp
        gmx_lapack/dsytrd.cpp
        gmx_lapack/dtrti2.cpp
        gmx_lapack/dtrtri.cpp
        gmx_lapack/ilasrt2.cpp
        gmx_lapack/sbdsdc.cpp
        gmx_lapack/sbdsqr.cpp
        gmx_lapack/sgebd2.cpp
        gmx_lapack/sgebrd.cpp
        gmx_lapack/sgelq2.cpp
        gmx_lapack/sgelqf.cpp
        gmx_lapack/sgeqr2.cpp
        gmx_lapack/sgeqrf.cpp
        gmx_lapack/sgesdd.cpp
        gmx_lapack/sgetf2.cpp
        gmx_lapack/sgetrf.cpp
        gmx_lapack/sgetri.cpp
        gmx_lapack/sgetrs.cpp
        gmx_lapack/slabrd.cpp
        gmx_lapack/slacpy.cpp
        gmx_lapack/slae2.cpp
        gmx_lapack/slaebz.cpp
        gmx_lapack/slaed6.cpp
        gmx_lapack/slaev2.cpp
        gmx_lapack/slagtf.cpp
        gmx_lapack/slagts.cpp
        gmx_lapack/slamrg.cpp
        gmx_lapack/slange.cpp
        gmx_lapack/slanst.cpp
        gmx_lapack/slansy.cpp
        gmx_lapack/slapy2.cpp
        gmx_lapack/slar1vx.cpp
        gmx_lapack/slarf.cpp
        gmx_lapack/slarfb.cpp
        gmx_lapack/slarfg.cpp
        gmx_lapack/slarft.cpp
        gmx_lapack/slarnv.cpp
        gmx_lapack/slarrbx.cpp
        gmx_lapack/slarrex.cpp
        gmx_lapack/slarrfx.cpp
        gmx_lapack/slarrvx.cpp
        gmx_lapack/slartg.cpp
        gmx_lapack/slaruv.cpp
        gmx_lapack/slas2.cpp
        gmx_lapack/slascl.cpp
        gmx_lapack/slasd0.cpp
        gmx_lapack/slasd1.cpp
        gmx_lapack/slasd2.cpp
        gmx_lapack/slasd3.cpp
        gmx_lapack/slasd4.cpp
        gmx_lapack/slasd5.cpp
        gmx_lapack/slasd6.cpp
        gmx_lapack/slasd7.cpp
        gmx_lapack/slasd8.cpp
        gmx_lapack/slasda.cpp
        gmx_lapack/slasdq.cpp
        gmx_lapack/slasdt.cpp
        gmx_lapack/slaset.cpp
        gmx_lapack/slasq1.cpp
        gmx_lapack/slasq2.cpp
        gmx_lapack/slasq3.cpp
        gmx_lapack/slasq4.cpp
        gmx_lapack/slasq5.cpp
        gmx_lapack/slasq6.cpp
        gmx_lapack/slasr.cpp
        gmx_lapack/slasrt.cpp
        gmx_lapack/slasrt2.cpp
        gmx_lapack/slassq.cpp
        gmx_lapack/slasv2.cpp
        gmx_lapack/slaswp.cpp
        gmx_lapack/slatrd.cpp
        gmx_lapack/sorg2r.cpp
        gmx_lapack/sorgbr.cpp
        gmx_lapack/sorgl2.cpp
        gmx_lapack/sorglq.cpp
        gmx_lapack/sorgqr.cpp
        gmx_lapack/sorm2l.cpp
        gmx_lapack/sorm2r.cpp
        gmx_lapack/sormbr.cpp
        gmx_lapack/sorml2.cpp
        gmx_lapack/sormlq.cpp
        gmx_lapack/sormql.cpp
        gmx_lapack/sormqr.cpp
        gmx_lapack/sormtr.cpp
        gmx_lapack/sstebz.cpp
        gmx_lapack/sstegr.cpp
        gmx_lapack/sstein.cpp
        gmx_lapack/ssteqr.cpp
        gmx_lapack/ssterf.cpp
        gmx_lapack/sstevr.cpp
        gmx_lapack/ssyevr.cpp
        gmx_lapack/ssytd2.cpp
        gmx_lapack/ssytrd.cpp
        gmx_lapack/strti2.cpp
        gmx_lapack/strtri.cpp
        )
    target_sources(linearalgebra PRIVATE ${LAPACK_SOURCES})
endif()

if (BUILD_SHARED_LIBS)
    set_target_properties(linearalgebra PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(linearalgebra PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories(linearalgebra PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
if (GMX_EXTERNAL_BLAS AND GMX_EXTERNAL_LAPACK AND GMX_FFT_MKL)
    target_link_libraries(linearalgebra PRIVATE MKL::MKL)
endif()

gmx_target_compile_options(linearalgebra)
target_compile_definitions(linearalgebra PRIVATE HAVE_CONFIG_H)
# The linearalgebra code is all considered external, and we will
# not keep it free of warnings. Any compiler suppressions required
# should be added here.
if (MSVC)
    gmx_target_warning_suppression(linearalgebra /wd4305 HAS_NO_MSVC_INITIALIZING_FLOAT_FROM_DOUBLE)
    gmx_target_warning_suppression(linearalgebra /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION_DOUBLE_TO_REAL)
    gmx_target_warning_suppression(linearalgebra /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
else()
    # We usually enable this warning, but here it is being used
    # correctly to construct strings for use with an API that does
    # not expect null termination of C strings.
    gmx_target_warning_suppression(linearalgebra -Wno-stringop-truncation HAS_NO_STRINGOP_TRUNCATION)
endif()
target_link_libraries(linearalgebra PRIVATE legacy_api)
target_link_libraries(linearalgebra PRIVATE common)
# TODO: Remove `legacy_modules` once specific modules are explicitly linked.
target_link_libraries(linearalgebra PRIVATE legacy_modules utility)


if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
