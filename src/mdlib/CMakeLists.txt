#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
#
# The nbnxn kernel files take the longest time to compile, so we sneak
# them in first to take best advantage of make-time parallelisation.
file(GLOB MDLIB_SOURCES nbnxn_kernels/simd_4xn/*.c nbnxn_kernels/simd_2xnn/*.c nbnxn_kernels/*.c *.c)

if(GMX_GPU)
    add_subdirectory(nbnxn_cuda)
    set(GMX_GPU_LIBRARIES nbnxn_cuda)
endif()

# apply gcc 4.4.x bug workaround
if(GMX_USE_GCC44_BUG_WORKAROUND)
   include(gmxGCC44O3BugWorkaround)
   gmx_apply_gcc44_bug_workaround("force.c")
   gmx_apply_gcc44_bug_workaround("constr.c")
endif()

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library

if(NOT GMX_FFT_FFTPACK)
list(REMOVE_ITEM MDLIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/fftpack.c)
endif()

add_library(md ${MDLIB_SOURCES})

target_link_libraries(md ${GMX_GPU_LIBRARIES} gmx)
if(GMX_BUILD_OWN_FFTW)
    # Only needed for cmake 2.8.7, otherwise this should be automatic
    # This dependency has to be made here rather than the CMakeLists.txt that
    # does the FFTW build, because of the order in which
    # add_subdirectory() calls are made in the top-level CMakeLists.txt; the
    # md library target does not necessarily exist yet. Also enabling and
    # disabling GMX_BUILD_OWN_FFTW changes dependencies correctly.
    add_dependencies(md gmxfftw)
endif()
set_target_properties(md PROPERTIES OUTPUT_NAME "${MD_PREFIX}md${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}"
    COMPILE_FLAGS "${OpenMP_C_FLAGS}")

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmd.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "lib${MD_PREFIX}md${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
