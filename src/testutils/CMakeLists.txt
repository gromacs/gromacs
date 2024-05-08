#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2011- The GROMACS Authors
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

if (NOT GMX_BUILD_UNITTESTS)
    gmx_add_missing_tests_notice("Unit tests have not been run. You need to set GMX_BUILD_UNITTESTS=ON if you want to build and run them.")
    return()
endif()

if (GMX_GPU_CUDA)
    add_library(testutils)

    #we cannot link with both cufftmp and cufft
    if(GMX_USE_cuFFTMp)
        target_link_libraries(testutils PRIVATE ${cuFFTMp_LIBRARY})
    else()
        target_link_libraries(testutils PRIVATE CUDA::cufft)
    endif()
    if (GMX_NVSHMEM)
        target_link_libraries(testutils PUBLIC nvshmem_host_lib)
    endif()
elseif (GMX_GPU_HIP)
    set_source_files_properties(test_device.cpp PROPERTIES LANGUAGE HIP)
    gmx_hip_add_library(testutils test_device.cpp)

    target_link_libraries(testutils PRIVATE hip::host)
    if (GMX_GPU_FFT_VKFFT)
        target_link_libraries(testutils PRIVATE VkFFT)
    elseif(GMX_GPU_FFT_HIPFFT)
        target_link_libraries(testutils PRIVATE roc::rocfft hip::hipfft)
    elseif(GMX_GPU_FFT_ROCFFT)
        target_link_libraries(testutils PRIVATE roc::rocfft)
    endif()
else()
    add_library(testutils STATIC ${UNITTEST_TARGET_OPTIONS})
endif()

# Module interface / provided facilities
target_include_directories(testutils PUBLIC include)

# Executable targets for tests based on `testutils` acquire the source for
# their entry point from unittest_main.cpp when linking to the `testutils` target.
target_sources(testutils INTERFACE unittest_main.cpp)


target_sources(testutils PRIVATE
               cmdlinetest.cpp
               conftest.cpp
               filematchers.cpp
               interactivetest.cpp
               loggertest.cpp
               mpi_printer.cpp
               mpitest.cpp
               refdata.cpp
               refdata_xml.cpp
               simulationdatabase.cpp
               stdiohelper.cpp
               stringtest.cpp
               testasserts.cpp
               testfilemanager.cpp
               testfileredirector.cpp
               test_device.cpp
               test_hardware_environment.cpp
               testinit.cpp
               testmatchers.cpp
               testoptions.cpp
               textblockmatchers.cpp
               topologyhelpers.cpp
               tprfilegenerator.cpp
               trajectoryreader.cpp
               xvgtest.cpp
               )


if(HAVE_TINYXML2)
    target_include_directories(testutils SYSTEM PRIVATE ${TinyXML2_INCLUDE_DIR})
    target_link_libraries(testutils PRIVATE ${TinyXML2_LIBRARIES})
else()
    target_include_directories(testutils SYSTEM BEFORE PRIVATE ${CMAKE_SOURCE_DIR}/src/external/tinyxml2)
    target_sources(testutils PRIVATE ${CMAKE_SOURCE_DIR}/src/external/tinyxml2/tinyxml2.cpp)
endif()

target_include_directories(testutils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (GMX_GPU_SYCL)
    add_sycl_to_target(
        TARGET testutils
        SOURCES test_device.cpp
        )
endif()

gmx_target_compile_options(testutils)
target_compile_definitions(testutils PRIVATE HAVE_CONFIG_H)
target_include_directories(testutils SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
target_compile_definitions(testutils PRIVATE TMPI_USE_VISIBILITY)
target_link_libraries(testutils PRIVATE libgromacs ${GMX_COMMON_LIBRARIES} gmock)
target_link_libraries(testutils PUBLIC common)

# GROMACS module dependencies.
# Note that testutils conveys transitive dependencies on some modules.
# TODO: Remove `legacy_modules` once specific modules are explicitly linked.
target_link_libraries(testutils PUBLIC legacy_modules)
target_link_libraries(testutils PUBLIC utility math)
target_link_libraries(testutils PRIVATE topology)

add_subdirectory(tests)
