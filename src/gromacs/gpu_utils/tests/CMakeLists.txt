#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2017- The GROMACS Authors
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

# Arrange to compile files with test infrastructure as normal C++, and
# other files that sometimes have platform-specific symbols using the
# required infrastructure to compile and link, without needing three
# different files and strategies for the different GPU implementation
# flavours.

gmx_add_unit_test(GpuUtilsUnitTests gpu_utils-test HARDWARE_DETECTION
    CPP_SOURCE_FILES
        # Tests of code
        clfftinitializer.cpp
        device_availability.cpp
        device_stream_manager.cpp
        pinnedmemorychecker.cpp
        
    GPU_CPP_SOURCE_FILES
        device_buffer.cpp
        gpueventsynchronizer.cpp
        hostallocator.cpp
        typecasts.cpp

    CUDA_CU_SOURCE_FILES
        devicetransfers.cu
        typecasts_runner.cu

    HIP_CPP_SOURCE_FILES
        devicetransfers_hip.cpp
        typecasts_runner_hip.cpp

    OPENCL_CPP_SOURCE_FILES
        devicetransfers_ocl.cpp

    SYCL_CPP_SOURCE_FILES
        devicetransfers_sycl.cpp

    NON_GPU_CPP_SOURCE_FILES
        devicetransfers.cpp
    )
target_link_libraries(gpu_utils-test PRIVATE gpu_utils)
target_link_libraries(gpu_utils-test PRIVATE math)

if (GMX_NVSHMEM)
    gmx_add_mpi_unit_test(NvshmemTests nvshmem-test 2 HARDWARE_DETECTION NVSHMEM
    CUDA_CU_SOURCE_FILES
        nvshmem_simple_put.cu
    GPU_CPP_SOURCE_FILES
        nvshmem_test.cpp)

    target_include_directories(nvshmem-test PRIVATE ${PROJECT_SOURCE_DIR}/src/gromacs/math/tests)
    target_link_libraries(nvshmem-test PRIVATE gpu_utils)
endif()
