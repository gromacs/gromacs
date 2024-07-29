#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
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

# Set up the module library
add_library(ewald INTERFACE)
gmx_add_libgromacs_sources(
    calculate_spline_moduli.cpp
    ewald.cpp
    ewald_utils.cpp
    long_range_correction.cpp
    pme.cpp
    pme_gather.cpp
    pme_grid.cpp
    pme_load_balancing.cpp
    pme_only.cpp
    pme_pp.cpp
    pme_redistribute.cpp
    pme_solve.cpp
    pme_spline_work.cpp
    pme_spread.cpp
    # Files that implement stubs
    pme_gpu_program.cpp
    pme_pp_comm_gpu_impl.cpp
    pme_coordinate_receiver_gpu_impl.cpp
    pme_force_sender_gpu_impl.cpp
    )
if (GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
        # CUDA-specific sources
        pme_gather.cu
        pme_solve.cu
        pme_spread.cu
        pme_gpu_program_impl.cu
        pme_pp_comm_gpu_impl_gpu.cpp
        pme_pp_comm_gpu_impl_gpu.cu
        pme_force_sender_gpu_impl_gpu.cpp
        pme_force_sender_gpu_impl_gpu.cu
        pme_coordinate_receiver_gpu_impl_gpu.cpp
        pme_gpu_grid.cu
        # GPU-specific sources
        pme_gpu.cpp
        pme_gpu_internal.cpp
        pme_gpu_timings.cpp
        )
    _gmx_add_files_to_property(CUDA_SOURCES
        # Must add these files so they can include device_information.h
        pme_gpu_internal.cpp
        pme_gpu_timings.cpp
        pme_pp_comm_gpu_impl_gpu.cpp
        pme_force_sender_gpu_impl_gpu.cpp
        pme_coordinate_receiver_gpu_impl_gpu.cpp
        )
elseif (GMX_GPU_OPENCL)
    gmx_add_libgromacs_sources(
        # OpenCL-specific sources
        pme_gpu_program_impl_ocl.cpp
        pme_gpu_ocl_stubs.cpp
        # GPU-specific sources
        pme_gpu.cpp
        pme_gpu_internal.cpp
        pme_gpu_timings.cpp
        )
elseif (GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        # GPU-specific sources
        pme_gather_sycl.cpp
        pme_gpu.cpp
        pme_gpu_grid_sycl.cpp
        pme_gpu_internal.cpp
        pme_gpu_program_impl_sycl.cpp
        pme_gpu_timings.cpp
        pme_solve_sycl.cpp
        pme_spread_sycl.cpp
        pme_pp_comm_gpu_impl_gpu_sycl.cpp
        pme_force_sender_gpu_impl_gpu_sycl.cpp
        pme_pp_comm_gpu_impl_gpu.cpp
        pme_force_sender_gpu_impl_gpu.cpp
        pme_force_sender_gpu_impl_gpu_sycl.cpp
        pme_coordinate_receiver_gpu_impl_gpu.cpp
        )
    _gmx_add_files_to_property(SYCL_SOURCES
        pme_gather_sycl.cpp
        pme_gpu_grid_sycl.cpp
        pme_gpu_internal.cpp
        pme_gpu_program.cpp
        pme_gpu_program_impl_sycl.cpp
        pme_gpu_3dfft_sycl.cpp
        pme_gpu_timings.cpp
        pme_pp_comm_gpu_impl_gpu.cpp
        pme_pp_comm_gpu_impl_gpu_sycl.cpp
        pme_force_sender_gpu_impl_gpu.cpp
        pme_force_sender_gpu_impl_gpu_sycl.cpp
        pme_coordinate_receiver_gpu_impl_gpu.cpp
        pme_solve_sycl.cpp
        pme_spread_sycl.cpp
      )
else()
    gmx_add_libgromacs_sources(
        # Files that implement stubs
        pme_gpu_program_impl.cpp
        )
endif()

# Source files have the following private module dependencies.
target_link_libraries(ewald PRIVATE
                      #                      gmxlib
                      #                      math
                      #                      mdtypes
                      #                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(ewald PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(ewald PUBLIC
target_link_libraries(ewald INTERFACE
                      legacy_api
                      )

# TODO: when fileio is an OBJECT target
#target_link_libraries(ewald PUBLIC legacy_api)
#target_link_libraries(ewald PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(ewald PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(ewald PRIVATE legacy_modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()


set(PME_OCL_KERNEL_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/pme_gpu_calculate_splines.clh"
    "${CMAKE_CURRENT_SOURCE_DIR}/pme_solve.clh"
    "${CMAKE_CURRENT_SOURCE_DIR}/pme_gather.clh"
    "${CMAKE_CURRENT_SOURCE_DIR}/pme_spread.clh")

if(CLANG_TIDY_EXE)
   set(OCL_COMPILER "${CLANG_TIDY_EXE}")
   set(CLANG_TIDY_ARGS "-quiet;-checks=*,-readability-implicit-bool-conversion,-llvm-header-guard,-hicpp-signed-bitwise,-clang-analyzer-deadcode.DeadStores,-google-readability-todo,-clang-diagnostic-padded,-fcomment-block-commands=internal;--;${CMAKE_C_COMPILER}")
else()
   set(OCL_COMPILER "${CMAKE_C_COMPILER}")
endif()

# TODO: test all warp sizes on all vendor targets?
foreach(VENDOR AMD NVIDIA INTEL APPLE)
    foreach(WARPSIZE 16 32 64)
        math(EXPR SPREAD_WG_SIZE "8*${WARPSIZE}")
        math(EXPR SOLVE_WG_SIZE "8*${WARPSIZE}")
        math(EXPR GATHER_WG_SIZE "4*${WARPSIZE}")
        set(OBJ_FILE pme_ocl_kernel_warpSize${WARPSIZE}_${VENDOR}.o)
        add_custom_command(OUTPUT ${OBJ_FILE} COMMAND ${OCL_COMPILER}
        ${CMAKE_CURRENT_SOURCE_DIR}/pme_program.cl ${CLANG_TIDY_ARGS}
        -Xclang -finclude-default-header  -D_${VENDOR}_SOURCE_
        -Dwarp_size=${WARPSIZE}
        -Dorder=4
        -DthreadsPerAtom=16
        -Dc_pmeMaxUnitcellShift=2
        -Dc_skipNeutralAtoms=false
        -Dc_virialAndEnergyCount=7
        -Dc_spreadWorkGroupSize=${SPREAD_WG_SIZE}
        -Dc_solveMaxWorkGroupSize=${SOLVE_WG_SIZE}
        -Dc_gatherWorkGroupSize=${GATHER_WG_SIZE}
        -DDIM=3 -DXX=0 -DYY=1 -DZZ=2
        -DwrapX=true -DwrapY=true
        -c -I ${CMAKE_SOURCE_DIR}/src -std=cl1.2
        -Weverything  -Wno-conversion -Wno-missing-variable-declarations -Wno-used-but-marked-unused
        -Wno-cast-align -Wno-incompatible-pointer-types
        # to avoid  "warning: unknown command tag name" for \internal
        -Wno-documentation-unknown-command
        # to avoid pme_gpu_types.h:100:52: warning: padding struct 'struct PmeGpuConstParams' with 4 bytes to align 'd_virialAndEnergy'
        -Wno-padded
        -o${OBJ_FILE}
        )
        list(APPEND PME_OCL_KERNELS ${OBJ_FILE})
    endforeach()
endforeach()

add_custom_target(ocl_pme_kernels DEPENDS ${PME_OCL_KERNELS} )
gmx_set_custom_target_output(ocl_pme_kernels ${PME_OCL_KERNELS})
