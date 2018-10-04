#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2018, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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

if(GMX_USE_OPENCL)
    file(GLOB OPENCL_NB_SOURCES *.cpp)
    set(MDLIB_SOURCES ${MDLIB_SOURCES} ${OPENCL_NB_SOURCES} PARENT_SCOPE)
    file(GLOB MDLIB_OPENCL_KERNELS *.cl *.clh)
    set(MDLIB_OPENCL_KERNELS ${MDLIB_OPENCL_KERNELS} PARENT_SCOPE)
endif()

set(ELEC_DEFS
    "-DEL_CUTOFF\;-DEELNAME=_ElecCut"
    "-DEL_RF\;-DEELNAME=_ElecRF"
    "-DEL_EWALD_TAB\;-DEELNAME=_ElecEwQSTab"
    "-DEL_EWALD_TAB\;-DVDW_CUTOFF_CHECK\;-DEELNAME=_ElecEwQSTabTwinCut"
    "-DEL_EWALD_ANA\;-DEELNAME=_ElecEw"
    "-DEL_EWALD_ANA\;-DVDW_CUTOFF_CHECK\;-DEELNAME=_ElecEwTwinCut")
set(VDW_DEFS
    "-DVDWNAME=_VdwLJ"
    "-DLJ_COMB_GEOM\;-DVDWNAME=_VdwLJCombGeom"
    "-DLJ_COMB_LB\;-DVDWNAME=_VdwLJCombLB"
    "-DLJ_FORCE_SWITCH\;-DVDWNAME=_VdwLJFsw"
    "-DLJ_POT_SWITCH\;-DVDWNAME=_VdwLJPsw"
    "-DLJ_EWALD_COMB_GEOM\;-DVDWNAME=_VdwLJEwCombGeom"
    "-DLJ_EWALD_COMB_LB\;-DVDWNAME=_VdwLJEwCombLB")
set(EXTRA_ARGS "-std=cl1.2;-Xclang;-finclude-default-header;-c")
if(GMX_OCL_SPIRV)
   set(OCL_COMPILER "${OCL_COMPILER};-cc1")
   set(EXTRA_ARGS "-emit-spirv;-triple=spir-unknown-unknown")
   if (GMX_OCL_SPIRV_USEC) #No need for both C and C++. Can be removed after testing.
      set (EXTRA_ARGS "${EXTRA_ARGS};-cl-std=CL1.2;-include;opencl.h")
   else()
      set (EXTRA_ARGS "${EXTRA_ARGS};-cl-std=c++;-isystem${OCL_LIBCLCXX}")
   endif()
elseif(CLANG_TIDY_EXE)
   set(OCL_COMPILER "${CLANG_TIDY_EXE}")
   set(EXTRA_ARGS "-quiet;-checks=*,-readability-implicit-bool-conversion,-llvm-header-guard,-hicpp-signed-bitwise,-clang-analyzer-deadcode.DeadStores,-google-readability-todo;--;${CMAKE_C_COMPILER};${EXTRA_ARGS}")
else()
   set(OCL_COMPILER "${CMAKE_C_COMPILER}")
endif()
foreach(ELEC_DEF IN LISTS ELEC_DEFS)
    foreach(VDW_DEF IN LISTS VDW_DEFS)
        foreach(VENDOR AMD NVIDIA INTEL)
            if(VENDOR STREQUAL INTEL)
               set(CLUSTER_SIZE 4)
            else()
               set(CLUSTER_SIZE 8)
            endif()
            string(REGEX REPLACE ".*=" "" ELEC_NAME "${ELEC_DEF}")
            string(REGEX REPLACE ".*=" "" VDW_NAME "${VDW_DEF}")
            set(OBJ_FILE ${CMAKE_CURRENT_BINARY_DIR}/nbnxn_ocl_kernel${ELEC_NAME}${VDW_NAME}_${VENDOR}.o)
            add_custom_command(OUTPUT ${OBJ_FILE} COMMAND ${OCL_COMPILER}
                ${CMAKE_CURRENT_SOURCE_DIR}/nbnxn_ocl_kernels.cl ${EXTRA_ARGS}
                -D_${VENDOR}_SOURCE_
                -DGMX_OCL_FASTGEN ${ELEC_DEF} ${VDW_DEF}
                -DNBNXN_GPU_CLUSTER_SIZE=${CLUSTER_SIZE} -DIATYPE_SHMEM
                -I ${CMAKE_SOURCE_DIR}/src
                -Weverything  -Wno-conversion -Wno-missing-variable-declarations -Wno-used-but-marked-unused
                -Wno-cast-align -Wno-c++98-compat
                -Wno-old-style-cast #TODO: Remove if support for C is removed
                -o${OBJ_FILE}
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/nbnxn_ocl_kernels.cl #TODO rebuild if headers are changed
                )
            list(APPEND NBNXN_OCL_KERNELS ${OBJ_FILE})
        endforeach()
    endforeach()
endforeach()
add_custom_target(ocl_kernel DEPENDS ${NBNXN_OCL_KERNELS})
