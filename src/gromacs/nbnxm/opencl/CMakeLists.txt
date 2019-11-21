#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2018,2019, by the GROMACS development team, led by
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
    set(NBNXM_SOURCES ${NBNXM_SOURCES} ${OPENCL_NB_SOURCES} PARENT_SCOPE)
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
if(CLANG_TIDY_EXE)
   set(OCL_COMPILER "${CLANG_TIDY_EXE}")
   # TODO: remove readability-isolate-declaration when the nbnxm OpenCL kernels get modernized
   set(CLANG_TIDY_ARGS "-quiet;-checks=*,-readability-isolate-declaration,-readability-implicit-bool-conversion,-llvm-header-guard,-hicpp-signed-bitwise,-clang-analyzer-deadcode.DeadStores,-google-readability-todo;--;${CMAKE_C_COMPILER}")
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
            set(OBJ_FILE nbnxm_ocl_kernel${ELEC_NAME}${VDW_NAME}_${VENDOR}.o)
            add_custom_command(OUTPUT ${OBJ_FILE} COMMAND ${OCL_COMPILER}
                ${CMAKE_CURRENT_SOURCE_DIR}/nbnxm_ocl_kernels.cl ${CLANG_TIDY_ARGS}
                -Xclang -finclude-default-header  -D_${VENDOR}_SOURCE_
                -DGMX_OCL_FASTGEN ${ELEC_DEF} ${VDW_DEF}
                -DNBNXN_GPU_CLUSTER_SIZE=${CLUSTER_SIZE} -DIATYPE_SHMEM
                -c -I ${CMAKE_SOURCE_DIR}/src -std=cl1.2
                -Weverything  -Wno-conversion -Wno-missing-variable-declarations -Wno-used-but-marked-unused
                -Wno-cast-align -Wno-incompatible-pointer-types
                -o${OBJ_FILE}
                )
            list(APPEND NBNXM_OCL_KERNELS ${OBJ_FILE})
        endforeach()
    endforeach()
endforeach()
add_custom_target(ocl_nbnxm_kernels DEPENDS ${NBNXM_OCL_KERNELS})
