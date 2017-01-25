#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

function (gmx_test_clang_cuda_support)

    if ((NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR
        (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.9"))
        message(FATAL_ERROR "clang 3.9 or later required with GMX_CLANG_CUDA=ON!")
    endif()

    # NOTE: we'd ideally like to use a compile check here, but the link-stage
    # fails as the clang invocation generated seems to not handle well some
    # (GPU code) in the object file generated during compilation.
    # SET(CMAKE_REQUIRED_FLAGS ${FLAGS})
    # SET(CMAKE_REQUIRED_LIBRARIES ${LIBS})
    # CHECK_CXX_SOURCE_COMPILES("int main() { int c; cudaGetDeviceCount(&c); return 0; }" _CLANG_CUDA_COMPILES)
endfunction ()


if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.0" AND
    NOT CUDA_VERSION VERSION_LESS "8.0")
    message(FATAL_ERROR "clang ${CMAKE_CXX_COMPILER_VERSION} for CUDA is only compatible with CUDA version <8.0")
endif()

if (GMX_CUDA_TARGET_COMPUTE)
    message(WARNING "Values passed in GMX_CUDA_TARGET_COMPUTE will be ignored; clang will by default include PTX in the binary.")
endif()

if (GMX_CUDA_TARGET_SM)
    set(_CUDA_CLANG_GENCODE_FLAGS)
    set(_target_sm_list ${GMX_CUDA_TARGET_SM})
    foreach(_target ${_target_sm_list})
        list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_${_target}")
    endforeach()
else()
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_20")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_30")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_35")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_37")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_50")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_52")
    if (NOT CUDA_VERSION VERSION_LESS 8.0)
        list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_60")
        list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_61")
    endif()
    # TODO: test CUDA 9.0 and figure out which clang releases support it
    #       and the sm_70 arch.
endif()
if (GMX_CUDA_TARGET_SM)
    set_property(CACHE GMX_CUDA_TARGET_SM PROPERTY HELPSTRING "List of CUDA GPU architecture codes to compile for (without the sm_ prefix)")
    set_property(CACHE GMX_CUDA_TARGET_SM PROPERTY TYPE STRING)
endif()

# default flags
list(APPEND _CUDA_CLANG_FLAGS "-x cuda" "-ffast-math")
# CUDA toolkit
list(APPEND _CUDA_CLANG_FLAGS "--cuda-path=${CUDA_TOOLKIT_ROOT_DIR}")
# codegen flags
list(APPEND _CUDA_CLANG_FLAGS "${_CUDA_CLANG_GENCODE_FLAGS}")
foreach(_flag ${_CUDA_CLANG_FLAGS})
    set(GMX_CUDA_CLANG_FLAGS "${GMX_CUDA_CLANG_FLAGS} ${_flag}")
endforeach()

if (CUDA_USE_STATIC_CUDA_RUNTIME)
    set(GMX_CUDA_CLANG_LINK_LIBS "cudart_static")
else()
    set(GMX_CUDA_CLANG_LINK_LIBS "cudart")
endif()
set(GMX_CUDA_CLANG_LINK_LIBS "${GMX_CUDA_CLANG_LINK_LIBS}" "dl" "rt")
if (CUDA_64_BIT_DEVICE_CODE)
    set(GMX_CUDA_CLANG_LINK_DIRS "${CUDA_TOOLKIT_ROOT_DIR}/lib64")
else()
    set(GMX_CUDA_CLANG_LINK_DIRS "${CUDA_TOOLKIT_ROOT_DIR}/lib")
endif()

gmx_test_clang_cuda_support()
