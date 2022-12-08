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

function (gmx_test_clang_cuda_support)

    if ((NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR
        (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6))
        message(FATAL_ERROR "clang 6 or later required with GMX_CLANG_CUDA=ON!")
    endif()

    # NOTE: we'd ideally like to use a compile check here, but the link-stage
    # fails as the clang invocation generated seems to not handle well some
    # (GPU code) in the object file generated during compilation.
    # SET(CMAKE_REQUIRED_FLAGS ${FLAGS})
    # SET(CMAKE_REQUIRED_LIBRARIES ${LIBS})
    # CHECK_CXX_SOURCE_COMPILES("int main() { int c; cudaGetDeviceCount(&c); return 0; }" _CLANG_CUDA_COMPILES)
endfunction ()

if (GMX_CUDA_TARGET_COMPUTE)
    message(WARNING "Values passed in GMX_CUDA_TARGET_COMPUTE will be ignored; clang will by default include PTX in the binary.")
endif()

# At the time of writing, the latest released versions are Clang 15 and CUDA 11.8.
# Clang <14 support only CUDA 7.0-10.1; Clang 14-16 support CUDA 7.0-11.5.
# GROMACS requires CUDA 11.0, so no need to check for earlier versions
if ((CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0) OR (CUDA_VERSION VERSION_GREATER 11.5))
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 15.0)
        # We don't know about the future Clang versions, but so far Clang 16 docs state that only CUDA 7.0-11.5 are supported.
        set(_support_status "likely incompatible")
    elseif(CUDA_VERSION VERSION_GREATER 11.8)
        # No idea about future CUDA versions.
        set(_support_status "officially incompatible")
    else()
        # Our experience and multiple reports on the internet indicate that it works just fine.
        set(_support_status "officially incompatible (but generally working)")
    endif()
    message(NOTICE "Using ${_support_status} version of CUDA with Clang.")
    message(NOTICE "If Clang fails to recognize CUDA version, consider creating doing "
      "`echo \"CUDA Version ${CUDA_VERSION}\" | sudo tee \"${CUDA_TOOLKIT_ROOT_DIR}/version.txt\"`")
    list(APPEND _CUDA_CLANG_FLAGS "-Wno-unknown-cuda-version")
endif()

if (GMX_CUDA_TARGET_SM)
    set(_CUDA_CLANG_GENCODE_FLAGS)
    set(_target_sm_list ${GMX_CUDA_TARGET_SM})
    foreach(_target ${_target_sm_list})
        list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_${_target}")
    endforeach()
else()
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_35")
    # clang 6.0 + CUDA 9.0 seems to have issues generating code for sm_37
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0 OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0.999)
        list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_37")
    endif()
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_50")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_52")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_60")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_61")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_70")
    list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_75")
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0) # Clang 13 and earlier fail to recognize CUDA 11
        if(NOT CUDA_VERSION VERSION_LESS 11.0)
            list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_80")
        endif()
        if(NOT CUDA_VERSION VERSION_LESS 11.1)
            list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_86")
        endif()
    endif()
    # Enable once Clang recognizes sm_87 (https://reviews.llvm.org/D135306)
    # if(NOT CUDA_VERSION VERSION_LESS 11.4)
    #       list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_87")
    #    endif()
    # Enable once Clang recognizes CUDA 11.8 and newer (https://reviews.llvm.org/D135306)
    # if(NOT CUDA_VERSION VERSION_LESS 11.8)
    #     list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_89")
    # endif()
    # if(NOT CUDA_VERSION VERSION_LESS 12.0)
    #     list(APPEND _CUDA_CLANG_GENCODE_FLAGS "--cuda-gpu-arch=sm_90")
    # endif()
endif()
if (GMX_CUDA_TARGET_SM)
    set_property(CACHE GMX_CUDA_TARGET_SM PROPERTY HELPSTRING "List of CUDA GPU architecture codes to compile for (without the sm_ prefix)")
    set_property(CACHE GMX_CUDA_TARGET_SM PROPERTY TYPE STRING)
endif()

# default flags
list(APPEND _CUDA_CLANG_FLAGS "-x cuda" "-ffast-math" "-fcuda-flush-denormals-to-zero")
# Workaround for clang>=9 (Bug 45533). No CUDA file uses OpenMP.
list(APPEND _CUDA_CLANG_FLAGS "-fno-openmp")
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
