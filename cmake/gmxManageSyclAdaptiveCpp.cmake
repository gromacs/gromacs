#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
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

include(gmxFindFlagsForSource)

# This file has a lot of duplication in order to handle both hipSYCL 0.9.4 and future versions of
# AdaptiveCpp. It all can be simplified once backward compatibility is no loger required.
# See #4720

# Return all current CMake variables with name starting with "hipsycl" or "acpp" (case-insensitive).
# Result is in the form of a list of flags ("-Dfoo=bar;-Dbaz=true").
# Semicolons in values are escaped (needed for ACPP_TARGETS/HIPSYCL_TARGETS).
function(_getACppCmakeFlags RETURN_VAR)
    get_cmake_property(_VARS VARIABLES)
    list (SORT _VARS)
    set(RESULT "")
    foreach (_VARNAME ${_VARS})
        string(TOLOWER "${_VARNAME}" _VARNAME_LOWER)
        if (${_VARNAME_LOWER} MATCHES "^hipsycl" OR ${_VARNAME_LOWER} MATCHES "^acpp")
            # Escape semicolon. The number of backslashes was determined empirically.
            string(REPLACE ";" "\\\\\\;" _VARVALUE "${${_VARNAME}}")
            list(APPEND
                RESULT
                -D${_VARNAME}=${_VARVALUE}
                )
        endif()
    endforeach()
    set("${RETURN_VAR}" ${RESULT} PARENT_SCOPE)
endfunction()

if(NOT GMX_GPU_SYCL OR GMX_SYCL_DPCPP OR NOT GMX_SYCL_ACPP)
    message(FATAL_ERROR "Internal error: AdaptiveCpp/hipSYCL configuration script was included when it should not")
endif()

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    message(FATAL_ERROR "AdaptiveCpp/hipSYCL build requires Clang compiler, but ${CMAKE_CXX_COMPILER_ID} is used")
endif()
set(ACPP_CLANG "${CMAKE_CXX_COMPILER}")
set(HIPSYCL_CLANG "${ACPP_CLANG}")

# -Wno-unknown-cuda-version because Clang often complains about the newest CUDA, despite working fine with it.
# -Wno-unknown-attributes because AdaptiveCpp does not support reqd_sub_group_size (because it can only do some sub group sizes).
#    The latter can be added to ACPP_SYCLCC_EXTRA_COMPILE_OPTIONS
set(ACPP_EXTRA_ARGS "-Wno-unknown-cuda-version -Wno-unknown-attributes ${SYCL_CXX_FLAGS_EXTRA}")

# -ffast-math for performance
set(ACPP_EXTRA_COMPILE_OPTIONS -ffast-math)

# We want to inline aggressively, but only Clang 13 or newer supports this flag.
# Likely not needed on AMD, since AdaptiveCpp by default sets AMD-specific flags to force inlining, but no harm either.
check_cxx_compiler_flag("-fgpu-inline-threshold=1" HAS_GPU_INLINE_THRESHOLD)
if(${HAS_GPU_INLINE_THRESHOLD})
    list(APPEND ACPP_EXTRA_COMPILE_OPTIONS -fgpu-inline-threshold=99999)
endif()

# Backward-compatibility with hipSYCL 0.9.4
if (DEFINED HIPSYCL_TARGETS AND NOT DEFINED ACPP_TARGETS)
    set(ACPP_TARGETS "${HIPSYCL_TARGETS}")
endif()

# Must be called before find_package to capture all user-set CMake variables, but not those set automatically
_getACppCmakeFlags(_ALL_ACPP_CMAKE_FLAGS)

find_package(adaptivecpp QUIET)
if(NOT adaptivecpp_FOUND)
    set(HIPSYCL_SYCLCC_EXTRA_ARGS "${ACPP_EXTRA_ARGS}")
    set(HIPSYCL_SYCLCC_EXTRA_COMPILE_OPTIONS ${ACPP_EXTRA_COMPILE_OPTIONS})
    find_package(hipsycl REQUIRED)
    # Ensure the try_compile projects below find the same hipsycl)
    list(APPEND _ALL_ACPP_CMAKE_FLAGS -Dhipsycl_DIR=${hipsycl_DIR})
else()
    # Ensure the try_compile projects below find the same adaptivecpp)
    list(APPEND _ALL_ACPP_CMAKE_FLAGS -Dadaptivecpp_DIR=${adaptivecpp_DIR})
endif()

# If the user-set CMake variables change (e.g. because the user
# changed HIPSYCL_TARGETS/ACPP_TARGETS), then the try_compile tests below need
# to be re-run. Set and use an internal cache variable to detect
# the change and set a flag to rerun the tests.
if (DEFINED GMX_ALL_ACPP_CMAKE_FLAGS_COPY AND "${GMX_ALL_ACPP_CMAKE_FLAGS_COPY}" STREQUAL "${_ALL_ACPP_CMAKE_FLAGS}")
    set(_rerun_acpp_try_compile_tests FALSE)
else()
    # The new value should over-write the previous copy
    set(GMX_ALL_ACPP_CMAKE_FLAGS_COPY ${_ALL_ACPP_CMAKE_FLAGS} CACHE INTERNAL "Store the list of CMake variables needed for AdaptiveCpp/hipSYCL compilation test projects")
    set(_rerun_acpp_try_compile_tests TRUE)
endif()

# Does the AdaptiveCpp compiler work at all for the given targets?
if (NOT DEFINED GMX_ACPP_COMPILATION_WORKS OR _rerun_acpp_try_compile_tests)
    message(STATUS "Checking for valid AdaptiveCpp/hipSYCL compiler")
    try_compile(GMX_ACPP_COMPILATION_WORKS "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest" "${CMAKE_SOURCE_DIR}/cmake/AdaptiveCppTest/" "AdaptiveCppTest"
        OUTPUT_VARIABLE _ACPP_COMPILATION_OUTPUT
        CMAKE_FLAGS
        ${_ALL_ACPP_CMAKE_FLAGS})
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest")
    if(GMX_ACPP_COMPILATION_WORKS)
        message(STATUS "Checking for valid AdaptiveCpp/hipSYCL compiler - Success")
    endif()
endif()
if (NOT GMX_ACPP_COMPILATION_WORKS)
    message(FATAL_ERROR "AdaptiveCpp/hipSYCL compiler not working:\n${_ACPP_COMPILATION_OUTPUT}")
endif()

# Does AdaptiveCpp compilation target CUDA devices?
if(NOT DEFINED GMX_ACPP_HAVE_CUDA_TARGET OR _rerun_acpp_try_compile_tests)
    message(STATUS "Checking for AdaptiveCpp/hipSYCL CUDA target")
    try_compile(GMX_ACPP_HAVE_CUDA_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest" "${CMAKE_SOURCE_DIR}/cmake/AdaptiveCppTest/" "AdaptiveCppTest"
        CMAKE_FLAGS
        -DCHECK_CUDA_TARGET=ON
        ${_ALL_ACPP_CMAKE_FLAGS})
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest")
    if(GMX_ACPP_HAVE_CUDA_TARGET)
        message(STATUS "Checking for AdaptiveCpp/hipSYCL CUDA target - Enabled")
    else()
        message(STATUS "Checking for AdaptiveCpp/hipSYCL CUDA target - Disabled")
    endif()
endif()

# Does AdaptiveCpp compilation target HIP devices?
if(NOT DEFINED GMX_ACPP_HAVE_HIP_TARGET OR _rerun_acpp_try_compile_tests)
    message(STATUS "Checking for AdaptiveCpp/hipSYCL HIP target")
    try_compile(GMX_ACPP_HAVE_HIP_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest" "${CMAKE_SOURCE_DIR}/cmake/AdaptiveCppTest/" "AdaptiveCppTest"
        CMAKE_FLAGS
        -DCHECK_HIP_TARGET=ON
        ${_ALL_ACPP_CMAKE_FLAGS})
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest")
    if(GMX_ACPP_HAVE_HIP_TARGET)
        message(STATUS "Checking for AdaptiveCpp/hipSYCL HIP target - Enabled")
    else()
        message(STATUS "Checking for AdaptiveCpp/hipSYCL HIP target - Disabled")
    endif()
endif()

# Does AdaptiveCpp compilation target Intel LevelZero devices?
if(NOT DEFINED GMX_ACPP_HAVE_LEVELZERO_TARGET OR _rerun_acpp_try_compile_tests)
    message(STATUS "Checking for AdaptiveCpp/hipSYCL LevelZero target")
    try_compile(GMX_ACPP_HAVE_LEVELZERO_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest" "${CMAKE_SOURCE_DIR}/cmake/AdaptiveCppTest/" "AdaptiveCppTest"
        CMAKE_FLAGS
        -DCHECK_LEVELZERO_TARGET=ON
        ${_ALL_ACPP_CMAKE_FLAGS})
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest")
    if(GMX_ACPP_HAVE_LEVELZERO_TARGET)
        message(STATUS "Checking for AdaptiveCpp/hipSYCL LevelZero target - Enabled")
        message(FATAL_ERROR "GROMACS does not support LevelZero backend of AdaptiveCpp/hipSYCL")
    else()
        message(STATUS "Checking for AdaptiveCpp/hipSYCL LevelZero target - Disabled")
    endif()
endif()

# Does AdaptiveCpp compilation target generic (SSCP) compulation flow?
if(NOT DEFINED GMX_ACPP_HAVE_GENERIC_TARGET OR _rerun_acpp_try_compile_tests)
    message(STATUS "Checking for AdaptiveCpp/hipSYCL generic target")
    try_compile(GMX_ACPP_HAVE_GENERIC_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest" "${CMAKE_SOURCE_DIR}/cmake/AdaptiveCppTest/" "AdaptiveCppTest"
        CMAKE_FLAGS
        -DCHECK_GENERIC_TARGET=ON
        ${_ALL_ACPP_CMAKE_FLAGS})
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest")
    if(GMX_ACPP_HAVE_GENERIC_TARGET)
        message(STATUS "Checking for AdaptiveCpp/hipSYCL generic target - Enabled")
        message(FATAL_ERROR "GROMACS does not support generic/SSCP compilation flow of AdaptiveCpp/hipSYCL")
    else()
        message(STATUS "Checking for AdaptiveCpp/hipSYCL generic target - Disabled")
    endif()
endif()

if(NOT GMX_ACPP_HAVE_CUDA_TARGET AND NOT GMX_ACPP_HAVE_HIP_TARGET)
    message(WARNING "AdaptiveCpp/hipSYCL has no GPU targets set! Please, specify target hardware with -DHIPSYCL_TARGETS CMake option")
endif()
if(GMX_ACPP_HAVE_CUDA_TARGET AND GMX_ACPP_HAVE_HIP_TARGET)
    message(FATAL_ERROR "AdaptiveCpp/hipSYCL cannot have both CUDA and HIP targets active! This would require explicit multipass mode which both decreases performance on NVIDIA devices and has been removed in clang 12. Compile only for either CUDA or HIP targets.")
endif()
unset(_rerun_acpp_try_compile_tests)

if(GMX_GPU_FFT_VKFFT)
    include(gmxManageVkFft)
    if (GMX_ACPP_HAVE_CUDA_TARGET)
        gmx_manage_vkfft("CUDA")
    elseif (GMX_ACPP_HAVE_HIP_TARGET)
        gmx_manage_vkfft("HIP")
    else()
        message(FATAL_ERROR "VkFFT can only be used with HIP or CUDA backends")
    endif()
    set(_sycl_has_valid_fft TRUE)
endif()

# Try to detect if we need RDNA support. Not very robust, but should cover the most common use.
if (GMX_ACPP_HAVE_HIP_TARGET AND ${ACPP_TARGETS} MATCHES "gfx1[0-9][0-9][0-9]")
    set(_enable_rdna_support_automatically ON)
else()
    set(_enable_rdna_support_automatically OFF)
    # We assume that any GCN2-5 architecture (gfx7/8) and CDNA1-3 (gfx9 series) up until the time of writing of this conditional is 64-wide
    if (${ACPP_TARGETS} MATCHES "gfx[7-8][0-9][0-9]|gfx9[0-4][0-9ac]")
        option(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
            "Disable NBNXM GPU cluster pair splitting. Only supported with SYCL and 64-wide GPU architectures (like AMD GCN/CDNA)."
            ON)
        mark_as_advanced(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT)
    endif()
endif()
option(GMX_ACPP_ENABLE_AMD_RDNA_SUPPORT
    "Enable compiling kernels for AMD RDNA GPUs (gfx1xxx). When OFF, only CDNA and GCN are supported. Only used with AdaptiveCpp/hipSYCL."
    ${_enable_rdna_support_automatically})
mark_as_advanced(GMX_ACPP_ENABLE_AMD_RDNA_SUPPORT)

# Find a suitable rocFFT when AdaptiveCpp is targeting AMD devices
if (GMX_ACPP_HAVE_HIP_TARGET AND GMX_GPU_FFT_ROCFFT)
    # For consistency, we prefer to find rocFFT as part of the
    # default ROCm distribution that supports the version of
    # AdaptiveCpp that is being used. Other installations of rocFFT
    # might work, but could lead to problems that are hard to
    # trace.
    #
    # The AdaptiveCpp find package sets HIPSYCL_SYCLCC/ACPP_COMPILER which we can
    # use to find the JSON configuration file that points to the
    # default ROCm installation used by hipSYCL, which can be used
    # to find rocFFT.
    #
    # If this is unavailable or does not work, the user will need to
    # set CMAKE_PREFIX_PATH so CMake is able to find the dependencies
    # of rocFFT (namely hip, AMDDeviceLibs, amd_comgr, hsa-runtime64,
    # ROCclr).
    if (HIPSYCL_SYCLCC)
        get_filename_component(HIPSYCL_SYCLCC_DIR ${HIPSYCL_SYCLCC} DIRECTORY)
        find_file(ACPP_JSON syclcc.json
            HINTS ${HIPSYCL_SYCLCC_DIR}/../etc/hipSYCL
            DOC "location of hipSYCL/AdaptiveCpp JSON configuration file"
        )
    elseif(ACPP_COMPILER)
        get_filename_component(ACPP_COMPILER_DIR ${ACPP_COMPILER} DIRECTORY)
        # As of AdaptiveCpp 23.10.0, the file is still called "etc/hipSYCL/syclcc.json"
        find_file(ACPP_JSON syclcc.json
            HINTS ${ACPP_COMPILER_DIR}/../etc/hipSYCL
            DOC "location of hipSYCL/AdaptiveCpp JSON configuration file"
        )
    endif()
    if(HIPSYCL_SYCLCC_ROCM_PATH AND NOT ACPP_ROCM_PATH)
        set(ACPP_ROCM_PATH ${HIPSYCL_SYCLCC_ROCM_PATH})
    endif()
    if (ACPP_JSON AND NOT ACPP_ROCM_PATH)
        file(READ "${ACPP_JSON}" ACPP_JSON_CONTENTS)
        if (CMAKE_VERSION VERSION_LESS 3.19)
            # We want the value encoded by the line
            # "default-rocm-path" : "/opt/rocm",
            # so we use regular expressions to remove everything before
            # and after the relevant quotation marks.
            #
            # Remove this when GROMACS requires CMake 3.19 or higher, as the
            # proper JSON parsing below is more robust.
            string(REGEX REPLACE ".*\"default-rocm-path\" *: * \"" "" ACPP_ROCM_PATH_VALUE ${ACPP_JSON_CONTENTS})
            string(REGEX REPLACE "\",.*" "" ACPP_ROCM_PATH_VALUE ${HIPSYCL_SYCLCC_ROCM_PATH_VALUE})
        else()
            string(JSON ACPP_ROCM_PATH_VALUE GET ${ACPP_JSON_CONTENTS} "default-rocm-path")
        endif()
        set(ACPP_ROCM_PATH ${ACPP_ROCM_PATH_VALUE} CACHE FILEPATH "The default ROCm used by AdaptiveCpp/hipSYCL")
    endif()

    if(ACPP_ROCM_PATH)
        # Teach the rocFFT find package how to find the necessary components
        # from the ROCm distribution used by hipSYCL.
        set(hip_DIR ${ACPP_ROCM_PATH}/hip/lib/cmake/hip)
        set(AMDDeviceLibs_DIR ${ACPP_ROCM_PATH}/lib/cmake/AMDDeviceLibs)
        set(amd_comgr_DIR ${ACPP_ROCM_PATH}/lib/cmake/amd_comgr)
        set(hsa-runtime64_DIR ${ACPP_ROCM_PATH}/lib/cmake/hsa-runtime64)
        set(ROCclr_DIR ${ACPP_ROCM_PATH}/rocclr/lib/cmake/rocclr)
        set(rocfft_DIR ${ACPP_ROCM_PATH}/rocfft/lib/cmake/rocfft)
    endif()

    # Find rocFFT, either from the ROCm used by AdaptiveCpp, or as otherwise found on the system
    find_package(rocfft ${FIND_ROCFFT_QUIETLY} CONFIG HINTS ${ACPP_ROCM_PATH} PATHS /opt/rocm)
    if (NOT rocfft_FOUND)
        message(FATAL_ERROR "rocFFT is required for the AdaptiveCpp/hipSYCL build, but was not found")
    endif()
    set(FIND_ROCFFT_QUIETLY "QUIET")
    set(_sycl_has_valid_fft TRUE)
endif()

# Set new veriables for use in buildinfo.h.cmakein
if (hipsycl_FOUND)
    set(ACPP_COMPILER_LAUNCHER "${HIPSYCL_SYCLCC_LAUNCHER}")
    set(ACPP_EXTRA_ARGS "${HIPSYCL_SYCLCC_EXTRA_ARGS}")
    set(ACPP_EXTRA_COMPILE_OPTIONS "${HIPSYCL_SYCLCC_EXTRA_COMPILE_OPTIONS}")
    set(ACPP_TARGETS "${HIPSYCL_TARGETS}")
endif()

# Mark AdaptiveCpp-related CMake options as "advanced"
get_cmake_property(_VARS VARIABLES)
foreach (_VARNAME ${_VARS})
    if (_VARNAME MATCHES "^HIPSYCL" OR _VARNAME MATCHES "^ACPP")
        mark_as_advanced(${_VARNAME})
    endif()
endforeach()
mark_as_advanced(CLEAR HIPSYCL_TARGETS)
mark_as_advanced(CLEAR ACPP_TARGETS)
