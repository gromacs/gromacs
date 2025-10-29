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

# Return all current CMake variables with name starting with "acpp" (case-insensitive).
# Result is in the form of a list of flags ("-Dfoo=bar;-Dbaz=true").
# Semicolons in values are escaped (needed for ACPP_TARGETS).
function(_getACppCmakeFlags RETURN_VAR)
    get_cmake_property(_VARS VARIABLES)
    list (SORT _VARS)
    set(RESULT "")
    foreach (_VARNAME ${_VARS})
        string(TOLOWER "${_VARNAME}" _VARNAME_LOWER)
        if (${_VARNAME_LOWER} MATCHES "^acpp")
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
    message(FATAL_ERROR "Internal error: AdaptiveCpp configuration script was included when it should not")
endif()

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    message(FATAL_ERROR "AdaptiveCpp build requires Clang compiler, but ${CMAKE_CXX_COMPILER_ID} is used")
endif()
set(ACPP_CLANG "${CMAKE_CXX_COMPILER}")

# -Wno-unknown-cuda-version because Clang often complains about the newest CUDA, despite working fine with it.
# -Wno-unknown-attributes because AdaptiveCpp does not support reqd_sub_group_size (because it can only do some sub group sizes).
#    The latter can be added to ACPP_SYCLCC_EXTRA_COMPILE_OPTIONS
set(ACPP_EXTRA_ARGS "-Wno-unknown-cuda-version -Wno-unknown-attributes ${SYCL_CXX_FLAGS_EXTRA}")

# -ffast-math for performance
set(ACPP_EXTRA_COMPILE_OPTIONS -ffast-math)

# Enable instant submission unless _ALLOW_INSTANT_SUBMISSION or _FORCE_INSTANT_SUBMISSION is already set (to 0 or 1)
# Not necessary since https://github.com/AdaptiveCpp/AdaptiveCpp/pull/179, so should be possible to remove this once we require ACpp 25.10
if (NOT SYCL_CXX_FLAGS_EXTRA MATCHES "_INSTANT_SUBMISSION")
    list(APPEND ACPP_EXTRA_COMPILE_OPTIONS -DHIPSYCL_ALLOW_INSTANT_SUBMISSION=1) # ACpp 24.02 and earlier
    list(APPEND ACPP_EXTRA_COMPILE_OPTIONS -DACPP_ALLOW_INSTANT_SUBMISSION=1) # ACpp 24.06 and newer
endif()

# We want to inline aggressively where the compiler supports this flag.
# Likely not needed on AMD, since AdaptiveCpp by default sets AMD-specific flags to force inlining, but no harm either.
check_cxx_compiler_flag("-fgpu-inline-threshold=1" HAS_GPU_INLINE_THRESHOLD)
if(${HAS_GPU_INLINE_THRESHOLD})
    list(APPEND ACPP_EXTRA_COMPILE_OPTIONS -fgpu-inline-threshold=99999)
endif()

# We cannot use local_accessor::get_multi_ptr with ACpp (https://github.com/AdaptiveCpp/AdaptiveCpp/issues/1230),
# but `get_pointer` is already marked deprecated since ACpp 24.06
check_cxx_compiler_flag("-Wno-deprecated-declarations" CXXFLAGS_NO_DEPRECATED_DECLARATIONS)
if(${CXXFLAGS_NO_DEPRECATED_DECLARATIONS})
    list(APPEND ACPP_EXTRA_COMPILE_OPTIONS -Wno-deprecated-declarations)
endif()

# Must be called before find_package to capture all user-set CMake variables, but not those set automatically
_getACppCmakeFlags(_ALL_ACPP_CMAKE_FLAGS)

if (adaptivecpp_FIND_QUETLY_AFTER_FIRST_RUN)
    set (adaptivecpp_FIND_QUETLY ON)
else ()
    set (adaptivecpp_FIND_QUETLY OFF)
endif ()
find_package(adaptivecpp REQUIRED)
set (adaptivecpp_FIND_QUETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to find AdaptiveCpp")
# Ensure the try_compile projects below find the same adaptivecpp
list(APPEND _ALL_ACPP_CMAKE_FLAGS -Dadaptivecpp_DIR=${adaptivecpp_DIR})

# If the user-set CMake variables change (e.g. because the user changed ACPP_TARGETS), then the try_compile tests
# below need to be re-run. Set and use an internal cache variable to detect the change and set a flag to rerun the tests.
string(SHA1 _ALL_ACPP_CMAKE_FLAGS_HASH "${_ALL_ACPP_CMAKE_FLAGS}")
if (DEFINED GMX_ALL_ACPP_CMAKE_FLAGS_HASH AND "${GMX_ALL_ACPP_CMAKE_FLAGS_HASH}" STREQUAL "${_ALL_ACPP_CMAKE_FLAGS_HASH}")
    set(_rerun_acpp_try_compile_tests FALSE)
else()
    # The new value should over-write the previous copy
    set(GMX_ALL_ACPP_CMAKE_FLAGS_HASH ${_ALL_ACPP_CMAKE_FLAGS_HASH} CACHE INTERNAL "Store the hash of list of CMake variables needed for AdaptiveCpp compilation test projects")
    set(_rerun_acpp_try_compile_tests TRUE)
endif()

# Does the AdaptiveCpp compiler work at all for the given targets?
if (NOT DEFINED GMX_ACPP_COMPILATION_WORKS OR _rerun_acpp_try_compile_tests)
    if (NOT ACPP_CHECK_QUIETLY_AFTER_FIRST_RUN)
        message(STATUS "Checking for valid AdaptiveCpp compiler")
    endif()
    try_compile(GMX_ACPP_COMPILATION_WORKS "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest" "${CMAKE_SOURCE_DIR}/cmake/AdaptiveCppTest/" "AdaptiveCppTest"
        OUTPUT_VARIABLE _ACPP_COMPILATION_OUTPUT
        CMAKE_FLAGS ${_ALL_ACPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpAdaptiveCppTest")
    if(GMX_ACPP_COMPILATION_WORKS)
        if (NOT ACPP_CHECK_QUIETLY_AFTER_FIRST_RUN)
            message(STATUS "Checking for valid AdaptiveCpp compiler - Success")
        endif()
        foreach (target IN ITEMS CUDA HIP HIP_WAVE32 HIP_WAVE64 SPIRV GENERIC)
            if (_ACPP_COMPILATION_OUTPUT MATCHES "GMX_SYCL_TEST_HAVE_${target}_TARGET")
                set(GMX_ACPP_HAVE_${target}_TARGET ON CACHE INTERNAL "AdaptiveCpp flags/configuration have ${target} target")
            else()
                set(GMX_ACPP_HAVE_${target}_TARGET OFF CACHE INTERNAL "AdaptiveCpp flags/configuration have ${target} target")
            endif()
            if (NOT ACPP_CHECK_QUIETLY_AFTER_FIRST_RUN)
                message(STATUS "AdaptiveCpp has ${target} target enabled: ${GMX_ACPP_HAVE_${target}_TARGET}")
            endif()
        endforeach()
    endif()
endif()
set(ACPP_CHECK_QUIETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to check ACpp compilation")

if (NOT GMX_ACPP_COMPILATION_WORKS)
    message(FATAL_ERROR "AdaptiveCpp compiler not working:\n${_ACPP_COMPILATION_OUTPUT}")
endif()
if (GMX_ACPP_HAVE_SPIRV_TARGET)
    message(FATAL_ERROR "GROMACS does not support LevelZero (SPIR-V) backend of AdaptiveCpp")
endif()
if(GMX_ACPP_HAVE_GENERIC_TARGET)
    message(WARNING "The generic/SSCP compilation flow of AdaptiveCpp is experimental")
endif()
if(NOT GMX_ACPP_HAVE_CUDA_TARGET AND NOT GMX_ACPP_HAVE_HIP_TARGET AND NOT GMX_ACPP_HAVE_GENERIC_TARGET)
    message(WARNING "AdaptiveCpp has no GPU targets set! Please, specify target hardware with -DACPP_TARGETS CMake option")
endif()
if(GMX_ACPP_HAVE_CUDA_TARGET AND GMX_ACPP_HAVE_HIP_TARGET)
    message(FATAL_ERROR "AdaptiveCpp cannot have both CUDA and HIP targets active! This would require explicit multipass mode which both decreases performance on NVIDIA devices and is no longer supported in clang. Compile only for either CUDA or HIP targets.")
endif()
if(GMX_ACPP_HAVE_HIP_TARGET AND NOT GMX_ACPP_HAVE_HIP_WAVE64_TARGET AND NOT GMX_ACPP_HAVE_HIP_WAVE32_TARGET)
    message(FATAL_ERROR "AdaptiveCpp has a HIP target, but cannot determine its wave size")
endif()
unset(_rerun_acpp_try_compile_tests)

if(GMX_GPU_FFT_VKFFT)
    include(gmxManageVkFft)
    set(_sycl_has_valid_fft TRUE)
endif()

if (GMX_ACPP_HAVE_HIP_TARGET AND GMX_ACPP_HAVE_HIP_WAVE32_TARGET)
    set(_enable_rdna_support_automatically ON)
else()
    set(_enable_rdna_support_automatically OFF)
    if (GMX_ACPP_HAVE_HIP_WAVE64_TARGET)
        option(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
            "Disable NBNXM GPU cluster pair splitting. Only supported with SYCL and 64-wide GPU architectures (like AMD GCN/CDNA)."
            ON)
        mark_as_advanced(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT)
    endif()
endif()
option(GMX_ENABLE_AMD_RDNA_SUPPORT
    "Enable compiling kernels for AMD RDNA GPUs (gfx1xxx). When OFF, only CDNA and GCN are supported. Only used with AdaptiveCpp."
    ${_enable_rdna_support_automatically})
mark_as_advanced(GMX_ENABLE_AMD_RDNA_SUPPORT)

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
    # default ROCm installation used by AdaptiveCpp, which can be used
    # to find rocFFT.
    #
    # If this is unavailable or does not work, the user will need to
    # set CMAKE_PREFIX_PATH so CMake is able to find the dependencies
    # of rocFFT (namely hip, AMDDeviceLibs, amd_comgr, hsa-runtime64,
    # ROCclr).
    get_filename_component(ACPP_COMPILER_DIR ${ACPP_COMPILER} DIRECTORY)
    find_file(ACPP_JSON
        NAMES acpp-rocm.json
        HINTS ${ACPP_COMPILER_DIR}/../etc/AdaptiveCpp ${ACPP_COMPILER_DIR}/../etc/hipSYCL
        DOC "location of AdaptiveCpp JSON configuration file"
    )
    if (ACPP_JSON AND NOT ACPP_ROCM_PATH)
        file(READ "${ACPP_JSON}" ACPP_JSON_CONTENTS)
        string(JSON ACPP_ROCM_PATH_VALUE GET ${ACPP_JSON_CONTENTS} "default-rocm-path")
        set(ACPP_ROCM_PATH "${ACPP_ROCM_PATH_VALUE}" CACHE PATH "The default ROCm used by AdaptiveCpp" FORCE)
    endif()

    if(ACPP_ROCM_PATH)
        # Teach the rocFFT find package how to find the necessary components
        # from the ROCm distribution used by AdaptiveCpp.
        set(hip_DIR ${ACPP_ROCM_PATH}/lib/cmake/hip)
        set(AMDDeviceLibs_DIR ${ACPP_ROCM_PATH}/lib/cmake/AMDDeviceLibs)
        set(amd_comgr_DIR ${ACPP_ROCM_PATH}/lib/cmake/amd_comgr)
        set(hsa-runtime64_DIR ${ACPP_ROCM_PATH}/lib/cmake/hsa-runtime64)
        set(ROCclr_DIR ${ACPP_ROCM_PATH}/lib/cmake/rocclr)
        set(rocfft_DIR ${ACPP_ROCM_PATH}/lib/cmake/rocfft)
    endif()

    # Find rocFFT, either from the ROCm used by AdaptiveCpp, or as otherwise found on the system
    find_package(rocfft ${FIND_ROCFFT_QUIETLY} CONFIG HINTS ${ACPP_ROCM_PATH} PATHS /opt/rocm)
    if (NOT rocfft_FOUND)
        message(FATAL_ERROR "rocFFT is requested but was not found")
    endif()
    set(FIND_ROCFFT_QUIETLY "QUIET")
    set(_sycl_has_valid_fft TRUE)
endif()

string(REPLACE ";" " " ACPP_EXTRA_COMPILE_OPTIONS_STR "${ACPP_EXTRA_COMPILE_OPTIONS}")
string(STRIP "${ACPP_EXTRA_COMPILE_OPTIONS_STR}" ACPP_EXTRA_COMPILE_OPTIONS_STR)

# Mark AdaptiveCpp-related CMake options as "advanced"
get_cmake_property(_VARS VARIABLES)
foreach (_VARNAME ${_VARS})
    if (_VARNAME MATCHES "^ACPP")
        mark_as_advanced(${_VARNAME})
    endif()
endforeach()
mark_as_advanced(CLEAR ACPP_TARGETS)
