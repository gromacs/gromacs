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

# OpenCL required version: 1.2 or newer
set(REQUIRED_SYCL_MIN_VERSION_MAJOR 1)
set(REQUIRED_SYCL_MIN_VERSION_MINOR 2)
set(REQUIRED_SYCL_MIN_VERSION ${REQUIRED_SYCL_MIN_VERSION_MAJOR}.${REQUIRED_SYCL_MIN_VERSION_MINOR})

set(GMX_GPU_SYCL ON)

# CMake issue tracking the efforts to make a universal upstream module:
# https://gitlab.kitware.com/cmake/cmake/-/issues/21711

option(GMX_SYCL_HIPSYCL "Use hipSYCL instead of Intel oneAPI for SYCL compilation" OFF)

if(GMX_DOUBLE)
    message(FATAL_ERROR "SYCL acceleration is not available in double precision")
endif()

include(gmxFindFlagsForSource)

set(_sycl_has_valid_fft FALSE)

# Return all current CMake variables with name starting with "hipsycl" (case-insensitive).
# Result is in the form of a list of flags ("-Dfoo=bar;-Dbaz=true").
# Semicolons in values are escaped (needed for HIPSYCL_TARGETS).
function(_getHipSyclCmakeFlags RETURN_VAR)
    get_cmake_property(_VARS VARIABLES)
    list (SORT _VARS)
    set(RESULT "")
    foreach (_VARNAME ${_VARS})
            string(TOLOWER "${_VARNAME}" _VARNAME_LOWER)
            if (${_VARNAME_LOWER} MATCHES "^hipsycl")
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

if(GMX_SYCL_HIPSYCL)
    if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        message(FATAL_ERROR "HipSYCL build requires Clang compiler, but ${CMAKE_CXX_COMPILER_ID} is used")
    endif()
    set(HIPSYCL_CLANG "${CMAKE_CXX_COMPILER}")

    # -Wno-unknown-cuda-version because Clang often complains about the newest CUDA, despite working fine with it.
    # -Wno-unknown-attributes because hipSYCL does not support reqd_sub_group_size (because it can only do some sub group sizes).
    #    The latter can be added to HIPSYCL_SYCLCC_EXTRA_COMPILE_OPTIONS
    set(HIPSYCL_SYCLCC_EXTRA_ARGS "-Wno-unknown-cuda-version -Wno-unknown-attributes ${SYCL_CXX_FLAGS_EXTRA}")

    # -ffast-math for performance
    set(HIPSYCL_SYCLCC_EXTRA_COMPILE_OPTIONS -ffast-math)

    # We want to inline aggressively, but only Clang 13 or newer supports this flag.
    # Likely not needed on AMD, since hipSYCL by default sets AMD-specific flags to force inlining, but no harm either.
    check_cxx_compiler_flag("-fgpu-inline-threshold=1" HAS_GPU_INLINE_THRESHOLD)
    if(${HAS_GPU_INLINE_THRESHOLD})
        list(APPEND HIPSYCL_SYCLCC_EXTRA_COMPILE_OPTIONS -fgpu-inline-threshold=99999)
    endif()

    # Must be called before find_package to capture all user-set CMake variables, but not those set automatically
    _getHipSyclCmakeFlags(_ALL_HIPSYCL_CMAKE_FLAGS)

    find_package(hipsycl REQUIRED)
    # Ensure the try_compile projects below find the same hipsycl)
    list(APPEND _ALL_HIPSYCL_CMAKE_FLAGS -Dhipsycl_DIR=${hipsycl_DIR})

    # If the user-set CMake variables change (e.g. because the user
    # changed HIPSYCL_TARGETS), then the try_compile tests below need
    # to be re-run. Set and use an internal cache variable to detect
    # the change and set a flag to rerun the tests.
    if (DEFINED GMX_ALL_HIPSYCL_CMAKE_FLAGS_COPY AND "${GMX_ALL_HIPSYCL_CMAKE_FLAGS_COPY}" STREQUAL "${_ALL_HIPSYCL_CMAKE_FLAGS}")
        set(_rerun_hipsycl_try_compile_tests FALSE)
    else()
        # The new value should over-write the previous copy
        set(GMX_ALL_HIPSYCL_CMAKE_FLAGS_COPY ${_ALL_HIPSYCL_CMAKE_FLAGS} CACHE INTERNAL "Store the list of CMake variables needed for hipSYCL compilation test projects")
        set(_rerun_hipsycl_try_compile_tests TRUE)
    endif()

    # Does the hipSYCL compiler work at all for the given targets?
    if (NOT DEFINED GMX_HIPSYCL_COMPILATION_WORKS OR _rerun_hipsycl_try_compile_tests)
        message(STATUS "Checking for valid hipSYCL compiler")
        try_compile(GMX_HIPSYCL_COMPILATION_WORKS "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest" "${CMAKE_SOURCE_DIR}/cmake/HipSyclTest/" "HipSyclTest"
          OUTPUT_VARIABLE _HIPSYCL_COMPILATION_OUTPUT
          CMAKE_FLAGS
            ${_ALL_HIPSYCL_CMAKE_FLAGS})
        file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest")
        if(GMX_HIPSYCL_COMPILATION_WORKS)
            message(STATUS "Checking for valid hipSYCL compiler - Success")
        endif()
    endif()
    if (NOT GMX_HIPSYCL_COMPILATION_WORKS)
        message(FATAL_ERROR "hipSYCL compiler not working:\n${_HIPSYCL_COMPILATION_OUTPUT}")
    endif()

    # Does hipSYCL compilation target CUDA devices?
    if(NOT DEFINED GMX_HIPSYCL_HAVE_CUDA_TARGET OR _rerun_hipsycl_try_compile_tests)
        message(STATUS "Checking for hipSYCL CUDA target")
        try_compile(GMX_HIPSYCL_HAVE_CUDA_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest" "${CMAKE_SOURCE_DIR}/cmake/HipSyclTest/" "HipSyclTest"
          CMAKE_FLAGS
            -DCHECK_CUDA_TARGET=ON
            ${_ALL_HIPSYCL_CMAKE_FLAGS})
        file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest")
        if(GMX_HIPSYCL_HAVE_CUDA_TARGET)
            message(STATUS "Checking for hipSYCL CUDA target - Success")
        else()
            message(STATUS "Checking for hipSYCL CUDA target - Failed")
        endif()
    endif()

    # Does hipSYCL compilation target HIP devices?
    if(NOT DEFINED GMX_HIPSYCL_HAVE_HIP_TARGET OR _rerun_hipsycl_try_compile_tests)
        message(STATUS "Checking for hipSYCL HIP target")
        try_compile(GMX_HIPSYCL_HAVE_HIP_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest" "${CMAKE_SOURCE_DIR}/cmake/HipSyclTest/" "HipSyclTest"
          CMAKE_FLAGS
            -DCHECK_HIP_TARGET=ON
            ${_ALL_HIPSYCL_CMAKE_FLAGS})
        file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest")
        if(GMX_HIPSYCL_HAVE_HIP_TARGET)
            message(STATUS "Checking for hipSYCL HIP target - Success")
        else()
            message(STATUS "Checking for hipSYCL HIP target - Failed")
        endif()
    endif()

    # Does hipSYCL compilation target Intel Level0 devices?
    if(NOT DEFINED GMX_HIPSYCL_HAVE_LEVELZERO_TARGET OR _rerun_hipsycl_try_compile_tests)
        message(STATUS "Checking for hipSYCL LevelZero target")
        try_compile(GMX_HIPSYCL_HAVE_LEVELZERO_TARGET "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest" "${CMAKE_SOURCE_DIR}/cmake/HipSyclTest/" "HipSyclTest"
          CMAKE_FLAGS
            -DCHECK_LEVELZERO_TARGET=ON
            ${_ALL_HIPSYCL_CMAKE_FLAGS})
        file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeTmpHipSyclTest")
        if(GMX_HIPSYCL_HAVE_LEVELZERO_TARGET)
            message(STATUS "Checking for hipSYCL LevelZero target - Success")
            message(WARNING "GROMACS does not support LevelZero backend of hipSYCL")
        else()
            message(STATUS "Checking for hipSYCL LevelZero target - Failed")
        endif()
    endif()

    if(NOT GMX_HIPSYCL_HAVE_CUDA_TARGET AND NOT GMX_HIPSYCL_HAVE_HIP_TARGET)
        message(WARNING "hipSYCL has no GPU targets set! Please, specify target hardware with -DHIPSYCL_TARGETS CMake option")
    endif()
    if(GMX_HIPSYCL_HAVE_CUDA_TARGET AND GMX_HIPSYCL_HAVE_HIP_TARGET)
        message(FATAL_ERROR "hipSYCL cannot have both CUDA and HIP targets active! This would require explicit multipass mode which both decreases performance on NVIDIA devices and has been removed in clang 12. Compile only for either CUDA or HIP targets.")
    endif()
    unset(_rerun_hipsycl_try_compile_tests)

    if(GMX_GPU_FFT_VKFFT)
        include(gmxManageVkFft)
        if (GMX_HIPSYCL_HAVE_CUDA_TARGET)
            gmx_manage_vkfft("CUDA")
        elseif (GMX_HIPSYCL_HAVE_HIP_TARGET)
            gmx_manage_vkfft("HIP")
        else()
            message(FATAL_ERROR "VkFFT can only be used with HIP or CUDA backends")
        endif()
        set(_sycl_has_valid_fft TRUE)
    endif()

    # Try to detect if we need RDNA support. Not very robust, but should cover the most common use.
    if (GMX_HIPSYCL_HAVE_HIP_TARGET AND ${HIPSYCL_TARGETS} MATCHES "gfx1[0-9][0-9][0-9]")
        set(_enable_rdna_support_automatically ON)
    else()
        set(_enable_rdna_support_automatically OFF)
        # We assume that any GCN2-5 architecture (gfx7/8) and CDNA1-3 (gfx9 series) up until the time of writing of this conditional is 64-wide
        if (${HIPSYCL_TARGETS} MATCHES "gfx[7-8][0-9][0-9]|gfx9[0-4][0-9ac]")
            option(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
                "Disable NBNXM GPU cluster pair splitting. Only supported with SYCL and 64-wide GPU architectures (like AMD GCN/CDNA)."
                ON)
            mark_as_advanced(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT)
        endif()
    endif()
    option(GMX_HIPSYCL_ENABLE_AMD_RDNA_SUPPORT
        "Enable compiling kernels for AMD RDNA GPUs (gfx1xxx). When OFF, only CDNA and GCN are supported. Only used with hipSYCL."
        ${_enable_rdna_support_automatically})
    mark_as_advanced(GMX_HIPSYCL_ENABLE_AMD_RDNA_SUPPORT)

    # Find a suitable rocFFT when hipSYCL is targeting AMD devices
    if (GMX_HIPSYCL_HAVE_HIP_TARGET AND GMX_GPU_FFT_ROCFFT)
        # For consistency, we prefer to find rocFFT as part of the
        # default ROCm distribution that supports the version of
        # hipSYCL that is being used. Other installations of rocFFT
        # might work, but could lead to problems that are hard to
        # trace.
        #
        # The hipSYCL find package sets HIPSYCL_SYCLCC which we can
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
            find_file(HIPSYCL_SYCLCC_JSON syclcc.json
                HINTS ${HIPSYCL_SYCLCC_DIR}/../etc/hipSYCL
                DOC "location of hipSYCL JSON configuration file"
                )
            if (HIPSYCL_SYCLCC_JSON)
                if(NOT HIPSYCL_SYCLCC_ROCM_PATH)
                    file(READ ${HIPSYCL_SYCLCC_JSON} HIPSYCL_SYCLCC_JSON_CONTENTS)
                    if (CMAKE_VERSION VERSION_LESS 3.19)
                        # We want the value encoded by the line
                        # "default-rocm-path" : "/opt/rocm",
                        # so we use regular expressions to remove everything before
                        # and after the relevant quotation marks.
                        #
                        # Remove this when GROMACS requires CMake 3.19 or higher, as the
                        # proper JSON parsing below is more robust.
                        string(REGEX REPLACE ".*\"default-rocm-path\" *: * \"" "" HIPSYCL_SYCLCC_ROCM_PATH_VALUE ${HIPSYCL_SYCLCC_JSON_CONTENTS})
                        string(REGEX REPLACE "\",.*" "" HIPSYCL_SYCLCC_ROCM_PATH_VALUE ${HIPSYCL_SYCLCC_ROCM_PATH_VALUE})
                    else()
                        string(JSON HIPSYCL_SYCLCC_ROCM_PATH_VALUE GET ${HIPSYCL_SYCLCC_JSON_CONTENTS} "default-rocm-path")
                    endif()
                    set(HIPSYCL_SYCLCC_ROCM_PATH ${HIPSYCL_SYCLCC_ROCM_PATH_VALUE} CACHE FILEPATH "The default ROCm used by syclcc from hipSYCL")
                endif()

                if(HIPSYCL_SYCLCC_ROCM_PATH)
                    # Teach the rocFFT find package how to find the necessary components
                    # from the ROCm distribution used by hipSYCL.
                    set(hip_DIR ${HIPSYCL_SYCLCC_ROCM_PATH}/hip/lib/cmake/hip)
                    set(AMDDeviceLibs_DIR ${HIPSYCL_SYCLCC_ROCM_PATH}/lib/cmake/AMDDeviceLibs)
                    set(amd_comgr_DIR ${HIPSYCL_SYCLCC_ROCM_PATH}/lib/cmake/amd_comgr)
                    set(hsa-runtime64_DIR ${HIPSYCL_SYCLCC_ROCM_PATH}/lib/cmake/hsa-runtime64)
                    set(ROCclr_DIR ${HIPSYCL_SYCLCC_ROCM_PATH}/rocclr/lib/cmake/rocclr)
                    set(rocfft_DIR ${HIPSYCL_SYCLCC_ROCM_PATH}/rocfft/lib/cmake/rocfft)
                endif()
            endif()
        endif()


        # Find rocFFT, either from the ROCm used by hipSYCL, or as otherwise found on the system
        find_package(rocfft ${FIND_ROCFFT_QUIETLY} CONFIG HINTS ${HIPSYCL_SYCLCC_ROCM_PATH} PATHS /opt/rocm)
        if (NOT rocfft_FOUND)
            message(FATAL_ERROR "rocFFT is required for the hipSYCL build, but was not found")
        endif()
        set(FIND_ROCFFT_QUIETLY "QUIET")
        set(_sycl_has_valid_fft TRUE)
    endif()

    # Mark hipsycl-related CMake options as "advanced"
    get_cmake_property(_VARS VARIABLES)
    foreach (_VARNAME ${_VARS})
        if (_VARNAME MATCHES "^HIPSYCL")
            mark_as_advanced(${_VARNAME})
        endif()
    endforeach()
    mark_as_advanced(CLEAR HIPSYCL_TARGETS)
else()
    if(WIN32)
        if(CMAKE_VERSION VERSION_LESS "3.23.0")
            message(FATAL_ERROR "SYCL with DPC++ on Windows requires cmake 3.23 or later.")
        endif()
        if(NOT BUILD_SHARED_LIBS)
            message(FATAL_ERROR "SYCL with DPC++ on Windows doesn't work with static libraries. Set BUILD_SHARED_LIBS=on.")
            # Tested up to 3.23.1 and icx 2022.1. Problem is order of exe link argument order. Works if gromacs.lib
            # and -fsycl both appear before -link. Not possible to change order from cmake script. cmake fix is WIP.
        endif()
    endif()
    if(CMAKE_CXX_COMPILER MATCHES "dpcpp")
        message(FATAL_ERROR "Intel's \"dpcpp\" compiler is deprecated; please use \"icpx\" for SYCL builds")
    endif()

    # Find the flags to enable (or re-enable) SYCL with Intel extensions. In case we turned it off above,
    # it's important that we check the combination of both flags, to make sure the second one re-enables SYCL.
    if(NOT CHECK_SYCL_CXX_FLAGS_QUIETLY)
        message(STATUS "Checking for flags to enable SYCL")
    endif()
    set(SAMPLE_SYCL_SOURCE
        "#include <CL/sycl.hpp>
         int main(){
             sycl::queue q(sycl::default_selector{});
             return 0;
         }")
    set(SYCL_CXX_FLAGS "-fsycl")
    gmx_check_source_compiles_with_flags(
        "${SAMPLE_SYCL_SOURCE}"
        "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_CXX_FLAGS}"
        "CXX"
        SYCL_CXX_FLAGS_RESULT
        )
    if (SYCL_CXX_FLAGS_RESULT)
        if(NOT CHECK_SYCL_CXX_FLAGS_QUIETLY)
            message(STATUS "Checking for flags to enable SYCL - ${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_CXX_FLAGS}")
        endif()
        set(CHECK_SYCL_CXX_FLAGS_QUIETLY 1 CACHE INTERNAL "Keep quiet on future calls to detect SYCL flags" FORCE)
        set(SYCL_TOOLCHAIN_CXX_FLAGS "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_CXX_FLAGS}")
        set(SYCL_TOOLCHAIN_LINKER_FLAGS "${SYCL_TOOLCHAIN_LINKER_FLAGS} ${SYCL_CXX_FLAGS}")
    else()
        message(FATAL_ERROR "Cannot compile a SYCL program with ${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_CXX_FLAGS}. Try a different compiler or disable SYCL.")
    endif()

    # Add kernel-splitting flag if available, both for compiling and linking
    set(SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS "-fsycl-device-code-split=per_kernel")
    gmx_check_source_compiles_with_flags(
        "${SAMPLE_SYCL_SOURCE}"
        "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS}"
        "CXX"
        SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS_RESULT
        )
    if (SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS_RESULT)
        set(SYCL_TOOLCHAIN_CXX_FLAGS "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS}")
        set(SYCL_TOOLCHAIN_LINKER_FLAGS "${SYCL_TOOLCHAIN_LINKER_FLAGS} ${SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS}")
    else()
        message(WARNING "Cannot compile SYCL with per-kernel device-code splitting. Simulations will work, but the first step will be much slower than it needs to be. Try a different compiler.")
    endif()

    # Add fast-math flag where available
    gmx_find_flag_for_source(
        SYCL_FAST_MATH_CXX_FLAGS_RESULT
        "${SAMPLE_SYCL_SOURCE}"
        "CXX"
        SYCL_TOOLCHAIN_CXX_FLAGS
        SYCL_FAST_MATH_CXX_FLAGS
        "-ffast-math" "/clang:-ffast-math")
    if (SYCL_FAST_MATH_CXX_FLAGS_RESULT)
        set(SYCL_TOOLCHAIN_CXX_FLAGS "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_FAST_MATH_CXX_FLAGS}")
    endif()

    if("${SYCL_CXX_FLAGS_EXTRA}" MATCHES "fsycl-targets=.*(nvptx64|amdgcn)")
        # When compiling for NVIDIA/AMD, Intel LLVM produces tons of harmless warnings, ignore them
        set(SYCL_WARNINGS_CXX_FLAGS "-Wno-linker-warnings -Wno-override-module -Wno-sycl-target")
        gmx_check_source_compiles_with_flags(
            "${SAMPLE_SYCL_SOURCE}"
            "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_WARNING_CXX_FLAGS}"
            "CXX"
            SYCL_WARNINGS_CXX_FLAGS_RESULT
            )
        if (SYCL_WARNINGS_CXX_FLAGS_RESULT)
            set(SYCL_TOOLCHAIN_CXX_FLAGS "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_WARNINGS_CXX_FLAGS}")
            set(SYCL_TOOLCHAIN_LINKER_FLAGS "${SYCL_TOOLCHAIN_LINKER_FLAGS} ${SYCL_WARNINGS_CXX_FLAGS}")
        endif()
    endif()

    if(GMX_GPU_FFT_VKFFT)
        include(gmxManageVkFft)
        if ("${SYCL_CXX_FLAGS_EXTRA}" MATCHES "fsycl-targets=.*nvptx64")
            gmx_manage_vkfft("CUDA")
        else()
            message(FATAL_ERROR "VkFFT can only be used with CUDA backend")
        endif()
        set(_sycl_has_valid_fft TRUE)
    endif()

    include(gmxManageFFTLibraries)

    if(GMX_GPU_FFT_MKL)
        #MKLROOT is set by gmxManageFFTLibraries.cmake
        find_library(mkl_sycl_PATH mkl_sycl PATHS "${MKLROOT}/lib/intel64" REQUIRED)
        mark_as_advanced(mkl_sycl_PATH)
        list(APPEND GMX_EXTRA_LIBRARIES "${mkl_sycl_PATH};OpenCL")

        set(CMAKE_REQUIRED_FLAGS "${SYCL_TOOLCHAIN_CXX_FLAGS}")
        set(CMAKE_REQUIRED_LIBRARIES "${GMX_EXTRA_LIBRARIES};${FFT_LIBRARIES}")
        check_cxx_source_compiles("
#include <oneapi/mkl/dfti.hpp>
int main() {
    oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE, oneapi::mkl::dft::domain::REAL> d({3,5,7});
    sycl::queue q;
    d.commit(q);
}"
          CAN_LINK_SYCL_MKL)
        unset(CMAKE_REQUIRED_FLAGS)
        unset(CMAKE_REQUIRED_LIBRARIES)
        unset(CMAKE_REQUIRED_INCLUDES)
        if (NOT CAN_LINK_SYCL_MKL)
            message(WARNING "Cannot link mkl_sycl. Make sure the MKL and compiler versions are compatible.")
        endif()

        set(_sycl_has_valid_fft TRUE)
    endif()

    if(GMX_GPU_FFT_BBFFT)
        # The double-batched FFT library is still called by its former
        # name bbfft in the implementation. For now, only the shared
        # libraries can link into GROMACS shared libraries.
        if (BUILD_SHARED_LIBS)
            find_package(bbfft-sycl 0.3.1 REQUIRED shared)
        else()
            find_package(bbfft-sycl 0.3.1 REQUIRED)
        endif()
        set(_sycl_has_valid_fft TRUE)
    endif()

    # convert the space-separated strings to lists
    separate_arguments(SYCL_TOOLCHAIN_CXX_FLAGS)
    list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS ${SYCL_CXX_FLAGS_EXTRA})
    separate_arguments(SYCL_TOOLCHAIN_LINKER_FLAGS)
    list(APPEND SYCL_TOOLCHAIN_LINKER_FLAGS ${SYCL_CXX_FLAGS_EXTRA})

    # Make strings for pretty-printing in gmx -version
    string(REPLACE ";" " " SYCL_TOOLCHAIN_CXX_FLAGS_STR "${SYCL_TOOLCHAIN_CXX_FLAGS}")
    string(STRIP "${SYCL_TOOLCHAIN_CXX_FLAGS_STR}" SYCL_TOOLCHAIN_CXX_FLAGS_STR)
    string(REPLACE ";" " " SYCL_TOOLCHAIN_LINKER_FLAGS_STR "${SYCL_TOOLCHAIN_LINKER_FLAGS}")
    string(STRIP "${SYCL_TOOLCHAIN_LINKER_FLAGS_STR}" SYCL_TOOLCHAIN_LINKER_FLAGS_STR)

    # Add function wrapper similar to the one used by ComputeCPP and hipSYCL
    function(add_sycl_to_target)
        cmake_parse_arguments(
            PARSE_ARGV 0 # No positional arguments
            ARGS # Prefix for the resulting variables
            "" # No options
            "TARGET" # One-value keyword
            "SOURCES" # Multi-value keyword
            )
        set_property(SOURCE ${ARGS_SOURCES} APPEND PROPERTY COMPILE_OPTIONS ${SYCL_TOOLCHAIN_CXX_FLAGS})
        target_link_options(${ARGS_TARGET} PRIVATE ${SYCL_TOOLCHAIN_LINKER_FLAGS})
    endfunction(add_sycl_to_target)
endif()

if (GMX_GPU_FFT_CUFFT AND GMX_USE_HEFFTE)
    set(_sycl_has_valid_fft TRUE)
    if (NOT DEFINED ENV{GITLAB_CI}) # Don't warn in CI builds
        message(WARNING "SYCL build with HeFFTe and cuFFT should only ever be used for testing")
    endif()
endif()

if(NOT ${_sycl_has_valid_fft} AND NOT GMX_GPU_FFT_LIBRARY STREQUAL "NONE")
    set(_hint "")
    if (GMX_GPU_FFT_CUFFT OR GMX_GPU_FFT_CLFFT)
        set(_hint " It is not supported with SYCL.")
    elseif (GMX_SYCL_HIPSYCL AND GMX_GPU_FFT_MKL)
        set(_hint " MKL is only supported with Intel DPC++ compiler, not with hipSYCL")
    endif()
    message(FATAL_ERROR "The selected GPU FFT library ${GMX_GPU_FFT_LIBRARY} is not compatible.${_hint}")
endif()

if(NOT ${_sycl_has_valid_fft} AND NOT DEFINED ENV{GITLAB_CI}) # Don't warn in CI builds
    message(WARNING "Building SYCL version without GPU FFT library.  Will not be able to perform FFTs on a GPU, which is not good for performance.")
endif()
