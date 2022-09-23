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

option(GMX_SYCL_HIPSYCL "Use hipSYCL instead of Intel/Clang for SYCL compilation" OFF)

if(GMX_DOUBLE)
    message(FATAL_ERROR "SYCL acceleration is not available in double precision")
endif()

include(gmxFindFlagsForSource)

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
    set(HIPSYCL_CLANG "${CMAKE_CXX_COMPILER}")
    # -Wno-unknown-cuda-version because Clang-11 complains about CUDA 11.0-11.2, despite working fine with them.
    # -Wno-unknown-attributes because hipSYCL does not support reqd_sub_group_size (because it can only do some sub group sizes).
    set(HIPSYCL_SYCLCC_EXTRA_ARGS "-Wno-unknown-cuda-version -Wno-unknown-attributes ${SYCL_CXX_FLAGS_EXTRA}")

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

    # Find a suitable rocFFT when hipSYCL is targeting AMD devices
    if (GMX_HIPSYCL_HAVE_HIP_TARGET)
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
    endif()
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
        # At least Intel dpcpp defaults to having SYCL enabled for all code. This leads to two problems:
        #
        # 1. Compiles take ~3x longer, since every file has to be compiled for multiple targets.
        # 2. We get a ton of warnings for the device-specific pass when the compiler sees our SIMD code.
        #
        # To avoid this, we attempt to find a flag to disable SYCL for non-SYCL files. Unfortunately,
        # when using gmx_find_flag_for_source() that includes calling check_cxx_compiler_flag(),
        # this in turn exposes a bug in dpcpp, where an object file compiles with -fno-sycl leads to
        # a failed link stage (when the same flag is not used). Since none of this is critical, we handle
        # it by merely checking if it works to compile a source fils with this flag, and choking if SYCL
        # is still enabled.
    
        if(NOT CHECK_DISABLE_SYCL_CXX_FLAGS_QUIETLY)
            message(STATUS "Checking for flags to disable SYCL")
        endif()
    
        gmx_check_source_compiles_with_flags(
            "int main() { return 0; }"
            "-fno-sycl"
            "CXX"
            DISABLE_SYCL_CXX_FLAGS_RESULT)
    
        if(DISABLE_SYCL_CXX_FLAGS_RESULT)
            set(SYCL_TOOLCHAIN_CXX_FLAGS "-fno-sycl")
        endif()
        if(NOT CHECK_DISABLE_SYCL_CXX_FLAGS_QUIETLY)
            if(DISABLE_SYCL_CXX_FLAGS_RESULT)
                message(STATUS "Checking for flags to disable SYCL - -fno-sycl")
            else()
                message(WARNING "Cannot find flags to disable SYCL for non-SYCL hardware-specific C++ code. Expect many warnings, but they are likely benign.")
            endif()
            set(CHECK_DISABLE_SYCL_CXX_FLAGS_QUIETLY 1 CACHE INTERNAL "Keep quiet on future calls to detect no-SYCL flags" FORCE)
        endif()
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
        set(SYCL_WARNINGS_CXX_FLAGS "-Wno-linker-warnings -Wno-override-module")
        gmx_check_source_compiles_with_flags(
            "${SAMPLE_SYCL_SOURCE}"
            "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_WARNING_CXX_FLAGS}"
            "CXX"
            SYCL_WARNINGS_CXX_FLAGS_RESULT
            )
        if (SYCL_WARNINGS_CXX_FLAGS_RESULT)
            set(SYCL_TOOLCHAIN_CXX_FLAGS "${SYCL_TOOLCHAIN_CXX_FLAGS} ${SYCL_WARNINGS_CXX_FLAGS}")
        endif()
    endif()

    include(gmxManageFFTLibraries)
    # Don't warn in CI builds
    if(NOT GMX_FFT_MKL AND NOT DEFINED ENV{GITLAB_CI})
        message(WARNING "Building SYCL version with ${GMX_FFT_LIBRARY} instead of MKL. GPU FFT is disabled!")
    endif()
    if(GMX_FFT_MKL)
	list(APPEND GMX_EXTRA_LIBRARIES "mkl_sycl;OpenCL")
    endif()

    # Add function wrapper similar to the one used by ComputeCPP and hipSYCL
    function(add_sycl_to_target)
        cmake_parse_arguments(
            PARSE_ARGV 0 # No positional arguments
            ARGS # Prefix for the resulting variables
            "" # No options
            "TARGET" # One-value keyword
            "SOURCES" # Multi-value keyword
            )
        # convert the space-separated string to a list
        separate_arguments(SYCL_TOOLCHAIN_CXX_FLAGS)
        set_property(SOURCE ${ARGS_SOURCES} APPEND PROPERTY COMPILE_OPTIONS
            ${SYCL_TOOLCHAIN_CXX_FLAGS}
            ${SYCL_CXX_FLAGS_EXTRA})
        if (WIN32) # Linking flags handling is not reliable on Windows, so we pass the bare minimum
            target_link_options(${ARGS_TARGET} PRIVATE ${SYCL_CXX_FLAGS})
        else()
            string(REPLACE " " ";" SYCL_TOOLCHAIN_LINKER_FLAGS_LIST "${SYCL_TOOLCHAIN_LINKER_FLAGS} ${SYCL_CXX_FLAGS_EXTRA}")
            target_link_options(${ARGS_TARGET} PRIVATE ${SYCL_TOOLCHAIN_LINKER_FLAGS_LIST})
        endif()
    endfunction(add_sycl_to_target)
endif()
