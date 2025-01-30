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

# CMake issue tracking the efforts to make a universal upstream module:
# https://gitlab.kitware.com/cmake/cmake/-/issues/21711

include(gmxFindFlagsForSource)

if(NOT GMX_GPU_SYCL OR GMX_SYCL_ACPP OR NOT GMX_SYCL_DPCPP)
    message(FATAL_ERROR "Internal error: OneAPI configuration script was included when it should not")
endif()

if(WIN32)
    if(NOT BUILD_SHARED_LIBS)
        message(FATAL_ERROR "SYCL with DPC++ on Windows doesn't work with static libraries. Set BUILD_SHARED_LIBS=on.")
        # Tested up to cmake 3.23.1 and icx 2022.1. Problem is order of exe link argument order. Works if gromacs.lib
        # and -fsycl both appear before -link. Not possible to change order from cmake script. cmake fix is WIP.
    endif()
endif()

# Prohibit the dpcpp compiler, but don't prohibit the use of e.g.
# /opt/dpcpp/oneapi/....../icpx etc.
if(CMAKE_CXX_COMPILER MATCHES "dpcpp$")
    message(FATAL_ERROR "Intel's \"dpcpp\" compiler is not supported; please use \"icpx\" for SYCL builds")
endif()

# Find the flags to enable SYCL with Intel compiler.
if(NOT CHECK_SYCL_CXX_FLAGS_QUIETLY)
    message(STATUS "Checking for flags to enable SYCL")
endif()
set(SAMPLE_SYCL_SOURCE
    "#include <sycl/sycl.hpp>
         int main(){
             sycl::queue q(sycl::default_selector_v);
             return 0;
         }")
set(SYCL_CXX_FLAGS "-fsycl")
gmx_check_source_compiles_with_flags(
    "${SAMPLE_SYCL_SOURCE}"
    "${SYCL_CXX_FLAGS}"
    "CXX"
    SYCL_CXX_FLAGS_RESULT
)
if (SYCL_CXX_FLAGS_RESULT)
    if(NOT CHECK_SYCL_CXX_FLAGS_QUIETLY)
        message(STATUS "Checking for flags to enable SYCL - ${SYCL_CXX_FLAGS}")
    endif()
    set(CHECK_SYCL_CXX_FLAGS_QUIETLY 1 CACHE INTERNAL "Keep quiet on future calls to detect SYCL flags" FORCE)
    set(SYCL_TOOLCHAIN_CXX_FLAGS ${SYCL_CXX_FLAGS})
    set(SYCL_TOOLCHAIN_LINKER_FLAGS ${SYCL_CXX_FLAGS})
else()
    message(FATAL_ERROR "Cannot compile a SYCL program with ${SYCL_CXX_FLAGS}. Try a different compiler or disable SYCL.")
endif()

# Try compiling an empty kernel to sniff out the enabled targets
set(SAMPLE_SYCL_KERNEL_PROBE_SOURCE
"#include <sycl/sycl.hpp>
int main()
{
    sycl::queue q;
    q.parallel_for<class K>(sycl::range<1>{ 16 }, [=](sycl::id<1> itemIdx) {
#if defined(__SYCL_DEVICE_ONLY__) && defined(__NVPTX__)
#    warning GMX_DPCPP_TEST_HAVE_CUDA_TARGET
#endif
#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
#    warning GMX_DPCPP_TEST_HAVE_HIP_TARGET
#    if __AMDGCN_WAVEFRONT_SIZE == 64
#        warning GMX_DPCPP_TEST_HAVE_HIP_WAVE64_TARGET
#    elif __AMDGCN_WAVEFRONT_SIZE == 32
#        warning GMX_DPCPP_TEST_HAVE_HIP_WAVE32_TARGET
#    endif
#endif
#if defined(__SYCL_DEVICE_ONLY__) && (defined(__SPIR__) || defined(__SPIRV__))
#    warning GMX_DPCPP_TEST_HAVE_INTEL_TARGET
#endif
    });
    return 0;
}"
)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
try_compile(
    SYCL_KERNEL_COMPILATION_WORKS
    SOURCES_TYPE NORMAL
    SOURCE_FROM_VAR "sample_sycl_kernel_probe.cpp" SAMPLE_SYCL_KERNEL_PROBE_SOURCE
    COMPILE_DEFINITIONS "${SYCL_TOOLCHAIN_CXX_FLAGS};${SYCL_CXX_FLAGS_EXTRA}"
    OUTPUT_VARIABLE SYCL_KERNEL_COMPILATION_LOG
)
unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
if (SYCL_KERNEL_COMPILATION_WORKS)
    foreach (target IN ITEMS CUDA HIP HIP_WAVE32 HIP_WAVE64 INTEL)
        if(SYCL_KERNEL_COMPILATION_LOG MATCHES "GMX_DPCPP_TEST_HAVE_${target}_TARGET")
            set(GMX_DPCPP_HAVE_${target}_TARGET ON CACHE INTERNAL "")
        else()
            set(GMX_DPCPP_HAVE_${target}_TARGET OFF CACHE INTERNAL "")
        endif()
        if (NOT GMX_DPCPP_TARGET_STATUS_REPORTED)
            message(STATUS "oneAPI DPC++ has ${target} target enabled: ${GMX_DPCPP_HAVE_${target}_TARGET}")
        endif()
    endforeach()
    set(GMX_DPCPP_TARGET_STATUS_REPORTED TRUE CACHE INTERNAL "")
else()
    message(FATAL_ERROR "Cannot compile sample SYCL kernel code:\n${SYCL_KERNEL_COMPILATION_LOG}")
endif()

if(NOT GMX_DPCPP_HAVE_CUDA_TARGET AND NOT GMX_DPCPP_HAVE_HIP_TARGET AND NOT GMX_DPCPP_HAVE_INTEL_TARGET)
    message(WARNING "SYCL oneAPI has no known GPU targets enabled!")
endif()
if(GMX_DPCPP_HAVE_HIP_TARGET AND NOT GMX_DPCPP_HAVE_HIP_WAVE32_TARGET AND NOT GMX_DPCPP_HAVE_HIP_WAVE64_TARGET)
    message(FATAL_ERROR "SYCL oneAPI detected the use of HIP target, but cannot determine wave size")
endif()

list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS ${SYCL_CXX_FLAGS_EXTRA})
list(APPEND SYCL_TOOLCHAIN_LINKER_FLAGS ${SYCL_CXX_FLAGS_EXTRA})

# Add kernel-splitting flag if available, both for compiling and linking
set(SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS "-fsycl-device-code-split=per_kernel")
list(JOIN SYCL_TOOLCHAIN_CXX_FLAGS " " SYCL_TOOLCHAIN_CXX_FLAGS_STR)
gmx_check_source_compiles_with_flags(
    "${SAMPLE_SYCL_SOURCE}"
    "${SYCL_TOOLCHAIN_CXX_FLAGS_STR} ${SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS}"
    "CXX"
    SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS_RESULT
)
if (SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS_RESULT)
    list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS ${SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS})
    list(APPEND SYCL_TOOLCHAIN_LINKER_FLAGS ${SYCL_DEVICE_CODE_SPLIT_CXX_FLAGS})
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
    list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS ${SYCL_FAST_MATH_CXX_FLAGS})
endif()

# We compile PME kernels for all possible sub-group sizes, so the warning is useless.
# Added after oneAPI 2024.0 (https://github.com/intel/llvm/pull/11991).
gmx_check_compiler_flag(
    "-Wno-incorrect-sub-group-size"
    "CXX"
    HAVE_W_NO_INCORRECT_SUB_GROUP_SIZE_RESULT
)
if (HAVE_W_NO_INCORRECT_SUB_GROUP_SIZE_RESULT)
    list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS "-Wno-incorrect-sub-group-size")
endif()

# Force small GRF size on PVC, see #5105
gmx_check_compiler_flag(
    "-ftarget-register-alloc-mode=pvc:small"
    "CXX"
    HAVE_TARGET_REGISTER_ALLOC_MODE_FLAG
)
if (HAVE_TARGET_REGISTER_ALLOC_MODE_FLAG)
    list(APPEND SYCL_TOOLCHAIN_LINKER_FLAGS "-ftarget-register-alloc-mode=pvc:small")
endif()

if(GMX_DPCPP_HAVE_CUDA_TARGET OR GMX_DPCPP_HAVE_HIP_TARGET)
    # When compiling for NVIDIA/AMD, Intel LLVM produces tons of harmless warnings, ignore them
    set(SYCL_WARNINGS_CXX_FLAGS "-Wno-linker-warnings;-Wno-override-module;-Wno-sycl-target")
    list(JOIN SYCL_TOOLCHAIN_CXX_FLAGS " " SYCL_TOOLCHAIN_CXX_FLAGS_STR)
    list(JOIN SYCL_WARNINGS_CXX_FLAGS " " SYCL_WARNINGS_CXX_FLAGS_STR)
    gmx_check_source_compiles_with_flags(
        "${SAMPLE_SYCL_SOURCE}"
        "${SYCL_TOOLCHAIN_CXX_FLAGS_STR} ${SYCL_WARNINGS_CXX_FLAGS_STR}"
        "CXX"
        SYCL_WARNINGS_CXX_FLAGS_RESULT
    )
    if (SYCL_WARNINGS_CXX_FLAGS_RESULT)
        list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS ${SYCL_WARNINGS_CXX_FLAGS})
        list(APPEND SYCL_TOOLCHAIN_LINKER_FLAGS ${SYCL_WARNINGS_CXX_FLAGS})
    endif()

    # Set GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT when targeting only devices with 64-wide execution
    set(_have_subgroup_not_64 OFF)
    set(_have_subgroup_64 OFF)
    if (GMX_DPCPP_HAVE_CUDA_TARGET OR GMX_DPCPP_HAVE_INTEL_TARGET OR GMX_DPCPP_HAVE_HIP_WAVE32_TARGET)
        set(_have_subgroup_not_64 ON) # We have AMD RDNA, NVIDIA, or Intel target(s)
    endif()
    # We assume that any GCN2-5 architecture (gfx7/8) and CDNA1-3 (gfx9 series) up until the time of writing of this conditional is 64-wide
    if (GMX_DPCPP_HAVE_HIP_WAVE64_TARGET)
        set(_have_subgroup_64 ON) # We have AMD GCN/CDNA target(s)
    endif()
    if (_have_subgroup_64 AND NOT _have_subgroup_not_64)
        option(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
            "Disable NBNXM GPU cluster pair splitting. Only supported with SYCL and 64-wide GPU architectures (like AMD GCN/CDNA)."
            ON)
        mark_as_advanced(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT)
    endif()
endif()

# The API is experimental and unstable, so we don't build it by default
option(GMX_SYCL_ENABLE_GRAPHS "Enable support for SYCL Graphs (experimental oneAPI feature)")
mark_as_advanced(GMX_SYCL_ENABLE_GRAPHS)

if(GMX_SYCL_ENABLE_GRAPHS)
    if(GMX_INTEL_LLVM AND GMX_INTEL_LLVM_VERSION LESS_EQUAL 20240002)
        message(FATAL_ERROR "SYCL Graph support in oneAPI 2024.0.2 is advertised but not sufficient")
    endif()
    list(JOIN SYCL_TOOLCHAIN_CXX_FLAGS " " CMAKE_REQUIRED_FLAGS)
    check_cxx_symbol_exists(SYCL_EXT_ONEAPI_GRAPH "sycl/sycl.hpp" HAVE_SYCL_EXT_ONEAPI_GRAPH)
    # Check if we have an aspect denoting partial support for graphs.
    # It is a temporary measure introduced in the open-source IntelLLVM around oneAPI 2024.3
    # to denote LevelZero devices supporting all graph features except for update,
    # with the intention to remove it in the future.
    check_cxx_source_compiles("#include <sycl/sycl.hpp>\n int main(){ auto a = sycl::aspect::ext_oneapi_limited_graph; }"
        HAVE_SYCL_ASPECT_EXT_ONEAPI_LIMITED_GRAPH)
    unset(CMAKE_REQUIRED_FLAGS)
    if(HAVE_SYCL_EXT_ONEAPI_GRAPH)
        set(GMX_HAVE_GPU_GRAPH_SUPPORT ON)
    else()
        message(FATAL_ERROR "SYCL Graph support requested, but SYCL_EXT_ONEAPI_GRAPH extension is not available")
    endif()
else()
    set(GMX_HAVE_GPU_GRAPH_SUPPORT OFF)
endif()

# The API is experimental and unstable, so we don't build it by default
option(GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API "Enable support for SYCL experimental submission API (experimental oneAPI feature)")
mark_as_advanced(GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API)
if(GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API)
    if(GMX_INTEL_LLVM AND GMX_INTEL_LLVM_VERSION LESS_EQUAL 20250001)
        message(FATAL_ERROR "SYCL submit support in oneAPI 2025.0 is broken, and earlier versions do not support the API; set GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API=OFF")
    endif()
endif()

if(GMX_GPU_FFT_VKFFT)
    include(gmxManageVkFft)
    set(_sycl_has_valid_fft TRUE)
endif()

include(gmxManageFFTLibraries)

if(GMX_GPU_FFT_MKL)
    #MKLROOT is set by gmxManageFFTLibraries.cmake
    find_library(mkl_sycl_PATH mkl_sycl PATHS "${MKLROOT}/lib" "${MKLROOT}/lib/intel64" REQUIRED)
    mark_as_advanced(mkl_sycl_PATH)
    list(APPEND GMX_EXTRA_LIBRARIES "${mkl_sycl_PATH};OpenCL")

    list(JOIN SYCL_TOOLCHAIN_CXX_FLAGS " " CMAKE_REQUIRED_FLAGS)
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

if(GMX_GPU_FFT_ONEMATH)
    include(gmxManageOneMath)
    gmx_manage_onemath()
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

# We disable warnings about functions deprecated in SYCL2020, because
# sometimes there is no widely-supported alternative.
# E.g., https://github.com/AdaptiveCpp/AdaptiveCpp/issues/1230.
# It would be good to occasionally remove this macro and clean
# up the code. But, as of November 2023, it produces too much noise
# without offering a targeted way to suppress specific warnings.
list(APPEND SYCL_TOOLCHAIN_CXX_FLAGS "-DSYCL2020_DISABLE_DEPRECATION_WARNINGS")

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
