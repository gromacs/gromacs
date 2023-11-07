#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2009- The GROMACS Authors
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

cmake_minimum_required(VERSION 3.18.4)

cmake_policy(SET CMP0022 NEW)
cmake_policy(SET CMP0048 NEW) # As of CMake 3.22, default is still "OLD"
cmake_policy(SET CMP0068 NEW) # From CMake 3.9
cmake_policy(SET CMP0074 NEW) # From CMake 3.12
if(POLICY CMP0131) # From CMake 3.24
    cmake_policy(SET CMP0131 NEW)
endif()
if(POLICY CMP0135) # From CMake 3.24
    # We don't care about DOWNLOAD_EXTRACT_TIMESTAMP; prefer the new behavior
    cmake_policy(SET CMP0135 NEW)
endif()
if(POLICY CMP0146) # From CMake 3.27
    # We still use FindCUDA
    cmake_policy(SET CMP0146 OLD)
endif()

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform
    )

if(APPLE)
    if(CMAKE_OSX_DEPLOYMENT_TARGET)
        # Providing a default value >=10.14 helps to find modern C++ compatibility,
        # such as by defaulting to the Clang libc++ instead of libstdc++.
        if(CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS "10.14")
            if((NOT DEFINED _gmx_osx_deployment_target_checked)
                    OR NOT (_gmx_osx_deployment_target_checked VERSION_EQUAL CMAKE_OSX_DEPLOYMENT_TARGET))
                message(WARNING "CMAKE_OSX_DEPLOYMENT_TARGET less than 10.14 may have compatibility problems.")
            endif()
        endif()
        set(_gmx_osx_deployment_target_checked ${CMAKE_OSX_DEPLOYMENT_TARGET} CACHE STRING
                "Last seen CMAKE_OSX_DEPLOYMENT_TARGET, if previously set." FORCE)
    endif()
    if(CMAKE_OSX_ARCHITECTURES)
        string(FIND ${CMAKE_OSX_ARCHITECTURES} i386 _substring_index)
        if(_substring_index GREATER_EQUAL 0)
            message(FATAL_ERROR "CMAKE_OSX_ARCHITECTURES includes i386, but GROMACS requires 64-bit architecture.")
        endif()
    endif()
endif()

# The GROMACS convention is that these are the version number of the next
# release that is going to be made from this branch.
project(Gromacs VERSION 2024.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(gmxManageCcache)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(LibStdCpp)

# Python is first referenced in gmxVersionInfo, so we perform the search early
# to find a suitable installation for all components.
include(gmxPythonDiscovery)
# Set up common version variables, as well as general information about
# the build tree (whether the build is from a source package or from a git
# repository).  Also declares a few functions that will be used for generating
# version info files later.
include(gmxBuildTreeInfo)
include(gmxVersionInfo)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND UNIX)
    set(CMAKE_INSTALL_PREFIX "/usr/local/gromacs" CACHE STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif()
if("${CMAKE_INSTALL_PREFIX}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "GROMACS cannot be installed into the build tree, choose a different location for CMAKE_INSTALL_PREFIX")
endif()

include(gmxBuildTypeReference)
include(gmxBuildTypeProfile)
include(gmxBuildTypeTSAN)
include(gmxBuildTypeASAN)
include(gmxBuildTypeMSAN)
include(gmxBuildTypeUBSAN)
include(gmxBuildTypeReleaseWithAssert)

set(valid_build_types "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "Reference" "RelWithAssert" "Profile" "TSAN" "ASAN" "MSAN" "UBSAN")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: ${valid_build_types}." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${valid_build_types})
endif()
if(NOT "${CMAKE_BUILD_TYPE}" IN_LIST valid_build_types)
    message(FATAL_ERROR "Unsupported value of -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}. "
        "Only ${valid_build_types} are allowed (case-sensitive)")
endif()
if(CMAKE_CONFIGURATION_TYPES)
    # Add appropriate GROMACS-specific build types for the Visual
    # Studio generator (Debug, Release, MinSizeRel and RelWithDebInfo
    # are already present by default).
    list(APPEND CMAKE_CONFIGURATION_TYPES "RelWithAssert" "Reference")
    list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
        "List of configuration types"
        FORCE)
endif()
set(build_types_with_explicit_flags RELEASE DEBUG RELWITHDEBINFO RELWITHASSERT MINSIZEREL PROFILE)

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

include(gmxCTestUtilities)
gmx_ctest_init()

include(gmxCPackUtilities)
gmx_cpack_init()

# Variables that accumulate stuff influencing the installed headers
set(INSTALLED_HEADER_INCLUDE_DIRS "")
set(INSTALLED_HEADER_DEFINITIONS "")

########################################################################
# Global non-cache variables for implementing the build system
########################################################################

# These variables collect libraries that GROMACS requires for
# linking. They should be appended to with list(APPEND ${name}
# new-library) calls. They are:
#  - Libraries that are required for libgromacs (only)
set(GMX_EXTRA_LIBRARIES "")
#  - Libraries that are required for all code in the repository
set(GMX_COMMON_LIBRARIES "")
#  - Libraries that all code linked against libgromacs needs
#    (i.e., something that is exposed in installed headers).
set(GMX_PUBLIC_LIBRARIES "")

########################################################################
# Check and warn if cache generated on a different host is being reused
########################################################################
if(CMAKE_HOST_UNIX)
    execute_process(COMMAND hostname
                    OUTPUT_VARIABLE TMP_HOSTNAME
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Only check for host name if not running in a CI environment, as the cache might
    # be reused there between different machines in different stages
    if(GMX_BUILD_HOSTNAME AND NOT "${GMX_BUILD_HOSTNAME}" STREQUAL "${TMP_HOSTNAME}"
            AND NOT DEFINED ENV{CI_JOB_ID})
        message(WARNING "
            The CMake cache, probably generated on a different host (${GMX_BUILD_HOSTNAME}),
            is being reused! This could lead to inconsistencies; therefore, it is
            recommended to regenerate the cache!")
    endif()
    set(GMX_BUILD_HOSTNAME "${TMP_HOSTNAME}" CACHE INTERNAL
            "Hostname of the machine where the cache was generated.")
endif()

########################################################################
# Detect architecture before setting options so we can alter defaults
########################################################################
# Detect the architecture the compiler is targetting, detect
# SIMD instructions possibilities on that hardware, suggest SIMD instruction set
# to use if none is specified, and populate the cache option for CPU
# SIMD.
include(gmxDetectTargetArchitecture)
gmx_detect_target_architecture()

########################################################################
# User input options                                                   #
########################################################################
include(gmxOptionUtilities)

set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")

option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" OFF)

option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREAD_MPI  "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" ON)

option(GMX_MIMIC "Enable MiMiC QM/MM interface (CPMD is required)" OFF)

option(GMX_CP2K "Enable CP2K QM/MM interface (CP2K 8.1 or later is required)" OFF)

# We need to enable Fortran, because CP2K will be linked
if(GMX_CP2K)
    enable_language(Fortran)
endif()

option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
mark_as_advanced(GMX_FAHCORE)

option(GMX_COOL_QUOTES "Enable GROMACS cool quotes" ON)
mark_as_advanced(GMX_COOL_QUOTES)
gmx_add_cache_dependency(GMX_COOL_QUOTES BOOL "NOT GMX_FAHCORE" OFF)

option(GMX_INSTALL_LEGACY_API "Install legacy headers" OFF)

gmx_option_multichoice(
    GMX_GPU
    "Framework for GPU acceleration"
    OFF
    OFF CUDA OpenCL SYCL)
if(GMX_GPU STREQUAL SYCL)
    gmx_option_multichoice(
        GMX_SYCL
        "SYCL implementation to use for GPU acceleration"
        AUTO
        DPCPP ACPP AUTO)
    if (GMX_SYCL STREQUAL "AUTO")
        if(GMX_SYCL_HIPSYCL) # Check legacy flag
            message(WARNING "GMX_SYCL_HIPSYCL option is deprecated; please use -DGMX_SYCL=ACPP")
            set(GMX_SYCL "ACPP")
        else()
            # Check direct compiler support (DPC++)
            include(CheckCXXCompilerFlag)
            check_cxx_compiler_flag("-fsycl" HAVE_SYCL_DPCPP)
            if(HAVE_SYCL_DPCPP)
                set(GMX_SYCL "DPCPP")
                message(STATUS "GROMACS SYCL backend: Setting GMX_SYCL=DPCPP (Intel oneAPI DPC++)")
            else()
                set(GMX_SYCL "ACPP")
                message(STATUS "GROMACS SYCL backend: Setting GMX_SYCL=ACPP (AdaptiveCpp/hipSYCL)")
            endif()
        endif ()
    endif ()
endif()

gmx_option_multichoice(
    GMX_SIMD
    "SIMD instruction set for CPU kernels and compiler optimization"
    "AUTO"
    AUTO None SSE2 SSE4.1 AVX_128_FMA AVX_256 AVX2_256 AVX2_128 AVX_512 AVX_512_KNL ARM_NEON_ASIMD ARM_SVE IBM_VSX Reference)

include(gmxTestIntelLLVM)

if (GMX_INTEL_LLVM)
    set(GMX_FFT_LIBRARY_DEFAULT "mkl")
else()
    set(GMX_FFT_LIBRARY_DEFAULT "fftw3")
endif()

gmx_option_multichoice(
    GMX_FFT_LIBRARY
    "FFT library"
    "${GMX_FFT_LIBRARY_DEFAULT}"
    fftw3 mkl "fftpack[built-in]")
gmx_dependent_option(
    GMX_BUILD_OWN_FFTW
    "Download and build FFTW 3 during the GROMACS build process, rather than fall back on the really slow fftpack."
    OFF
    "GMX_FFT_LIBRARY STREQUAL FFTW3")
gmx_dependent_option(
    GMX_DISABLE_FFTW_MEASURE
    "Do not optimize FFTW setups (not needed with SSE)"
    OFF
    "GMX_FFT_LIBRARY STREQUAL FFTW3")
mark_as_advanced(GMX_BUILD_OWN_FFTW)
mark_as_advanced(GMX_DISABLE_FFTW_MEASURE)

gmx_dependent_option(
    GMX_USE_HEFFTE
    "Use HeFFTe for distributed FFT support. Used with CUDA backend"
    OFF
    "GMX_GPU STREQUAL CUDA OR GMX_GPU STREQUAL SYCL;GMX_MPI")
gmx_dependent_option(
    GMX_USE_CUFFTMP
    "Use cuFFTMp for distributed FFT support. Used with CUDA backend"
    OFF
    "GMX_GPU STREQUAL CUDA;GMX_MPI")

# Here the default GPU FFT library is set up depending
# on the build configuration. Other choices are
# available, these are just the defaults.
#
# =================+======================================
# Configuration    |  Default GPU FFT library
# =================+======================================
# CUDA             | cuFFT
# OpenCL           | VkFFT or clFFT according to toolchain
# SYCL via hipSYCL | VkFFT
# SYCL via DPCPP   | MKL
# =================+======================================
if (GMX_GPU)
    if (GMX_GPU STREQUAL CUDA)
        set(GMX_GPU_FFT_LIBRARY_DEFAULT "cuFFT")
    elseif(GMX_GPU STREQUAL OPENCL)
        if (APPLE OR MSVC)
            set(GMX_GPU_FFT_LIBRARY_DEFAULT "VkFFT")
        else()
            set(GMX_GPU_FFT_LIBRARY_DEFAULT "clFFT")
        endif()
    elseif(GMX_GPU STREQUAL SYCL)
        if(GMX_SYCL STREQUAL ACPP)
            set(GMX_GPU_FFT_LIBRARY_DEFAULT "VkFFT")
        else()
            set(GMX_GPU_FFT_LIBRARY_DEFAULT "MKL")
        endif()
    endif()

    gmx_option_multichoice(
        GMX_GPU_FFT_LIBRARY
        "GPU FFT library"
        "${GMX_GPU_FFT_LIBRARY_DEFAULT}"
        cuFFT clFFT VkFFT MKL rocFFT BBFFT none)
    # The validity of chosen option is later verified by the selected backend

    foreach (_gpu_fft_library cuFFT clFFT VkFFT MKL rocFFT BBFFT)
        string(TOUPPER "${_gpu_fft_library}" _gpu_fft_library_upper)
        if (GMX_GPU_FFT_LIBRARY STREQUAL ${_gpu_fft_library_upper})
            if (NOT GMX_GPU_FFT_QUIET_AFTER_FIRST_RUN)
                message(STATUS "Selected GPU FFT library - ${_gpu_fft_library}")
            endif()
            set(GMX_GPU_FFT_QUIET_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future runs of cmake")
            set(_value TRUE)
        else()
            set(_value FALSE)
        endif()
        set("GMX_GPU_FFT_${_gpu_fft_library_upper}" "${_value}" CACHE INTERNAL "Use ${_gpu_fft_library} library for FFTs on GPUs")
    endforeach()
endif()
gmx_dependent_option(GMX_NVSHMEM "Build a parallel NVSHMEM multi-GPU code of GROMACS" OFF "GMX_GPU STREQUAL CUDA;GMX_MPI")

if (GMX_NVSHMEM)
    if (NVSHMEM_ROOT)
        set(GMX_NVSHMEM_HOME ${NVSHMEM_ROOT} CACHE PATH "path to NVSHMEM")
    elseif (DEFINED ENV{NVHPC_ROOT})
        set(GMX_NVSHMEM_HOME "$ENV{NVHPC_ROOT}/comm_libs/nvshmem" CACHE PATH "path to NVSHMEM")
    else()
        message(FATAL_ERROR "NVSHMEM_ROOT path must be defined if GMX_NVSHMEM is set")
    endif()
endif()

gmx_dependent_cache_variable(GMX_SIMD_REF_FLOAT_WIDTH  "Reference SIMD single precision width" STRING "4" "GMX_SIMD STREQUAL REFERENCE")
gmx_dependent_cache_variable(GMX_SIMD_REF_DOUBLE_WIDTH "Reference SIMD double precision width" STRING "4" "GMX_SIMD STREQUAL REFERENCE")

# This should be moved to a separate NBNXN cmake module when that code is cleaned up and modularized

option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
mark_as_advanced(GMX_BROKEN_CALLOC)

option(GMX_OPENMP "Enable OpenMP-based multithreading" ON)

option(GMX_USE_TNG "Use the TNG library for trajectory I/O" ON)

option(GMX_CYCLE_SUBCOUNTERS "Enable cycle subcounters to get a more detailed cycle timings" OFF)
mark_as_advanced(GMX_CYCLE_SUBCOUNTERS)

option(GMX_SKIP_DEFAULT_CFLAGS "Don't automatically add suggested/required Compiler flags." OFF)
mark_as_advanced(GMX_SKIP_DEFAULT_CFLAGS)

option(GMX_BUILD_FOR_COVERAGE
       "Tune build for better code coverage metrics (e.g., disable asserts)"
       OFF)
mark_as_advanced(GMX_BUILD_FOR_COVERAGE)

option(GMX_DEVELOPER_BUILD
    "Enable Developer convenience features: always build unit-tests"
    OFF)
mark_as_advanced(GMX_DEVELOPER_BUILD)

gmx_set_boolean(GMX_COMPILER_WARNINGS_DEFAULT "NOT SOURCE_IS_SOURCE_DISTRIBUTION")
option(GMX_COMPILER_WARNINGS
    "Enable a default set of compiler warnings"
    ${GMX_COMPILER_WARNINGS_DEFAULT})
mark_as_advanced(GMX_COMPILER_WARNINGS)
# Always turn on compiler warnings with a developer build.
gmx_add_cache_dependency(GMX_COMPILER_WARNINGS BOOL "NOT GMX_DEVELOPER_BUILD" ON)

option(GMX_BUILD_SHARED_EXE
    "Build executables as shared binaries. If not set, this disables rpath and dynamic linker flags in an attempt to build a static binary, but this may require setting up the toolchain properly and making appropriate libraries available."
    ON)
mark_as_advanced(GMX_BUILD_SHARED_EXE)

option(GMX_PHYSICAL_VALIDATION
       "Include physical validation tests in ctest environment. These can then be called using 'make check-phys' or
       'make check-all'. Warning: Running the physical validation tests takes significantly more time than other tests!"
       OFF)
mark_as_advanced(GMX_PHYSICAL_VALIDATION)

######################################################################
# Detect OpenMP support
######################################################################
# The OpenMP detection _must_ come before tests for other CFLAGS.
include(gmxManageOpenMP)



######################################################################
# Compiler tests
# These need to be done early (before further tests).
#####################################################################

include(gmxCFlags)
gmx_c_flags()

# These variables should be used for CMake-style lists (ie. separated
# by semicolons) of additional compiler flags which are not generated
# in gmxCFlags nor are SIMD or MPI related.
#
# TODO These variables should be consolidated into
# EXTRA_COMPILER_FLAGS so that we we don't perpetrate bugs where
# things that work in C compilation (e.g. merging from old branches)
# might not also work for C++ compilation.
set(EXTRA_C_FLAGS "")
set(EXTRA_CXX_FLAGS "")

# Run through a number of tests for buggy compilers and other issues
include(gmxTestCompilerProblems)
gmx_test_compiler_problems()

# Implement double-precision option. This is complicated because we
# need installed headers to use the precision mode of the build that
# produced the library, but cannot use config.h in that case. We also
# want such variables to always have a definition, because #if is more
# robust than #ifdef. So, we put this value on the compiler command
# line in all cases.
if(GMX_DOUBLE)
    set(GMX_DOUBLE_VALUE 1)
else()
    set(GMX_DOUBLE_VALUE 0)
endif()
add_definitions(-DGMX_DOUBLE=${GMX_DOUBLE_VALUE})
list(APPEND INSTALLED_HEADER_DEFINITIONS "-DGMX_DOUBLE=${GMX_DOUBLE_VALUE}")

option(GMX_IMD "Enable Interactive Molecular Dynamics (IMD) sessions, e.g. with VMD" ON)
if(GMX_IMD AND WIN32)
    list(APPEND GMX_EXTRA_LIBRARIES "wsock32")
endif()


########################################################################
# Parse present/past contributors and project leaders from AUTHORS     #
########################################################################
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS AUTHORS)
file(READ "${CMAKE_SOURCE_DIR}/AUTHORS" AUTHORS)
string(REGEX REPLACE " <[^<]+>" "" AUTHORS_NO_EMAIL "${AUTHORS}")
string(REGEX MATCH "CURRENT CONTRIBUTORS\n=+\n(.*)\n\n+PREVIOUS CONTRIBUTORS\n=+\n(.*)\n\n+CURRENT PROJECT LEADERS\n=+\n(.*)\n\n+" MATCH "${AUTHORS_NO_EMAIL}")

set(GMX_CURRENT_CONTRIBUTORS ${CMAKE_MATCH_1})
set(GMX_PREVIOUS_CONTRIBUTORS ${CMAKE_MATCH_2})
set(GMX_CURRENT_PROJECT_LEADERS ${CMAKE_MATCH_3})

string(REGEX REPLACE "\n" ", " GMX_CURRENT_CONTRIBUTORS_STRING "${GMX_CURRENT_CONTRIBUTORS}")
string(REGEX REPLACE "\n" ", " GMX_PREVIOUS_CONTRIBUTORS_STRING "${GMX_PREVIOUS_CONTRIBUTORS}")
string(REGEX REPLACE "\n" ", " GMX_CURRENT_PROJECT_LEADERS_STRING "${GMX_CURRENT_PROJECT_LEADERS}")

string(REGEX REPLACE "([^\n]+)" "    \"\\1\"," GMX_CURRENT_CONTRIBUTORS "${GMX_CURRENT_CONTRIBUTORS}")
string(REGEX REPLACE "([^\n]+)" "    \"\\1\"," GMX_PREVIOUS_CONTRIBUTORS "${GMX_PREVIOUS_CONTRIBUTORS}")
string(REGEX REPLACE "([^\n]+)" "    \"\\1\"," GMX_CURRENT_PROJECT_LEADERS "${GMX_CURRENT_PROJECT_LEADERS}")

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
check_include_files(unistd.h     HAVE_UNISTD_H)
check_include_files(pwd.h        HAVE_PWD_H)
check_include_files(dirent.h     HAVE_DIRENT_H)
check_include_files(time.h       HAVE_TIME_H)
check_include_files(sys/time.h   HAVE_SYS_TIME_H)
check_include_files(io.h         HAVE_IO_H)
check_include_files(sched.h      HAVE_SCHED_H)
check_include_files(xmmintrin.h  HAVE_XMMINTRIN_H)

include(CheckCXXSymbolExists)
check_cxx_symbol_exists(gettimeofday      sys/time.h   HAVE_GETTIMEOFDAY)
check_cxx_symbol_exists(sysconf           unistd.h     HAVE_SYSCONF)
check_cxx_symbol_exists(nice              unistd.h     HAVE_NICE)
check_cxx_symbol_exists(fsync             unistd.h     HAVE_FSYNC)
check_cxx_symbol_exists(_fileno           stdio.h      HAVE__FILENO)
check_cxx_symbol_exists(fileno            stdio.h      HAVE_FILENO)
check_cxx_symbol_exists(_commit           io.h         HAVE__COMMIT)
check_cxx_symbol_exists(sigaction         signal.h     HAVE_SIGACTION)

# We cannot check for the __builtins as symbols, but check if code compiles
check_cxx_source_compiles("int main(){ return __builtin_clz(1);}"   HAVE_BUILTIN_CLZ)
check_cxx_source_compiles("int main(){ return __builtin_clzll(1);}" HAVE_BUILTIN_CLZLL)
if(MSVC)
    check_cxx_source_compiles("#include <intrin.h>\n int main(){unsigned long r;unsigned long i=1;_BitScanReverse(&r,i);return r;}" HAVE_BITSCANREVERSE)
    check_cxx_source_compiles("#include <intrin.h>\n int main(){unsigned long r;unsigned __int64 i=1;_BitScanReverse64(&r,i);return r;}" HAVE_BITSCANREVERSE64)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")
    check_cxx_source_compiles("int main(){ return __cntlz4(1);}" HAVE_CNTLZ4)
    check_cxx_source_compiles("int main(){ return __cntlz8(1);}" HAVE_CNTLZ8)
endif()

include(CheckLibraryExists)
find_library(HAVE_LIBM m)
mark_as_advanced(HAVE_LIBM)
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)
check_library_exists(m feenableexcept "" HAVE_FEENABLEEXCEPT)
check_library_exists(m fedisableexcept "" HAVE_FEDISABLEEXCEPT)

include(TestSchedAffinity)
test_sched_affinity(HAVE_SCHED_AFFINITY)

# Aligned memory allocation. We need to check for both mm_malloc(),
# posix_memalign(), memalign(), and on windows also _aligned_malloc()
include(gmxTestMMMalloc)
gmx_test_mm_malloc(HAVE__MM_MALLOC)
check_cxx_symbol_exists(posix_memalign    stdlib.h     HAVE_POSIX_MEMALIGN)
check_cxx_symbol_exists(memalign          stdlib.h     HAVE_MEMALIGN)
if(MSVC)
    # No need to waste time on this test on platforms where it will never be true
    check_cxx_symbol_exists(_aligned_malloc   stdlib.h     HAVE__ALIGNED_MALLOC)
endif()

include(TestBigEndian)
test_big_endian(GMX_INTEGER_BIG_ENDIAN)

gmx_set_boolean(GMX_USE_NICE "HAVE_UNISTD_H AND HAVE_NICE")

########################################################################
#Process MPI settings
########################################################################
include(gmxManageMPI)

########################################################################
#Process MiMiC settings
########################################################################
include(gmxManageMimic)

########################################################################
#Process CP2K settings
########################################################################
include(gmxManageCP2K)

########################################################################
#Process shared/static library settings
########################################################################
include(gmxManageSharedLibraries)


########################################################################
# Specify install locations
########################################################################
# Use GNUInstallDirs to set paths on multiarch systems.
include(GNUInstallDirs)

set(GMX_INSTALL_DATASUBDIR "gromacs" CACHE STRING "Subdirectory for GROMACS data under CMAKE_INSTALL_DATADIR")
mark_as_advanced(GMX_INSTALL_DATASUBDIR)

# Internal convenience so we do not have to join two path segments in the code
set(GMX_INSTALL_GMXDATADIR ${CMAKE_INSTALL_DATADIR}/${GMX_INSTALL_DATASUBDIR})

# If the nesting level wrt. the installation root is changed,
# gromacs-config.cmake.cmakein needs to be adapted.
set(GMX_INSTALL_CMAKEDIR  ${CMAKE_INSTALL_DATAROOTDIR}/cmake)

# TODO: Make GMXRC adapt if this is changed
set(GMX_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

list(APPEND INSTALLED_HEADER_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR})

# Binary and library suffix options
include(gmxManageSuffixes)


########################################################################
# Find external packages                                               #
########################################################################

option(GMX_HWLOC "Use hwloc portable hardware locality library" OFF)

if (GMX_HWLOC)
    # Find quietly the second time.
    if (HWLOC_FIND_QUIETLY_AFTER_FIRST_RUN)
        set(HWLOC_FIND_QUIETLY TRUE)
    endif()
    find_package(HWLOC 1.5) 

    if (HWLOC_FOUND)
        if (HWLOC_LIBRARIES MATCHES ".a$")
            set(_STATIC_HWLOC TRUE)
        endif()

        gmx_check_if_changed(HWLOC_FOUND_CHANGED HWLOC_FOUND)
        if (_STATIC_HWLOC AND HWLOC_FOUND_CHANGED AND NOT GMX_HWLOC_FORCE)
            message(STATUS "Static hwloc library found, will not attempt using it as it could lead to link-time errors. To use the detected library, manually set GMX_HWLOC=ON and you will likely have to pass appropriate linker flags too to satisfy the link-time dependencies of your hwloc library. Try \"pkg-config --libs --static hwloc\" for suggestions on what you will need.")
            set(GMX_USE_HWLOC OFF)
        else()
            set(GMX_USE_HWLOC ON)
        endif()

        if (GMX_USE_HWLOC)
            include_directories(SYSTEM ${HWLOC_INCLUDE_DIRS})
            list(APPEND GMX_EXTRA_LIBRARIES ${HWLOC_LIBRARIES})
        endif()
    elseif(GMX_HWLOC_FORCE)
        message(FATAL_ERROR "HWLOC package support required, but not found.")
    endif()

    if (HWLOC_FOUND AND HWLOC_VERSION VERSION_LESS "2")
        message(STATUS "Support for hwloc versions 1.x is deprecated")
    endif()

    set(HWLOC_FIND_QUIETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to find HWLOC")
endif()

option(GMX_EXTERNAL_TINYXML2 "Use external TinyXML-2 instead of compiling the version bundled with GROMACS." OFF)
mark_as_advanced(GMX_EXTERNAL_TINYXML2)
if(GMX_EXTERNAL_TINYXML2)
    # Find an external TinyXML-2 library.
    find_package(TinyXML2 6.0.0)
    set(HAVE_TINYXML2 ${TinyXML2_FOUND})
    if(NOT HAVE_TINYXML2)
        message(FATAL_ERROR "External TinyXML-2 could not be found, please adjust your search paths")
    endif()
endif()

option(GMX_EXTRAE "Add support for tracing using EXTRAE" OFF)
mark_as_advanced(GMX_EXTRAE)

if (GMX_EXTRAE)
  find_package(EXTRAE)
  if(EXTRAE_FOUND)
    include_directories(SYSTEM ${EXTRAE_INCLUDE_DIR})
    set(HAVE_EXTRAE 1)
  else()
    message(FATAL_ERROR "EXTRAE library was not found. Please add the correct path to CMAKE_PREFIX_PATH")
  endif()
endif()

include(ThreadMPI)
# Enable core threading facilities
tmpi_enable_core("${CMAKE_SOURCE_DIR}/src/external/thread_mpi/include")
if(GMX_THREAD_MPI)
    # enable MPI functions
    tmpi_enable()
endif()
# If atomics are manually disabled a define is needed because atomics.h doesn't depend on config.h
if (TMPI_ATOMICS_DISABLED)
   add_definitions(-DTMPI_ATOMICS_DISABLED)
endif()

include(gmxManageTNG)

include(gmxManageLmfit)

include(gmxManageMuparser)

include(gmxManageColvars)

##################################################
# Process SIMD instruction settings
##################################################
# This checks what flags to add in order to
# support the SIMD instructions we need, it sets
# correct defines for the SIMD instructions supported,
# and adds advanced options to control accuracy
# for SIMD math operations.
include(gmxManageSimd)
gmx_manage_simd()

# The earliest version of the CUDA toolkit that supports c++17 is 11.0
set(REQUIRED_CUDA_VERSION 11.0)
set(REQUIRED_CUDA_COMPUTE_CAPABILITY 3.5)

if(GMX_GPU)

    string(TOUPPER "${GMX_GPU}" _gmx_gpu_uppercase)
    if(${_gmx_gpu_uppercase} STREQUAL "CUDA")
        include(gmxManageCuda)
    elseif(${_gmx_gpu_uppercase} STREQUAL "OPENCL")
        message(STATUS "GPU support with OpenCL is deprecated. It is still fully supported (and " 
            "recommended for AMD, Intel, and Apple GPUs). It may be replaced by different "
            "approaches in future releases of GROMACS.")
        include(gmxManageOpenCL)
    elseif(${_gmx_gpu_uppercase} STREQUAL "SYCL")
        include(gmxManageSycl)
    endif()
    if(NOT GMX_OPENMP)
        message(WARNING "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading, which is currently not enabled.")
    endif()

    if (GMX_OPENCL_NB_CLUSTER_SIZE)
        message(WARNING "GMX_OPENCL_NB_CLUSTER_SIZE is deprecated, use GMX_GPU_NB_CLUSTER_SIZE instead")
    endif()
    if (GMX_OPENCL_NB_CLUSTER_SIZE AND GMX_GPU_NB_CLUSTER_SIZE)
        if (NOT ${GMX_OPENCL_NB_CLUSTER_SIZE} EQUAL ${GMX_GPU_NB_CLUSTER_SIZE})
            message(FATAL_ERROR "Mismatching values passed to GMX_OPENCL_NB_CLUSTER_SIZE and GMX_GPU_NB_CLUSTER_SIZE; the former is deprecated, use only the latter!")
        endif()
    endif()
    # Only OpenCL and SYCL support changing the default cluster size
    if (${_gmx_gpu_uppercase} STREQUAL "CUDA")
        if (GMX_GPU_NB_CLUSTER_SIZE AND NOT "${GMX_GPU_NB_CLUSTER_SIZE}" EQUAL 8)
            message(FATAL_ERROR "Changing GMX_GPU_NB_CLUSTER_SIZE is not supported in CUDA (the default GMX_GPU_NB_CLUSTER_SIZE=8 is used)")
        endif()
        set(GMX_GPU_NB_CLUSTER_SIZE 8 CACHE STRING "Cluster size used by the nonbonded kernel.")
    else()
        # use the legacy GMX_OPENCL_NB_CLUSTER_SIZE variable if set, otherwise set the defaults
        if (GMX_OPENCL_NB_CLUSTER_SIZE)
            set(_gmx_gpu_nb_cluster_size_value ${GMX_OPENCL_NB_CLUSTER_SIZE})
        else()
            # default cluster size is 8 with OpenCL and 4 with SYCL for now
            if(${_gmx_gpu_uppercase} STREQUAL "OPENCL")
                set(_gmx_gpu_nb_cluster_size_value 8)
            elseif(GMX_GPU_SYCL)
                if (GMX_SYCL_ACPP AND NOT GMX_HIPSYCL_HAVE_LEVELZERO_TARGET)
                    set(_gmx_gpu_nb_cluster_size_value 8)
                else()
                    # Either DPCPP or hipSYCL targeting Intel Level0
                    set(_gmx_gpu_nb_cluster_size_value 4)
                endif()
            endif()
        endif()
        set(GMX_GPU_NB_CLUSTER_SIZE ${_gmx_gpu_nb_cluster_size_value} CACHE STRING "Cluster size used by the nonbonded kernel. Set to 4 for Intel GPUs.")
        mark_as_advanced(GMX_GPU_NB_CLUSTER_SIZE)
    endif()

    if (NOT (${_gmx_gpu_uppercase} STREQUAL "SYCL"))
        if (GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT)
            message(FATAL_ERROR "Disabling cluster pair splitting is only supported in SYCL, set GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT=off otherwise")
        endif()
    endif()
endif()
set(GMX_GPU_NB_NUM_CLUSTER_PER_CELL_X 2 CACHE STRING "Number of clusters along X in a pair-search grid cell for GPU lists")
set(GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Y 2 CACHE STRING "Number of clusters along Y in a pair-search grid cell for GPU lists")
set(GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Z 2 CACHE STRING "Number of clusters along Z in a pair-search grid cell for GPU lists")
mark_as_advanced(GMX_GPU_NB_NUM_CLUSTER_PER_CELL_X)
mark_as_advanced(GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Y)
mark_as_advanced(GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Z)

# For build with CUDA/SYCL and Lib-MPI, check if underlying MPI implementation is GPU-aware
# GPU-aware MPI allows direct GPU communication without staging data through host
if((GMX_GPU_CUDA OR GMX_GPU_SYCL) AND GMX_LIB_MPI)
    include(gmxManageGpuAwareMpi)
else()
    set(HAVE_MPI_EXT 0)
    set(MPI_SUPPORTS_CUDA_AWARE_DETECTION 0)
    set(MPI_SUPPORTS_HIP_AWARE_DETECTION 0)
    set(MPI_SUPPORTS_ZE_AWARE_DETECTION 0)
endif()

if(CYGWIN)
    set(GMX_CYGWIN 1)
endif()

if(GMX_GPU_CUDA AND NOT GMX_GPU_FFT_CUFFT)
    message(FATAL_ERROR "The CUDA build only supports cuFFT GPU FFT library")
endif()

if(GMX_USE_HEFFTE)
    if(NOT GMX_GPU_CUDA AND NOT GMX_GPU_SYCL)
        message(FATAL_ERROR "HeFFTe support requires a CUDA or SYCL build")
    endif()
    if(NOT GMX_LIB_MPI)
        message(FATAL_ERROR "HeFFTe support requires a library MPI build")
    endif()
    if(GMX_GPU_CUDA)
        find_package(Heffte 2.2.0 REQUIRED CUDA)
    elseif(GMX_GPU_SYCL)
        if(GMX_GPU_FFT_ROCFFT)
            find_package(Heffte 2.2.0 REQUIRED ROCM)
        elseif(GMX_GPU_FFT_CUFFT)
            find_package(Heffte 2.2.0 REQUIRED CUDA)
            if (NOT DEFINED ENV{GITLAB_CI}) # Don't warn in CI builds
                message(WARNING "HeFFTe with cuFFT backend should be used in a CUDA build")
            endif()
        elseif(GMX_GPU_FFT_MKL)
            find_package(Heffte 2.2.0 REQUIRED ONEAPI)
        else()
            message(FATAL_ERROR "Your GPU FFT library is incompatible for use in a GROMACS SYCL build with heFFTe. "
                    "Use -DGMX_GPU_FFT_LIBRARY=rocFFT or -DGMX_GPU_FFT_LIBRARY=MKL or -DGMX_GPU_FFT_LIBRARY=CUFFT")
        endif()
    endif()
    set(GMX_USE_Heffte ${Heffte_FOUND})
    # Find quietly in subsequent passes
    set(Heffte_FIND_QUIETLY TRUE CACHE INTERNAL "Find HeFFTe quietly in future CMake passes")
endif()

if(GMX_USE_CUFFTMP)
    if(NOT GMX_GPU_CUDA)
        message(FATAL_ERROR "cuFFTMp support requires a CUDA build")
    endif()
    if(NOT GMX_LIB_MPI)
        message(FATAL_ERROR "cuFFTMp support requires a library MPI build")
    endif()
    find_package(cuFFTMp REQUIRED CUDA)
    set(GMX_USE_cuFFTMp ${cuFFTMp_FOUND})
endif()

# Set NVSHMEM host and device libs
include(gmxManageNvshmem)

if(WIN32)
    set(GMX_NATIVE_WINDOWS 1)
    # This makes windows.h not declare min/max as macros that would break
    # C++ code using std::min/std::max.
    add_definitions(-DNOMINMAX)
endif()

option(GMX_BUILD_UNITTESTS "Build unit tests with BUILD_TESTING" ON)
mark_as_advanced(GMX_BUILD_UNITTESTS)
gmx_add_cache_dependency(GMX_BUILD_UNITTESTS BOOL BUILD_TESTING OFF)

########################################################################
# Our own GROMACS tests
########################################################################

include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/src/external)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

include(gmxSetBuildInformation)
gmx_set_build_information()

# Anything but truly ancient x86 hardware should support rdtscp, so we enable it by default.
# The inline assembly calling it is only ever compiled on x86, so defaulting to ON is OK.
option(GMX_USE_RDTSCP "Use low-latency RDTSCP instruction for x86 CPU-based timers for mdrun execution; might need to be off when compiling for heterogeneous environments" ON)
mark_as_advanced(GMX_USE_RDTSCP)

include(gmxTestLargeFiles)
gmx_test_large_files(GMX_LARGEFILES)

include(gmxTestSignal)
gmx_test_sigusr1(HAVE_SIGUSR1)

include(gmxTestPipes)
gmx_test_pipes(HAVE_PIPES)

include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
# Darwin has system XDR, but it uses a three-argument flavour of
# xdrproc_t that it guarantees will still work if you pass the normal
# two-argument xdr filters, but gcc 8 warns about the cast necessary
# to do that, so it's simpler to just use our own XDR library.
#
# TODO It would be better to craft a cmake test which fails if such
# XDR operations cause warnings, and succeeds otherwise, because it is
# generally preferable to use system libraries where possible.
if(NOT GMX_SYSTEM_XDR OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(GMX_INTERNAL_XDR 1)
endif()

##################################################
# Process FFT library settings
##################################################
include(gmxManageFFTLibraries)


include(gmxManageLinearAlgebraLibraries)

include(gmxManagePluginSupport)
gmx_manage_plugin_support()

if(GMX_USE_PLUGINS)
    if(NOT GMX_VMD_PLUGIN_PATH)
        find_package(VMD)
    endif()
endif()

# Link real-time library for POSIX timers. The check for clock_gettime
# confirms the linkability of rt.
if(HAVE_TIME_H AND HAVE_UNISTD_H AND HAVE_CLOCK_GETTIME)
    list(APPEND GMX_EXTRA_LIBRARIES rt)
endif()

# Math and thread libraries must often come after all others when linking...
if (HAVE_LIBM)
    list(APPEND GMX_PUBLIC_LIBRARIES m)
endif()

option(GMX_NACL "Configure for Native Client builds" OFF)
if (GMX_NACL)
  list(APPEND GMX_EXTRA_LIBRARIES nosys)
  set(GMX_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lnosys")
  # TODO: Is this still necessary with the check for its presence?
  set(GMX_USE_NICE 0)
  set(GMX_NO_RENAME 1)
endif()
mark_as_advanced(GMX_NACL)

if(GMX_FAHCORE)
  set(COREWRAP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../corewrap" CACHE STRING
      "Path to swindirect.h")
  include_directories(${COREWRAP_INCLUDE_DIR})
endif()

option(GMX_BUILD_HELP "Build completions (requires that compiled binaries can be executed on build host) and install man pages if built (requires building the 'man' target manually)" OFF)
mark_as_advanced(GMX_BUILD_HELP)
if (GMX_BUILD_HELP AND SOURCE_IS_SOURCE_DISTRIBUTION AND BUILD_IS_INSOURCE)
    message(FATAL_ERROR
        "Rebuilding shell completions or man pages is not supported for "
        "in-source builds from a source distribution. "
        "Set GMX_BUILD_HELP=OFF or do an out-of-source build to proceed.")
endif()

if (GMX_BUILD_FOR_COVERAGE)
    # Set flags for coverage build here instead having to do so manually
    set(CMAKE_C_FLAGS "-g -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "-g -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")
endif()

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (GMX_LIB_MPI)
    # TODO: Restrict the scope of MPI dependence.
    # Targets that actually need MPI headers and build tool flags should
    # manage their own `target_link_libraries` locally. Such a change is beyond
    # the scope of the bug fix for #4678.
    set(_gmx_mpi_cxx_link_flags "${MPI_CXX_LINK_FLAGS}")
    set(_gmx_mpi_c_compile_flags "${MPI_C_COMPILE_OPTIONS};${MPI_C_COMPILE_DEFINITIONS}")
    set(_gmx_mpi_cxx_compile_flags "${MPI_CXX_COMPILE_OPTIONS};${MPI_CXX_COMPILE_DEFINITIONS}")
endif ()
if (NOT GMX_SKIP_DEFAULT_CFLAGS)
    #TODO(#3672): Use target_link_libraries(... MPI::MPI_CXX) instead of ${MPI_CXX_LINK_FLAGS}
    set(CMAKE_EXE_LINKER_FLAGS "${FFT_LINKER_FLAGS} ${_gmx_mpi_cxx_link_flags} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${FFT_LINKER_FLAGS} ${_gmx_mpi_cxx_link_flags} ${CMAKE_SHARED_LINKER_FLAGS}")
else()
    message("Recommended flags which are not added because GMX_SKIP_DEFAULT_CFLAGS=yes:")
    message("CMAKE_C_FLAGS: ${SIMD_C_FLAGS};${_gmx_mpi_c_compile_flags};${EXTRA_C_FLAGS};${GMXC_CFLAGS}")
    foreach(build_type ${build_types_with_explicit_flags})
        message("CMAKE_C_FLAGS_${build_type}: ${GMXC_CFLAGS_${build_type}}")
    endforeach()
    message("CMAKE_CXX_FLAGS: ${SIMD_CXX_FLAGS};${_gmx_mpi_cxx_compile_flags};${EXTRA_CXX_FLAGS};${GMXC_CXXFLAGS}")
    foreach(build_type ${build_types_with_explicit_flags})
        message("CMAKE_CXX_FLAGS_${build_type}: ${GMXC_CXXFLAGS_${build_type}}")
    endforeach()
    message("CMAKE_EXE_LINKER_FLAGS: ${FFT_LINKER_FLAGS} ${_gmx_mpi_cxx_link_flags}")
    message("CMAKE_SHARED_LINKER_FLAGS: ${FFT_LINKER_FLAGS} ${_gmx_mpi_cxx_link_flags}")
endif()
# Allow `admin` directory to be easily conveyed to nested CMake commands.
set(GMX_ADMIN_DIR ${CMAKE_SOURCE_DIR}/admin)

################################################################
# Shared library load path settings
################################################################
if(NOT GMX_BUILD_SHARED_EXE)
    # No rpath
    set(CMAKE_SKIP_RPATH TRUE)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS) # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
else()
    # The build folder always has bin/ and lib/; if we are also going to
    # install to lib/, then the installation RPATH works also in the build
    # tree.  This makes installation slightly faster (no need to rewrite the
    # RPATHs), and makes the binaries in the build tree relocatable.
    if(CMAKE_INSTALL_LIBDIR STREQUAL "lib")
        set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
        set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR TRUE)
    endif()
    # Set the RPATH as relative to the executable location to make the
    # binaries relocatable.
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin") #Assume OS X >=10.5
        set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    else()
        set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    endif()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

#COPYING file: Only necessary for binary distributions.
#Simpler to always install.
install(FILES COPYING DESTINATION ${GMX_INSTALL_GMXDATADIR} COMPONENT data)

if (BUILD_TESTING)
    include(tests/CheckTarget.cmake)
endif()

option(GMX_PYTHON_PACKAGE
       "Configure gmxapi Python package for use in build tree. Requires pybind11 installed for project Python interpreter."
       OFF)
mark_as_advanced(GMX_PYTHON_PACKAGE)

find_package(ImageMagick QUIET COMPONENTS convert)
mark_as_advanced(ImageMagick_EXECUTABLE_DIR)
include(gmxTestImageMagick)
GMX_TEST_IMAGEMAGICK(IMAGE_CONVERT_POSSIBLE)
add_subdirectory(share)
add_subdirectory(scripts)
add_subdirectory(api)
add_subdirectory(src)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

if(GMX_PYTHON_PACKAGE)
    add_subdirectory(python_packaging)
endif()

add_subdirectory(docs)

gmx_cpack_write_config()

#######################
## uninstall target
#######################
CONFIGURE_FILE(   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
                  IMMEDIATE @ONLY)
###########################
ADD_CUSTOM_TARGET(uninstall
                  "${CMAKE_COMMAND}" -P
                  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")
###########################
set_directory_properties(PROPERTIES
            ADDITIONAL_MAKE_CLEAN_FILES "install_manifest.txt")
