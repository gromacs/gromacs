#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

cmake_minimum_required(VERSION 2.8.8)
# When we require cmake >= 2.8.12, it will provide
# CMAKE_MINIMUM_REQUIRED_VERSION automatically, but in the meantime we
# need to set a variable, and it must have a different name.
set(GMX_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform)

project(Gromacs)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set up common version variables, as well as general information about
# the build tree (whether the build is from a source package or from a git
# repository).  Also declares a few functions that will be used for generating
# version info files later.
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
include(gmxBuildTypeReleaseWithAssert)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Reference RelWithAssert Profile." FORCE)
    # There's no need to offer a user the choice of ThreadSanitizer
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo" "Reference" "RelWithAssert" "Profile")
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

# Set a default valgrind suppression file.
# This unfortunately needs to duplicate information from CTest to work as
# expected...
set(MEMORYCHECK_SUPPRESSIONS_FILE
    "${CMAKE_SOURCE_DIR}/cmake/legacy_and_external.supp"
    CACHE FILEPATH
    "File that contains suppressions for the memory checker")
include(CTest)

include(gmxCPackUtilities)
gmx_cpack_init()

# Variables that accumulate stuff influencing the installed headers
set(INSTALLED_HEADER_INCLUDE_DIRS "")
set(INSTALLED_HEADER_DEFINITIONS "")

########################################################################
# Check and warn if cache generated on a different host is being reused
########################################################################
if(CMAKE_HOST_UNIX)
    execute_process(COMMAND hostname
                    OUTPUT_VARIABLE TMP_HOSTNAME
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(GMX_BUILD_HOSTNAME AND NOT "${GMX_BUILD_HOSTNAME}" STREQUAL "${TMP_HOSTNAME}")
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

set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")

if(GMX_TARGET_FUJITSU_SPARC64)
    # Fujitsu only has SIMD in double precision, so this will be faster
    set(GMX_DOUBLE_DEFAULT ON)
else()
    set(GMX_DOUBLE_DEFAULT OFF)
endif()
option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" ${GMX_DOUBLE_DEFAULT})
option(GMX_RELAXED_DOUBLE_PRECISION "Accept single precision 1/sqrt(x) when using Fujitsu HPC-ACE SIMD" OFF)
mark_as_advanced(GMX_RELAXED_DOUBLE_PRECISION)

option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREAD_MPI  "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" ON)
gmx_dependent_option(
    GMX_MPI_IN_PLACE
    "Enable MPI_IN_PLACE for MPIs that have it defined"
    ON
    GMX_MPI)
mark_as_advanced(GMX_MPI_IN_PLACE)
option(GMX_SOFTWARE_INVSQRT "Use GROMACS software 1/sqrt" ON)
mark_as_advanced(GMX_SOFTWARE_INVSQRT)
option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
mark_as_advanced(GMX_FAHCORE)

option(GMX_COOL_QUOTES "Enable GROMACS cool quotes" ON)
mark_as_advanced(GMX_COOL_QUOTES)
gmx_add_cache_dependency(GMX_COOL_QUOTES BOOL "NOT GMX_FAHCORE" OFF)

option(GMX_USE_OPENCL "Enable OpenCL acceleration" OFF)

# Decide on GPU settings based on user-settings and GPU/CUDA detection.
# We support CUDA >=v4.0 on *nix, but <= v4.1 doesn't work with MSVC
if(MSVC)
    set(REQUIRED_CUDA_VERSION 4.1)
else()
    set(REQUIRED_CUDA_VERSION 4.0)
endif()
set(REQUIRED_CUDA_COMPUTE_CAPABILITY 2.0)

# OpenCL required version: 1.1 or newer
set(REQUIRED_OPENCL_MIN_VERSION 1.1)

if(NOT GMX_USE_OPENCL)
    # CUDA detection is done only if GMX_USE_OPENCL is OFF
    include(gmxManageGPU)
else()
    #Now the OpenCL path
    if(GMX_GPU)
        include(gmxManageOpenCL)
    else(GMX_GPU)
        message(FATAL_ERROR "OpenCL requested but GPU option is not enabled (try -DGMX_GPU=on) ")
    endif(GMX_GPU)
endif()

include(gmxDetectSimd)
gmx_detect_simd(GMX_SUGGESTED_SIMD)

gmx_option_multichoice(
    GMX_SIMD
    "SIMD instruction set for CPU kernels and compiler optimization"
    "${GMX_SUGGESTED_SIMD}"
    None SSE2 SSE4.1 AVX_128_FMA AVX_256 AVX2_256 AVX_512F AVX_512ER MIC ARM_NEON ARM_NEON_ASIMD IBM_QPX IBM_VMX IBM_VSX Sparc64_HPC_ACE Reference)

if(GMX_TARGET_MIC)
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

gmx_option_multichoice(
    GMX_QMMM_PROGRAM
    "QM package for QM/MM"
    None
    none gaussian mopac gamess orca)

gmx_dependent_cache_variable(GMX_SIMD_REF_FLOAT_WIDTH  "Reference SIMD single precision width" STRING "4" "GMX_SIMD STREQUAL REFERENCE")
gmx_dependent_cache_variable(GMX_SIMD_REF_DOUBLE_WIDTH "Reference SIMD double precision width" STRING "2" "GMX_SIMD STREQUAL REFERENCE")

# This should be moved to a separate NBNXN cmake module when that code is cleaned up and modularized

option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
mark_as_advanced(GMX_BROKEN_CALLOC)
option(GMX_LOAD_PLUGINS "Compile with plugin support, needed to read VMD supported file formats" ON)
mark_as_advanced(GMX_LOAD_PLUGINS)

option(GMX_GPU  "Enable GPU acceleration" ON)
option(GMX_OPENMP "Enable OpenMP-based multithreading" ON)

option(GMX_USE_TNG "Use the TNG library for trajectory I/O" ON)

option(GMX_BUILD_MDRUN_ONLY "Build and install only the mdrun binary" OFF)

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

option(GMX_BUILD_SHARED_EXE
    "Build exectuables as shared binaries. If not set, this disables rpath and dynamic linker flags in an attempt to build a static binary, but this may require setting up the toolchain properly and making appropriate libraries available."
    ON)
mark_as_advanced(GMX_BUILD_SHARED_EXE)

######################################################################
# Detect OpenMP support
######################################################################
# The OpenMP detection _must_ come before tests for other CFLAGS.
include(gmxManageOpenMP)



######################################################################
# Compiler tests
# These need to be done early (before further tests).
#####################################################################

# The cmake/Check{C,CXX}CompilerFlag.cmake files in the GROMACS distribution
# are used with permission from CMake v3.0.0 so that GROMACS can detect
# invalid options with the Intel Compilers, and we have added a line
# to detect warnings with the Fujitsu compilers on K computer and ICC.
# CMake-3.0 also has a bug where the FAIL_REGEX pattern for AIX contains
# a semicolon. Since this is also used as a separator in lists inside CMake,
# that string ends up being split into two separate patterns, and the last
# part is just a single word that also matches other messages. We solved this
# by replacing the semicolon with a period that matches any character.
#
# These files should be removed from the source tree when a CMake version that
# includes the features in question becomes required for building GROMACS.
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

include(gmxCFlags)
gmx_c_flags()

# This variable should be used for additional compiler flags which are not
# generated in gmxCFlags nor are SIMD or MPI related.
set(EXTRA_C_FLAGS "")
set(EXTRA_CXX_FLAGS "")

# Run through a number of tests for buggy compilers and other issues
include(gmxTestCompilerProblems)
gmx_test_compiler_problems()
# GMX_SIMD will not be set automatically until the second
# pass (which is not strictly guaranteed to occur), so putting this
# check here among logically-related tests is inefficient, but the
# potential loss is likely zero.
if(GMX_SIMD STREQUAL "AVX_256"
        AND CMAKE_COMPILER_IS_GNUCC
        AND (C_COMPILER_VERSION VERSION_EQUAL "4.6.1"
            OR CXX_COMPILER_VERSION VERSION_EQUAL "4.6.1"))
    message(FATAL_ERROR "gcc 4.6.1 has buggy support for AVX, and GROMACS mdrun will not work. If you want simulation performance, use a more recent compiler. Otherwise, use GMX_SIMD=SSE4.1")
    # See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=49002
endif()


if(GMX_DOUBLE)
    add_definitions(-DGMX_DOUBLE)
    list(APPEND INSTALLED_HEADER_DEFINITIONS "-DGMX_DOUBLE")
    if(GMX_RELAXED_DOUBLE_PRECISION)
        add_definitions(-DGMX_RELAXED_DOUBLE_PRECISION)
    endif()
endif()
if(GMX_SOFTWARE_INVSQRT)
    list(APPEND INSTALLED_HEADER_DEFINITIONS "-DGMX_SOFTWARE_INVSQRT")
endif()

if(WIN32)
    list(APPEND GMX_EXTRA_LIBRARIES "wsock32")
    add_definitions(-DGMX_HAVE_WINSOCK)
endif()



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

check_include_files(regex.h      HAVE_POSIX_REGEX)
# TODO: It could be nice to inform the user if no regex support is found,
# as selections won't be fully functional.

include(CheckCXXSymbolExists)
check_cxx_symbol_exists(posix_memalign    stdlib.h     HAVE_POSIX_MEMALIGN)
check_cxx_symbol_exists(memalign          stdlib.h     HAVE_MEMALIGN)
check_cxx_symbol_exists(_aligned_malloc   stdlib.h     HAVE__ALIGNED_MALLOC)
check_cxx_symbol_exists(gettimeofday      sys/time.h   HAVE_GETTIMEOFDAY)
check_cxx_symbol_exists(sysconf           unistd.h     HAVE_SYSCONF)
check_cxx_symbol_exists(nice              unistd.h     HAVE_NICE)
check_cxx_symbol_exists(fsync             unistd.h     HAVE_FSYNC)
check_cxx_symbol_exists(_fileno           stdio.h      HAVE__FILENO)
check_cxx_symbol_exists(fileno            stdio.h      HAVE_FILENO)
check_cxx_symbol_exists(_commit           io.h         HAVE__COMMIT)
check_cxx_symbol_exists(sigaction         signal.h     HAVE_SIGACTION)
check_cxx_symbol_exists(rsqrt             math.h       HAVE_RSQRT)
check_cxx_symbol_exists(rsqrtf            math.h       HAVE_RSQRTF)
check_cxx_symbol_exists(sqrtf             math.h       HAVE_SQRTF)

include(CheckLibraryExists)
find_library(HAVE_LIBM m)
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)
check_library_exists(m feenableexcept "" HAVE_FEENABLEEXCEPT)

include(TestSchedAffinity)
test_sched_affinity(HAVE_SCHED_AFFINITY)

include(TestBigEndian)
test_big_endian(GMX_INTEGER_BIG_ENDIAN)

set(GMX_USE_NICE 0)
if (HAVE_UNISTD_H AND HAVE_NICE)
    set(GMX_USE_NICE 1)
endif()

# Management of GROMACS options for specific toolchains should go
# here. Because the initial settings for some of the main options have
# already happened, but things like library detection and MPI compiler
# feature detection have not, the docstrings for any over-rides of
# GROMACS defaults or user settings will make sense. Also, any
# toolchain-related reasons for choosing whether to detect various
# things can be sorted out now, before the detection takes place.
if(${CMAKE_SYSTEM_NAME} MATCHES BlueGene)
    include(gmxManageBlueGene)
endif()

if(GMX_TARGET_FUJITSU_SPARC64)
    include(gmxManageFujitsuSparc64)
endif()

########################################################################
#Process MPI settings
########################################################################
include(gmxManageMPI)


########################################################################
#Process shared/static library settings
########################################################################
include(gmxManageSharedLibraries)


########################################################################
# Find external packages                                               #
########################################################################

# TNG wants zlib if it is available. And static libxml2 might have a dependency
find_package(ZLIB QUIET)
include(gmxTestZLib)
gmx_test_zlib(HAVE_ZLIB)

# Unconditionally find the package, as it is also required for unit
# tests. This exports LIBXML2_FOUND, which we should not use because
# it does not tell us that linking will succeed. Instead, we test that
# next.
if(DEFINED LIBXML2_LIBRARIES)
  set(LibXml2_FIND_QUIETLY TRUE)
endif()
find_package(LibXml2)
include(gmxTestLibXml2)
gmx_test_libxml2(HAVE_LIBXML2)
option(GMX_XML "Use libxml2 to parse xml files (currently has no effect)" ${HAVE_LIBXML2})
set(PKG_XML "")
mark_as_advanced(GMX_XML)
# Don't actually do anything, since libxml2 is currently not used by libgromacs
#if(GMX_XML AND NOT HAVE_LIBXML2)
#    message(FATAL_ERROR "libxml2 not found. Set GMX_XML=OFF to compile without XML support")
#endif()
#if(GMX_XML)
#    include_directories(${LIBXML2_INCLUDE_DIR})
#    set(PKG_XML libxml-2.0)
#    set(XML_LIBRARIES ${LIBXML2_LIBRARIES})
#endif()

option(GMX_EXTRAE "Add support for tracing using EXTRAE" OFF)
mark_as_advanced(GMX_EXTRAE)

if (GMX_EXTRAE)
  find_package(EXTRAE)
  if(EXTRAE_FOUND)
    include_directories(${EXTRAE_INCLUDE_DIR})
    set(HAVE_EXTRAE 1)
  else()
    message(FATAL_ERROR "EXTRAE library was not found. Please add the correct path to CMAKE_PREFIX_PATH")
  endif()
endif()

option(GMX_X11 "Use X window system" OFF)
if (GMX_X11)
    find_package(X11)
    # X11 includes/libraries are only set in the ngmx subdirectory!
    if(NOT X11_FOUND)
        message(FATAL_ERROR
                "X11 include files and/or libraries were not found. "
                "Set GMX_X11=OFF to compile without X11 support. "
                "gmx view will not be available.")
    endif()
    include_directories(${X11_INCLUDE_DIR})
endif()

include(ThreadMPI)
# Enable core threading facilities
tmpi_enable_core("${CMAKE_SOURCE_DIR}/src/external/thread_mpi/include")
# Enable tMPI C++ support
tmpi_enable_cxx()
if(GMX_THREAD_MPI)
    # enable MPI functions
    tmpi_enable()
    set(MPI_IN_PLACE_EXISTS 1)
endif()
# If atomics are manually disabled a define is needed because atomics.h doesn't depend on config.h
if (TMPI_ATOMICS_DISABLED)
   add_definitions(-DTMPI_ATOMICS_DISABLED)
endif()

if(GMX_GPU)
    # now that we have detected the dependencies, do the second configure pass
    gmx_gpu_setup()
endif()

if(CYGWIN)
    set(GMX_CYGWIN 1)
endif()

if(WIN32)
    set(GMX_NATIVE_WINDOWS 1)
    # This makes windows.h not declare min/max as macros that would break
    # C++ code using std::min/std::max.
    add_definitions(-DNOMINMAX)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin") #Work-around for cmake bug #10837
    if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID MATCHES "Intel" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")
    endif()
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
    endif()
endif()

# Detect boost unless GMX_EXTERNAL_BOOST is explicitly OFF
# Used for default if GMX_EXTERNAL_BOOST is not defined (first CMake pass)
if(NOT DEFINED GMX_EXTERNAL_BOOST OR GMX_EXTERNAL_BOOST)
    if(DEFINED Boost_INCLUDE_DIR)
        set(Boost_FIND_QUIETLY TRUE)
    endif()
    find_package(Boost 1.44.0)
    if(Boost_FOUND AND Boost_VERSION VERSION_LESS "104400")
        set(Boost_FOUND FALSE)
    endif()
    # Print the notification only on first run, when determining the default
    if(NOT DEFINED GMX_EXTERNAL_BOOST AND NOT Boost_FOUND)
        message("Boost >= 1.44 not found. Using minimal internal version. "
                "This may cause trouble if you plan on compiling/linking other "
                "software that uses Boost against GROMACS.")
    endif()
endif()
option(GMX_EXTERNAL_BOOST "Use external Boost instead of minimal built-in version"
       ${Boost_FOUND})
if(GMX_EXTERNAL_BOOST AND NOT Boost_FOUND)
    message(FATAL_ERROR
        "Boost >= 1.44 not found. "
        "You can set GMX_EXTERNAL_BOOST=OFF to compile against minimal "
        "version of Boost included with GROMACS.")
endif()

if(NOT DEFINED GMX_BUILD_UNITTESTS AND NOT HAVE_LIBXML2)
    message(WARNING "libxml2 not found. Will build GROMACS without unit-tests. This is not recommended, because the unit-tests help to verify that GROMACS functions correctly. Most likely you are missing the libxml2-dev(el) package. After you installed it, set GMX_BUILD_UNITTESTS=ON.")
endif()
option(GMX_BUILD_UNITTESTS "Build unit tests with BUILD_TESTING (uses Google C++ Testing and Mocking Frameworks, requires libxml2)" ${HAVE_LIBXML2})
mark_as_advanced(GMX_BUILD_UNITTESTS)
gmx_add_cache_dependency(GMX_BUILD_UNITTESTS BOOL BUILD_TESTING OFF)
if (GMX_BUILD_UNITTESTS AND NOT HAVE_LIBXML2)
    message(FATAL_ERROR
        "Cannot build unit tests without libxml2. "
        "Either set GMX_BUILD_UNITTESTS=OFF or tell CMake how to find a working version of libxml2.")
endif()

########################################################################
# Our own GROMACS tests
########################################################################

add_definitions( -DHAVE_CONFIG_H )
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/src)
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/src/external/thread_mpi/include)
# Required for config.h, maybe should only be set in src/CMakeLists.txt
include_directories(BEFORE ${CMAKE_BINARY_DIR}/src)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

include(gmxSetBuildInformation)
gmx_set_build_information()
# Turn on RDTSCP if:
# - the build system's CPU supports it
# - the acceleration is set to AVX as all AVX-capable CPUs support AVX (which
#   at this point means that the user set it).
# Note: it's better to not use the later set value of GMX_SIMD because
# it reflects the system's capability of both compiling and running AVX code.
# TODO: After merge with 5.0 one could implement a cache variable dependency
# such that GMX_USE_RDTSCP can change if GMX_SIMD is changed to AVX
# after the first cmake pass.
if (BUILD_CPU_FEATURES MATCHES "rdtscp" OR GMX_SIMD MATCHES "AVX")
    set(GMX_USE_RDTSCP_DEFAULT_VALUE ON)
else()
    set(GMX_USE_RDTSCP_DEFAULT_VALUE OFF)
endif()
option(GMX_USE_RDTSCP "Use RDTSCP for better CPU-based timers (available on recent x86 CPUs; might need to be off when compiling for heterogeneous environments)" ${GMX_USE_RDTSCP_DEFAULT_VALUE})
mark_as_advanced(GMX_USE_RDTSCP)
if(GMX_USE_RDTSCP)
    set(HAVE_RDTSCP 1)
endif()

include(gmxTestFloatFormat)
gmx_test_float_format(GMX_FLOAT_FORMAT_IEEE754
                      GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER
                      GMX_IEEE754_BIG_ENDIAN_WORD_ORDER)

include(gmxTestLargeFiles)
gmx_test_large_files(GMX_LARGEFILES)

include(gmxTestSignal)
gmx_test_sigusr1(HAVE_SIGUSR1)

include(gmxTestPipes)
gmx_test_pipes(HAVE_PIPES)

include(gmxTestIsfinite)
gmx_test_isfinite(HAVE_ISFINITE)
gmx_test__isfinite(HAVE__ISFINITE)
gmx_test__finite(HAVE__FINITE)

include(gmxTestCXX11)
gmx_test_cxx11(GMX_CXX11_SUPPORTED GMX_CXX11_FLAGS)
include(CMakeDependentOption)
# clang with some older versions of boost fail to work with C++11
if(CMAKE_CXX_COMPILER_ID MATCHES Clang AND
    (Boost_VERSION VERSION_EQUAL 104601 OR Boost_VERSION VERSION_EQUAL 104700))
    set(CLANG_AND_OLD_EXTERNAL_BOOST TRUE)
endif()
# nvcc does not support C++11 flags, so with GPUs we prefer to skip C++11 flags
# entirely to keep the compilation environment uniform.
cmake_dependent_option(GMX_CXX11
    "Use C++11 features"
    ON "GMX_CXX11_SUPPORTED AND NOT GMX_GPU AND NOT CLANG_AND_OLD_EXTERNAL_BOOST" OFF)
mark_as_advanced(GMX_CXX11)
if(GMX_CXX11)
    set(EXTRA_CXX_FLAGS "${EXTRA_CXX_FLAGS} ${GMX_CXX11_FLAGS}")
endif()

# Now we can test for CXX11_REGEX include file
check_include_file_cxx(regex     HAVE_CXX11_REGEX)

include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
if(NOT GMX_SYSTEM_XDR)
    set(GMX_INTERNAL_XDR 1)
endif()


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

# Process QM/MM Settings
if(${GMX_QMMM_PROGRAM} STREQUAL "GAUSSIAN")
    set(GMX_QMMM_GAUSSIAN 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "MOPAC")
    set(GMX_QMMM_MOPAC 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "GAMESS")
    set(GMX_QMMM_GAMESS 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "ORCA")
    set(GMX_QMMM_ORCA 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "NONE")
    # nothing to do
else()
    gmx_invalid_option_value(GMX_QMMM_PROGRAM)
endif()


##################################################
# Process FFT library settings
##################################################
include(gmxManageFFTLibraries)


include(gmxManageLinearAlgebraLibraries)

# Whether GROMACS will really try to compile support for VMD plugins
set(GMX_USE_PLUGINS OFF)

if(GMX_LOAD_PLUGINS)
  if(NOT WIN32)
    # Native Windows does not have, nor need dlopen
    include(gmxTestdlopen)
    gmx_test_dlopen(HAVE_DLOPEN)
  endif()

  # so, should we use plug-ins?
  if(WIN32 OR (HAVE_DLOPEN AND BUILD_SHARED_LIBS))
    if(NOT VMD_QUIETLY)
      MESSAGE(STATUS "Using dynamic plugins (e.g VMD-supported file formats)")
    endif()
    if(NOT GMX_VMD_PLUGIN_PATH)
      find_package(VMD)
    endif()
    set(GMX_USE_PLUGINS ON)
    list(APPEND GMX_EXTRA_LIBRARIES ${CMAKE_DL_LIBS}) # magic cross-platform pre-set variable for dlopen library
    set(PKG_DL_LIBS "-l${CMAKE_DL_LIBS}")
  else()
    set(PKG_DL_LIBS)
  endif()
endif()
set(VMD_QUIETLY TRUE CACHE INTERNAL "")

# Link real-time library for POSIX timers. The check for clock_gettime
# confirms the linkability of rt.
if(HAVE_TIME_H AND HAVE_UNISTD_H AND HAVE_CLOCK_GETTIME)
    list(APPEND GMX_EXTRA_LIBRARIES rt)
endif()

# Math and thread libraries must often come after all others when linking...
if(HAVE_LIBM)
    list(APPEND GMX_EXTRA_LIBRARIES m)
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

# Value of GMX_BUILD_HELP=AUTO tries to generate things, but will only
# produce warnings if that fails.
set(build_help_default AUTO)
if (SOURCE_IS_SOURCE_DISTRIBUTION OR CMAKE_CROSSCOMPILING)
    set(build_help_default OFF)
endif()
gmx_option_trivalue(GMX_BUILD_HELP "Build completions automatically (requires that compiled binaries can be executed on the build host) and install man pages if built (requires building the 'man' target manually)" ${build_help_default})
mark_as_advanced(GMX_BUILD_HELP)
if (GMX_BUILD_HELP AND SOURCE_IS_SOURCE_DISTRIBUTION AND BUILD_IS_INSOURCE)
    message(FATAL_ERROR
        "Rebuilding shell completions or man pages is not supported for "
        "in-source builds from a source distribution. "
        "Set GMX_BUILD_HELP=OFF or do an out-of-source build to proceed.")
endif()

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (NOT GMX_SKIP_DEFAULT_CFLAGS)
    set(CMAKE_C_FLAGS "${SIMD_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${SIMD_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${FFT_LINKER_FLAGS} ${MPI_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${FFT_LINKER_FLAGS} ${MPI_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
else()
    message("Recommended flags which are not added because GMX_SKIP_DEFAULT_CFLAGS=yes:")
    message("CMAKE_C_FLAGS: ${SIMD_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${GMXC_CFLAGS}")
    message("CMAKE_C_FLAGS_RELEASE: ${GMXC_CFLAGS_RELEASE}")
    message("CMAKE_C_FLAGS_DEBUG: ${GMXC_CFLAGS_DEBUG}")
    message("CMAKE_CXX_FLAGS: ${SIMD_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_CXX_FLAGS} ${GMXC_CXXFLAGS}")
    message("CMAKE_CXX_FLAGS_RELEASE: ${GMXC_CXXFLAGS_RELEASE}")
    message("CMAKE_CXX_FLAGS_DEBUG: ${GMXC_CXXFLAGS_DEBUG}")
    message("CMAKE_EXE_LINKER_FLAGS: ${FFT_LINKER_FLAGS} ${MPI_LINKER_FLAGS}")
    message("CMAKE_SHARED_LINKER_FLAGS: ${MPI_LINKER_FLAGS}")
endif()

if(NOT GMX_OPENMP)
    #Unset all OpenMP flags in case OpenMP was disabled either by the user
    #or because it was only partially detected (e.g. only for C but not C++ compiler)
    unset(OpenMP_C_FLAGS CACHE)
    unset(OpenMP_CXX_FLAGS CACHE)
else()
    set(GMX_EXE_LINKER_FLAGS ${GMX_EXE_LINKER_FLAGS} ${OpenMP_LINKER_FLAGS})
    set(GMX_SHARED_LINKER_FLAGS ${GMX_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS})
endif()

########################################################################
# Specify install locations
########################################################################
# Use GNUInstallDirs to set paths on multiarch systems.
include(GNUInstallDirs)

# Customization for the installation tree paths.
set(GMX_LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
    "Library installation directory (default: ${CMAKE_INSTALL_LIBDIR})")
set(GMX_DATA_INSTALL_DIR gromacs CACHE STRING
    "Data installation directory under share/ (default: gromacs)")
mark_as_advanced(GMX_LIB_INSTALL_DIR GMX_DATA_INSTALL_DIR)

# These variables are used internally to provide a central location for
# customizing the install locations.
set(LIB_INSTALL_DIR       ${GMX_LIB_INSTALL_DIR})
set(BIN_INSTALL_DIR       bin)
# This variable also gets written into config.h for use in finding the data
# directories.
set(DATA_INSTALL_DIR      share/${GMX_DATA_INSTALL_DIR})
set(MAN_INSTALL_DIR       share/man)
# If the nesting level wrt. the installation root is changed,
# gromacs-config.cmake.cmakein needs to be adapted.
set(CMAKE_INSTALL_DIR     share/cmake)
# TODO: Make GMXRC adapt if this is changed
set(PKGCONFIG_INSTALL_DIR ${LIB_INSTALL_DIR}/pkgconfig)
set(OCL_INSTALL_DIR       ${DATA_INSTALL_DIR}/opencl)
set(INCL_INSTALL_DIR      include)

list(APPEND INSTALLED_HEADER_INCLUDE_DIRS ${INCL_INSTALL_DIR})

# Binary and library suffix options
include(gmxManageSuffixes)

################################################################
# Shared library load path settings
################################################################
# CMake supports RPATH on OS X only from 2.8.12 upwards.
# CMAKE_SYSTEM_VERSION > 8.0 matches OS X 10.5 and above, where RPATH support
# was added.

if(NOT GMX_BUILD_SHARED_EXE)
    # No rpath
    set(CMAKE_SKIP_RPATH TRUE)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS) # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
elseif((NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin") OR
   ((CMAKE_SYSTEM_VERSION VERSION_GREATER 8.0) AND (NOT CMAKE_VERSION VERSION_LESS 2.8.12)))
    # The build folder always has bin/ and lib/; if we are also going to
    # install to lib/, then the installation RPATH works also in the build
    # tree.  This makes installation slightly faster (no need to rewrite the
    # RPATHs), and makes the binaries in the build tree relocatable.
    if(GMX_LIB_INSTALL_DIR STREQUAL "lib")
        set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    endif()
    # Set the RPATH as relative to the executable location to make the
    # binaries relocatable.
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${GMX_LIB_INSTALL_DIR}")
    else()
        set(CMAKE_INSTALL_RPATH "@executable_path/../${GMX_LIB_INSTALL_DIR}")
    endif()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_MACOSX_RPATH 1)
else()
    # We are on Darwin/OSX, and CMake cannot handle RPATHs automatically.
    if(CMAKE_SYSTEM_VERSION VERSION_GREATER 8.0)
        # Set the RPATH options manually.
        set(CMAKE_INSTALL_NAME_DIR "@rpath")
        set(GMX_EXE_LINKER_FLAGS ${GMX_EXE_LINKER_FLAGS} "-Wl,-rpath,@executable_path/../${GMX_LIB_INSTALL_DIR}")
    else()
        # Use the old INSTALL_NAME_DIR mechanism if RPATH is not supported.
        set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
    endif()
endif()

#COPYING file: Only necessary for binary distributions.
#Simpler to always install.
install(FILES COPYING DESTINATION ${DATA_INSTALL_DIR} COMPONENT data)

if(GMX_EXTERNAL_BOOST)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    list(APPEND INSTALLED_HEADER_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
else()
    include_directories(BEFORE SYSTEM ${CMAKE_SOURCE_DIR}/src/external/boost)
    list(APPEND INSTALLED_HEADER_INCLUDE_DIRS ${INCL_INSTALL_DIR}/gromacs/external/boost)
    list(APPEND INSTALLED_HEADED_DEFINITIONS "-DBOOST_NO_TYPEID")
    # typeid not supported for minimal internal version
    # (would add significant amount of code)
    add_definitions(-DBOOST_NO_TYPEID)
    if (NOT GMX_BUILD_MDRUN_ONLY)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/external/boost/boost
                DESTINATION ${INCL_INSTALL_DIR}/gromacs/external/boost
                COMPONENT development)
    endif()
endif()

if (GMX_BUILD_FOR_COVERAGE)
    # Code heavy with asserts makes conditional coverage close to useless metric,
    # as by design most of the false branches are impossible to trigger in
    # correctly functioning code.  And the benefit of testing those that could
    # be triggered by using an API against its specification isn't usually
    # worth the effort.
    add_definitions(-DNDEBUG -DBOOST_DISABLE_ASSERTS -DGMX_DISABLE_ASSERTS)
endif()

if (BUILD_TESTING)
    # "tests" target builds all the separate test binaries.
    add_custom_target(tests)
    # "run-ctest" is an internal target that actually runs the tests.
    # This is necessary to be able to add separate targets that execute as part
    # of 'make check', but are ensured to be executed after the actual tests.
    add_custom_target(run-ctest
                      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                      COMMENT "Running all tests"
                      VERBATIM)
    add_dependencies(run-ctest tests)
    # "check" target builds and runs all tests.
    add_custom_target(check DEPENDS run-ctest)
endif()

if (NOT GMX_BUILD_MDRUN_ONLY)
    add_subdirectory(docs)
    add_subdirectory(share)
    add_subdirectory(scripts)
endif()
add_subdirectory(src)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

gmx_cpack_write_config()

# Issue a warning if NVIDIA GPUs were detected, but CUDA was not found.
# Don't bother the user after the first configure pass.
if ((CUDA_NOTFOUND_AUTO AND GMX_DETECT_GPU_AVAILABLE) AND NOT GMX_GPU_DETECTION_DONE)
    message(WARNING "${CUDA_NOTFOUND_MESSAGE}")
endif()
set(GMX_GPU_DETECTION_DONE TRUE CACHE INTERNAL "Whether GPU detection has already been done")

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
