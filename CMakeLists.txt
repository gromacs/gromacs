#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform)

project(Gromacs)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# PROJECT_VERSION should have the following structure:
# VERSION-dev[-SUFFIX] where the VERSION should have the for: vMajor.vMinor.vPatch
#
# The "-dev" suffix is important to keep because it makes possible to distinguish
# between a build from official release and a build from git release branch on a
# machine with no git.
#
# NOTE: when releasing the "-dev" suffix needs to be stripped off!
set(PROJECT_VERSION "5.0-beta2-dev")
# The version number of the regressiontest tarball against which this
# git branch can be tested. Normally, this will be the version of the
# last patch release. Comment the next line out for branches leading
# to a major/minor release.
set(REGRESSIONTEST_VERSION "5.0-beta1")
set(CUSTOM_VERSION_STRING ""
    CACHE STRING "Custom version string (if empty, use hard-coded default)")
mark_as_advanced(CUSTOM_VERSION_STRING)
if (CUSTOM_VERSION_STRING)
    set(PROJECT_VERSION ${CUSTOM_VERSION_STRING})
endif (CUSTOM_VERSION_STRING)
set(LIBRARY_SOVERSION 0)
set(LIBRARY_VERSION ${LIBRARY_SOVERSION}.0.0)
# It is a bit irritating, but this has to be set separately for now!
SET(CPACK_PACKAGE_VERSION_MAJOR "5")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
#SET(CPACK_PACKAGE_VERSION_PATCH "0")

# The numerical gromacs version. It is 40600 for 4.6.0.
# The #define GMX_VERSION in gromacs/version.h is set to this value.
math(EXPR NUM_VERSION
    "${CPACK_PACKAGE_VERSION_MAJOR}*10000 + ${CPACK_PACKAGE_VERSION_MINOR}*100")
if(CPACK_PACKAGE_VERSION_PATCH)
    math(EXPR NUM_VERSION
         "${NUM_VERSION} + ${CPACK_PACKAGE_VERSION_PATCH}")
endif()

# The API version tracks the numerical Gromacs version (for now).
# It is potentially different from the Gromacs version in the future, if
# the programs/libraries diverge from the presumably more stable API.
# The #define GMX_API_VERSION in version.h is set to this value to
# provide backward compatibility of software written against the Gromacs API.
set(API_VERSION ${NUM_VERSION})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND UNIX)
    set(CMAKE_INSTALL_PREFIX "/usr/local/gromacs" CACHE STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif()

set(GMX_INSTALL_PREFIX "" CACHE STRING "Prefix gets appended to CMAKE_INSTALL_PREFIX. For cpack it sets the root folder of the archive.")
mark_as_advanced(GMX_INSTALL_PREFIX)

include(gmxBuildTypeReference)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Reference." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

enable_language(C)
enable_language(CXX)

set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_VENDOR "gromacs.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gromacs - a toolkit for high-performance molecular simulation")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/admin/InstallWelcome.txt")
# Its GPL/LGPL, so they do not have to agree to a license for mere usage, but some installers require this...
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/admin/InstallInfo.txt")
set(CPACK_SOURCE_IGNORE_FILES "\\\\.isreposource$;\\\\.git/;\\\\.gitignore$")
set(CPACK_PROJECT_CONFIG_FILE "${CMAKE_SOURCE_DIR}/CPackInit.cmake")
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_SOURCE_DIR};/;${CMAKE_BINARY_DIR}/share/man/man1;share/man/man1;${CMAKE_BINARY_DIR}/share/man/man7;share/man/man7;${CMAKE_BINARY_DIR}/share/html/final;share/html/final")
set(CPACK_PACKAGE_CONTACT "gmx-users@gromacs.org")

#must come after all cpack settings!
include(CPack)

# Set a default valgrind suppression file.
# This unfortunately needs to duplicate information from CTest to work as
# expected...
set(MEMORYCHECK_SUPPRESSIONS_FILE
    "${CMAKE_SOURCE_DIR}/cmake/legacy_and_external.supp"
    CACHE FILEPATH
    "File that contains suppressions for the memory checker")
include(CTest)

set(SOURCE_IS_GIT_REPOSITORY OFF)
set(SOURCE_IS_SOURCE_DISTRIBUTION OFF)
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(SOURCE_IS_GIT_REPOSITORY ON)
endif()
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/admin/.isreposource")
    set(SOURCE_IS_SOURCE_DISTRIBUTION ON)
endif()

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
            is being reused! This could lead to inconsitencies; therefore, it is
            recommended to regenerate the cache!")
    endif()
    set(GMX_BUILD_HOSTNAME "${TMP_HOSTNAME}" CACHE INTERNAL
            "Hostname of the machine where the cache was generated.")
endif()

########################################################################
# User input options                                                   #
########################################################################
include(gmxOptionUtilities)

set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")

option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" OFF)
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

option(GMX_COOL_QUOTES "Enable Gromacs cool quotes" ON)
mark_as_advanced(GMX_COOL_QUOTES)
gmx_add_cache_dependency(GMX_COOL_QUOTES BOOL "NOT GMX_FAHCORE" OFF)

# decide on GPU settings based on user-settings and GPU/CUDA detection
include(gmxManageGPU)

# Detect the architecture the compiler is targetting, detect
# acceleration possibilities on that hardware, suggest an acceleration
# to use if none is specified, and populate the cache option for CPU
# accleration.
include(gmxDetectTargetArchitecture)
gmx_detect_target_architecture()
include(gmxDetectAcceleration)
gmx_detect_acceleration(GMX_SUGGESTED_CPU_ACCELERATION)

gmx_option_multichoice(
    GMX_CPU_ACCELERATION
    "Acceleration for CPU kernels and compiler optimization"
    "${GMX_SUGGESTED_CPU_ACCELERATION}"
    None SSE2 SSE4.1 AVX_128_FMA AVX_256 IBM_QPX Sparc64_HPC_ACE)

gmx_option_multichoice(
    GMX_FFT_LIBRARY
    "FFT library"
    "fftw3"
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

option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
mark_as_advanced(GMX_BROKEN_CALLOC)
option(GMX_LOAD_PLUGINS "Compile with plugin support, needed to read VMD supported file formats" ON)
mark_as_advanced(GMX_LOAD_PLUGINS)

option(GMX_GPU  "Enable GPU acceleration" ON)
option(GMX_OPENMP "Enable OpenMP-based multithreading" ON)

option(GMX_USE_TNG "Use the TNG library for trajectory I/O" ON)

option(GMX_GIT_VERSION_INFO "Generate git version information" ${SOURCE_IS_GIT_REPOSITORY})
mark_as_advanced(GMX_GIT_VERSION_INFO)

if(UNIX)
    option(GMX_SYMLINK_OLD_BINARY_NAMES "Create symbolic links for pre-5.0 binary names" ON)
endif()
option(GMX_BUILD_MDRUN_ONLY "Build and install only the mdrun binary" OFF)

option(GMX_CYCLE_SUBCOUNTERS "Enable cycle subcounters to get a more detailed cycle timings" OFF)
mark_as_advanced(GMX_CYCLE_SUBCOUNTERS)

option(GMX_SKIP_DEFAULT_CFLAGS "Don't automatically add suggested/required Compiler flags." OFF)
mark_as_advanced(GMX_SKIP_DEFAULT_CFLAGS)

option(GMX_BUILD_FOR_COVERAGE
       "Tune build for better code coverage metrics (e.g., disable asserts)"
       OFF)
mark_as_advanced(GMX_BUILD_FOR_COVERAGE)


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
# are used with permission from CMake v2.8.9 so that GROMACS can detect
# invalid options with the Intel Compilers.
# These files should be removed from the source tree when a CMake version that
# includes the features in question becomes required for building GROMACS.
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)


include(gmxCFlags)
gmx_c_flags()

# This variable should be used for additional compiler flags which are not
# generated in gmxCFlags nor are acceleration or MPI related.
set(EXTRA_C_FLAGS "")
set(EXTRA_CXX_FLAGS "")

# Run through a number of tests for buggy compilers and other issues
include(gmxTestCompilerProblems)
gmx_test_compiler_problems()


set(PKG_CFLAGS "")
if(GMX_DOUBLE)
    add_definitions(-DGMX_DOUBLE)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_DOUBLE")
endif(GMX_DOUBLE)
if(GMX_SOFTWARE_INVSQRT)
  set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_SOFTWARE_INVSQRT")
endif(GMX_SOFTWARE_INVSQRT)



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
check_include_file_cxx(regex     HAVE_CXX11_REGEX)
# TODO: It could be nice to inform the user if no regex support is found,
# as selections won't be fully functional.

include(CheckFunctionExists)
check_function_exists(posix_memalign    HAVE_POSIX_MEMALIGN)
check_function_exists(memalign          HAVE_MEMALIGN)
check_function_exists(_aligned_malloc   HAVE__ALIGNED_MALLOC)
check_function_exists(gettimeofday      HAVE_GETTIMEOFDAY)
check_function_exists(fsync             HAVE_FSYNC)
check_function_exists(_fileno           HAVE__FILENO)
check_function_exists(fileno            HAVE_FILENO)
check_function_exists(_commit           HAVE__COMMIT)
check_function_exists(sigaction         HAVE_SIGACTION)
check_function_exists(sysconf           HAVE_SYSCONF)
check_function_exists(sched_setaffinity HAVE_SCHED_SETAFFINITY)
check_function_exists(sched_getaffinity HAVE_SCHED_GETAFFINITY)
check_function_exists(rsqrt             HAVE_RSQRT)
check_function_exists(rsqrtf            HAVE_RSQRTF)
check_function_exists(sqrtf             HAVE_SQRTF)

include(CheckLibraryExists)
check_library_exists(m sqrt "" HAVE_LIBM)
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)

include(TestBigEndian)
test_big_endian(GMX_INTEGER_BIG_ENDIAN)

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

# Unconditionally find the package, as it is also required for unit tests
include(gmxTestLibXml2)
gmx_test_libxml2(HAVE_LIBXML2)
option(GMX_XML "Use libxml2 to parse xml files (currently, only used in testing)" ${HAVE_LIBXML2})
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
#endif(GMX_XML)

option(GMX_GSL "Add support for gsl" OFF)
if (GMX_GSL)
  find_package(GSL)
  set(PKG_GSL "")
  if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIR})
    set(PKG_GSL gsl)
    set(HAVE_LIBGSL 1)
  endif(GSL_FOUND)
endif (GMX_GSL)

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
endif(GMX_X11)

include(ThreadMPI)
set(THREAD_MPI_LIB thread_mpi)
# Enable core threading facilities
tmpi_enable_core("${CMAKE_SOURCE_DIR}/src/gromacs/legacyheaders")
# Enable tMPI C++ support
tmpi_enable_cxx()
if(GMX_THREAD_MPI)
    # enable MPI functions
    tmpi_enable()
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_THREAD_MPI")
    set(GMX_MPI 1)
    set(MPI_IN_PLACE_EXISTS 1)
endif(GMX_THREAD_MPI)
tmpi_get_source_list(THREAD_MPI_SRC)

if(GMX_GPU)
    # now that we have detected the dependencies, do the second configure pass
    gmx_gpu_setup()
endif(GMX_GPU)

if(CYGWIN)
    set(GMX_CYGWIN 1)
endif(CYGWIN)

if(WIN32 AND NOT CYGWIN)
    set(GMX_NATIVE_WINDOWS 1)
    # This makes windows.h not declare min/max as macros that would break
    # C++ code using std::min/std::max.
    add_definitions(-DNOMINMAX)
endif()

# only bother with finding git and using version.h if the source is a git repo
if(GMX_GIT_VERSION_INFO)
    if (NOT SOURCE_IS_GIT_REPOSITORY)
        message(FATAL_ERROR
            "Cannot generate git version information from source tree not under git. "
            "Set GMX_GIT_VERSION_INFO=OFF to proceed.")
    endif ()
    # We need at least git v1.5.3 be able to parse git's date output. If not
    # found or the version is too small, we can't generate version information.
    find_package(Git)

    if(NOT GIT_FOUND OR GIT_VERSION_STRING VERSION_LESS "1.5.3")
        message(FATAL_ERROR
            "No compatible git version found (>= 1.5.3 required). "
            "Won't be able to generate development version information. "
            "Set GMX_GIT_VERSION_INFO=OFF to proceed.")
    endif()
endif()

# Detect boost unless GMX_EXTERNAL_BOOST is explicitly OFF
# Used for default if GMX_EXTERNAL_BOOST is not defined (first CMake pass)
if(NOT DEFINED GMX_EXTERNAL_BOOST OR GMX_EXTERNAL_BOOST)
    find_package(Boost 1.44.0)
    if(Boost_FOUND AND Boost_VERSION VERSION_LESS "104400")
        set(Boost_FOUND FALSE)
    endif()
    # Print the notification only on first run, when determining the default
    if(NOT DEFINED GMX_EXTERNAL_BOOST AND NOT Boost_FOUND)
        message("Boost >= 1.44 not found. Using minimal internal version. "
                "This may cause trouble if you plan on compiling/linking other "
                "software that uses Boost against Gromacs.")
    endif()
endif()
option(GMX_EXTERNAL_BOOST "Use external Boost instead of minimal built-in version"
       ${Boost_FOUND})
if(GMX_EXTERNAL_BOOST AND NOT Boost_FOUND)
    message(FATAL_ERROR
        "Boost >= 1.44 not found. "
        "You can set GMX_EXTERNAL_BOOST=OFF to compile against minimal "
        "version of Boost included with Gromacs.")
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
# Generate development version info for cache
########################################################################
# set(GEN_VERSION_INFO_INTERNAL "ON")
# include(gmxGenerateVersionString)

########################################################################
# Our own GROMACS tests
########################################################################

add_definitions( -DHAVE_CONFIG_H )
include_directories(${CMAKE_SOURCE_DIR}/src)
# Required for config.h, maybe should only be set in src/CMakeLists.txt
include_directories(${CMAKE_BINARY_DIR}/src)
# Required for gmx_header_config_gen.h to be found before installation
include_directories(${CMAKE_BINARY_DIR}/src/gromacs/utility)
# Required for now to make old code compile
include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/legacyheaders)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

include(gmxSetBuildInformation)
gmx_set_build_information()
if(BUILD_CPU_FEATURES MATCHES "rdtscp" AND NOT GMX_DISTRIBUTABLE_BUILD)
    set(HAVE_RDTSCP 1)
endif(BUILD_CPU_FEATURES MATCHES "rdtscp" AND NOT GMX_DISTRIBUTABLE_BUILD)

include(gmxTestFloatFormat)
gmx_test_float_format(GMX_FLOAT_FORMAT_IEEE754
                      GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER
                      GMX_IEEE754_BIG_ENDIAN_WORD_ORDER)

include(gmxTestLargeFiles)
gmx_test_large_files(GMX_LARGEFILES)

include(gmxTestSignal)
gmx_test_sigusr1(HAVE_SIGUSR1)

include(gmxTestInline)
gmx_test_inline(INLINE_KEYWORD)

include(gmxTestRestrict)
gmx_test_restrict(RESTRICT_KEYWORD)

include(gmxTestPipes)
gmx_test_pipes(HAVE_PIPES)

include(gmxTestIsfinite)
gmx_test_isfinite(HAVE_ISFINITE)
gmx_test__isfinite(HAVE__ISFINITE)
gmx_test__finite(HAVE__FINITE)

include(gmxTestCXX11)
gmx_test_cxx11(GMX_CXX11_SUPPORTED GMX_CXX11_FLAGS)
include(CMakeDependentOption)
# nvcc does not support C++11 flags, so with GPUs we prefer to skip C++11 flags
# entirely to keep the compilation environment uniform.
cmake_dependent_option(GMX_CXX11
    "Use C++11 features"
    ON "GMX_CXX11_SUPPORTED AND NOT GMX_GPU" OFF)
mark_as_advanced(GMX_CXX11)
if(GMX_CXX11)
    set(EXTRA_CXX_FLAGS "${EXTRA_CXX_FLAGS} ${GMX_CXX11_FLAGS}")
endif()


include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
if(NOT GMX_SYSTEM_XDR)
    set(GMX_INTERNAL_XDR 1)
endif(NOT GMX_SYSTEM_XDR)

# include avx test source, used if the AVX flags are set below
include(gmxTestAVXMaskload)

# Process nonbonded accelerated kernels settings
#
# Note that for the backward-compatible x86 SIMD architectures, the
# GMX_CPU_ACCELERATION determines the maximum level of the instruction
# set used (e.g. GMX_CPU_ACCLERATION=SSE4.1 implies
# SSE2). Accordingly, there are a set of CMake variables
# GMX_<arch>_<feature-set> that are exported to the C code to specify
# CPU features that should be used. This means that the logic for
# requiring such backward compatibility is all located here.
if(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    # nothing to do
    set(ACCELERATION_STATUS_MESSAGE "CPU acceleration disabled")

elseif(${GMX_CPU_ACCELERATION} STREQUAL "SSE2")

    GMX_TEST_CFLAG(GNU_SSE2_CFLAG "-msse2" ACCELERATION_C_FLAGS)
    if(NOT GNU_SSE2_CFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CFLAG(MSVC_SSE2_CFLAG "/arch:SSE2" ACCELERATION_C_FLAGS)
    endif(NOT GNU_SSE2_CFLAG AND GMX_NATIVE_WINDOWS)

    GMX_TEST_CXXFLAG(GNU_SSE2_CXXFLAG "-msse2" ACCELERATION_CXX_FLAGS)
    if(NOT GNU_SSE2_CXXFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CXXFLAG(MSVC_SSE2_CXXFLAG "/arch:SSE2" ACCELERATION_CXX_FLAGS)
    endif(NOT GNU_SSE2_CXXFLAG AND GMX_NATIVE_WINDOWS)

    # We dont warn for lacking SSE2 flag support, since that is probably standard today.

    # Only test the include after we have tried to add the correct flag for SSE2 support
    check_include_file(emmintrin.h  HAVE_EMMINTRIN_H ${ACCELERATION_C_FLAGS})

    if(NOT HAVE_EMMINTRIN_H)
        message(FATAL_ERROR "Cannot find emmintrin.h, which is required for SSE2 intrinsics support.")
    endif(NOT HAVE_EMMINTRIN_H)

    set(GMX_CPU_ACCELERATION_X86_SSE2 1)
    # The user should not be able to set this orthogonally to the acceleration
    set(GMX_X86_SSE2 1)
    set(ACCELERATION_STATUS_MESSAGE
        "Enabling SSE2 Gromacs acceleration")

elseif(${GMX_CPU_ACCELERATION} STREQUAL "SSE4.1")

    GMX_TEST_CFLAG(GNU_SSE4_CFLAG "-msse4.1" ACCELERATION_C_FLAGS)
    if (NOT GNU_SSE4_CFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CFLAG(MSVC_SSE4_CFLAG "/arch:SSE4.1" ACCELERATION_C_FLAGS)
    endif(NOT GNU_SSE4_CFLAG AND GMX_NATIVE_WINDOWS)
    if (NOT GNU_SSE4_CFLAG AND NOT MSVC_SSE4_CFLAG)
        # Not surprising if we end up here! MSVC current does not support the SSE4.1 flag. However, it appears to accept SSE4.1
        # intrinsics when SSE2 support is enabled, so we try that instead first.
        if (GMX_NATIVE_WINDOWS)
            GMX_TEST_CFLAG(MSVC_SSE2_CFLAG "/arch:SSE2" ACCELERATION_C_FLAGS)
            message(WARNING "Neither SSE4.1 or SSE2 seems to be supported by your Windows compiler. Something is likely broken.")
        else()
            message(WARNING "No C SSE4.1 flag found. Consider a newer compiler, or use SSE2 for slightly lower performance")
        endif()
    endif(NOT GNU_SSE4_CFLAG AND NOT MSVC_SSE4_CFLAG)

    GMX_TEST_CXXFLAG(GNU_SSE4_CXXFLAG "-msse4.1" ACCELERATION_CXX_FLAGS)
    if (NOT GNU_SSE4_CXXFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CXXFLAG(MSVC_SSE4_CXXFLAG "/arch:SSE4.1" ACCELERATION_CXX_FLAGS)
    endif(NOT GNU_SSE4_CXXFLAG AND GMX_NATIVE_WINDOWS)
    if (NOT GNU_SSE4_CXXFLAG AND NOT MSVC_SSE4_CXXFLAG)
        message(WARNING "No C++ SSE4.1 flag found. Consider a newer compiler, or use SSE2 for slightly lower performance.")
        # Not surprising if we end up here! MSVC current does not support the SSE4.1 flag. However, it appears to accept SSE4.1
        # intrinsics when SSE2 support is enabled, so we try that instead.
        if (GMX_NATIVE_WINDOWS)
            GMX_TEST_CXXFLAG(MSVC_SSE2_CXXFLAG "/arch:SSE2" ACCELERATION_CXX_FLAGS)
        endif()
    endif(NOT GNU_SSE4_CXXFLAG AND NOT MSVC_SSE4_CXXFLAG)

    # This must come after we have added the -msse4.1 flag on some platforms.
    check_include_file(smmintrin.h  HAVE_SMMINTRIN_H ${ACCELERATION_C_FLAGS})

    if(NOT HAVE_SMMINTRIN_H)
        message(FATAL_ERROR "Cannot find smmintrin.h, which is required for SSE4.1 intrinsics support.")
    endif(NOT HAVE_SMMINTRIN_H)

    if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND CMAKE_C_COMPILER_VERSION VERSION_EQUAL "11.1")
        message(FATAL_ERROR "You are using Intel compiler version 11.1, and that compiler is known to produce incorrect results with SSE4.1 acceleration. You need to use another compiler (e.g. icc 12 or newer) or different acceleration (probably slower simulations).")
    endif()

    set(GMX_CPU_ACCELERATION_X86_SSE4_1 1)
    # The user should not be able to set this orthogonally to the acceleration
    set(GMX_X86_SSE4_1 1)
    set(GMX_X86_SSE2   1)
    set(ACCELERATION_STATUS_MESSAGE
        "Enabling SSE4.1 Gromacs acceleration")

elseif(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" OR ${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")

    # Set the AVX compiler flag for both these choices!

    GMX_TEST_CFLAG(GNU_AVX_CFLAG "-mavx" ACCELERATION_C_FLAGS)
    if (NOT GNU_AVX_CFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CFLAG(MSVC_AVX_CFLAG "/arch:AVX" ACCELERATION_C_FLAGS)
    endif (NOT GNU_AVX_CFLAG AND GMX_NATIVE_WINDOWS)
    if (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)
        message(WARNING "No C AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance) giving the -DGMX_CPU_ACCELERATION=SSE4.1 to cmake.")
    endif (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)

    GMX_TEST_CXXFLAG(GNU_AVX_CXXFLAG "-mavx" ACCELERATION_CXX_FLAGS)
    if (NOT GNU_AVX_CXXFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CXXFLAG(MSVC_AVX_CXXFLAG "/arch:AVX" ACCELERATION_CXX_FLAGS)
    endif (NOT GNU_AVX_CXXFLAG AND GMX_NATIVE_WINDOWS)
    if (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)
        message(WARNING "No C++ AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance) giving the -DGMX_CPU_ACCELERATION=SSE4.1 to cmake.")
    endif (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)

    # Set the FMA4 flags (MSVC doesn't require any)
    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" AND NOT MSVC)
        GMX_TEST_CFLAG(GNU_FMA_CFLAG "-mfma4" ACCELERATION_C_FLAGS)
        if (NOT GNU_FMA_CFLAG)
            message(WARNING "No C FMA4 flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
        endif(NOT GNU_FMA_CFLAG)
        GMX_TEST_CFLAG(GNU_XOP_CFLAG "-mxop" ACCELERATION_C_FLAGS)
        # No big deal if we do not have xop, so no point yelling warnings about it.
        GMX_TEST_CXXFLAG(GNU_FMA_CXXFLAG "-mfma4" ACCELERATION_CXX_FLAGS)
        if (NOT GNU_FMA_CXXFLAG)
            message(WARNING "No C++ FMA flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
        endif (NOT GNU_FMA_CXXFLAG)
        GMX_TEST_CXXFLAG(GNU_XOP_CXXFLAG "-mxop" ACCELERATION_CXX_FLAGS)
        # No big deal if we do not have xop, so no point yelling warnings about it.
    endif()

    # Only test the header after we have tried to add the flag for AVX support
    check_include_file(immintrin.h  HAVE_IMMINTRIN_H ${ACCELERATION_C_FLAGS})

    if(NOT HAVE_IMMINTRIN_H)
        message(FATAL_ERROR "Cannot find immintrin.h, which is required for AVX intrinsics support. Consider switching compiler.")
    endif(NOT HAVE_IMMINTRIN_H)

    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")
        try_compile(TEST_AVX ${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/cmake/TestAVX.c"
            COMPILE_DEFINITIONS "${ACCELERATION_C_FLAGS}")
        if(NOT TEST_AVX)
            message(FATAL_ERROR "Cannot compile AVX intrinsics. Consider switching compiler.")
        endif()
    endif()

    # GCC requires x86intrin.h for FMA support. MSVC 2010 requires intrin.h for FMA support.
    check_include_file(x86intrin.h HAVE_X86INTRIN_H ${ACCELERATION_C_FLAGS})
    check_include_file(intrin.h HAVE_INTRIN_H ${ACCELERATION_C_FLAGS})

    # The user should not be able to set this orthogonally to the acceleration
    set(GMX_X86_SSE4_1 1)
    set(GMX_X86_SSE2   1)

    # But just enable one of the choices internally...
    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA")
        # We don't have the full compiler version string yet (BUILD_C_COMPILER),
        # so we can't distinguish vanilla and Apple clang, but catering for AMD
        # hackintoshes is not worth the effort.
        if (APPLE AND (${CMAKE_C_COMPILER_ID} STREQUAL "Clang" OR
                    ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
            message(WARNING "Due to a known compiler bug, Clang up to version 3.2 (and Apple Clang up to version 4.1) produces incorrect code with AVX_128_FMA acceleration. As we can not work around this bug on OS X, you will have to select a different compiler or CPU acceleration.")
        endif()

        if (GMX_USE_CLANG_C_FMA_BUG_WORKAROUND)
            # we assume that we have an external assembler that supports AVX
            message(STATUS "Clang ${CMAKE_C_COMPILER_VERSION} detected, enabling FMA bug workaround")
            set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -no-integrated-as")
        endif()
        if (GMX_USE_CLANG_CXX_FMA_BUG_WORKAROUND)
            # we assume that we have an external assembler that supports AVX
            message(STATUS "Clang ${CMAKE_CXX_COMPILER_VERSION} detected, enabling FMA bug workaround")
            set(EXTRA_CXX_FLAGS "${EXTRA_CXX_FLAGS} -no-integrated-as")
        endif()

        set(GMX_CPU_ACCELERATION_X86_AVX_128_FMA 1)
        set(GMX_X86_AVX_128_FMA 1)
        set(ACCELERATION_STATUS_MESSAGE
            "Enabling 128-bit AVX Gromacs acceleration (with fused-multiply add)")

    else()
        # If we are not doing AVX_128, it must be AVX_256...
        set(GMX_CPU_ACCELERATION_X86_AVX_256 1)
        set(GMX_X86_AVX_256 1)
        set(ACCELERATION_STATUS_MESSAGE
            "Enabling 256-bit AVX Gromacs acceleration")
    endif()

    # Unfortunately gcc-4.5.2 and gcc-4.6.0 has a bug where they use the wrong datatype for the formal
    # parameter of the mask for maskload/maskstore arguments. Check if this is present, since we can work around it.
    gmx_test_avx_gcc_maskload_bug(${ACCELERATION_C_FLAGS} GMX_X86_AVX_GCC_MASKLOAD_BUG)

elseif(${GMX_CPU_ACCELERATION} STREQUAL "IBM_QPX")
    try_compile(TEST_QPX ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestQPX.c")

    if (TEST_QPX)
        message(WARNING "IBM QPX acceleration was selected. This will work, but SIMD-accelerated kernels are only available for the Verlet cut-off scheme. The plain C kernels that are used for the group cut-off scheme kernels will be slow, so please consider using the Verlet cut-off scheme.")
        set(GMX_CPU_ACCELERATION_IBM_QPX 1)
    else()
        message(FATAL_ERROR "Cannot compile the requested IBM QPX intrinsics. If you are compiling for BlueGene/Q with the XL compilers, use 'cmake .. -DCMAKE_TOOLCHAIN_FILE=Platform/BlueGeneQ-static-XL-C' to set up the tool chain.")
    endif()
elseif(${GMX_CPU_ACCELERATION} STREQUAL "SPARC64_HPC_ACE")
    set(GMX_CPU_ACCELERATION_SPARC64_HPC_ACE 1)
else()
    gmx_invalid_option_value(GMX_CPU_ACCELERATION)
endif()
gmx_check_if_changed(ACCELERATION_CHANGED GMX_CPU_ACCELERATION)
if (ACCELERATION_CHANGED AND DEFINED ACCELERATION_STATUS_MESSAGE)
    message(STATUS "${ACCELERATION_STATUS_MESSAGE}")
endif()

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
  if(CYGWIN OR NOT WIN32)
    # Native Windows does not have, nor need dlopen
    # Note that WIN32 is set with Cygwin, but Cygwin needs dlopen to use plug-ins
    include(gmxTestdlopen)
    gmx_test_dlopen(HAVE_DLOPEN)
  endif()

  # so, should we use plug-ins?
  if((WIN32 AND NOT CYGWIN) OR (HAVE_DLOPEN AND BUILD_SHARED_LIBS))
    if(NOT VMD_QUIETLY)
      MESSAGE(STATUS
          "Found the ability to use plug-ins when building shared libaries, "
          "so will compile to use plug-ins (e.g. to read VMD-supported file "
          "formats).")
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
endif(GMX_LOAD_PLUGINS)
set(VMD_QUIETLY TRUE CACHE INTERNAL "")

# Link real-time library for POSIX timers. The check for clock_gettime
# confirms the linkability of rt.
if(HAVE_TIME_H AND HAVE_UNISTD_H AND HAVE_CLOCK_GETTIME)
    list(APPEND GMX_EXTRA_LIBRARIES rt)
endif()

# Math and thread libraries must often come after all others when linking...
if(HAVE_LIBM)
    list(APPEND GMX_EXTRA_LIBRARIES m)
endif(HAVE_LIBM)
if (${CMAKE_SYSTEM_NAME} MATCHES "BlueGene")
    check_library_exists(mass_simd atan2f4 "" HAVE_MASS_SIMD)
    if(HAVE_MASS_SIMD)
        list(APPEND GMX_EXTRA_LIBRARIES mass_simd)
    else()
        message(FATAL_ERROR "Could not link to the SIMD version of the IBM MASS library. Please adjust your CMAKE_PREFIX_PATH to contain it")
    endif()
endif()


option(GMX_NACL "Configure for Native Client builds" OFF)
if (GMX_NACL)
  list(APPEND GMX_EXTRA_LIBRARIES nosys)
  set(GMX_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lnosys")
  set(GMX_NO_NICE 1)
  set(GMX_NO_RENAME 1)
endif()
mark_as_advanced(GMX_NACL)

if(GMX_FAHCORE)
  set(COREWRAP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../corewrap" CACHE STRING
      "Path to swindirect.h")
  include_directories(${COREWRAP_INCLUDE_DIR})
endif(GMX_FAHCORE)

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (NOT GMX_SKIP_DEFAULT_CFLAGS)
    set(CMAKE_C_FLAGS "${ACCELERATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${ACCELERATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${FFT_LINKER_FLAGS} ${MPI_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${MPI_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
else()
    message("Recommended flags which are not added because GMX_SKIP_DEFAULT_CFLAGS=yes:")
    message("CMAKE_C_FLAGS: ${ACCELERATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${GMXC_CFLAGS}")
    message("CMAKE_C_FLAGS_RELEASE: ${GMXC_CFLAGS_RELEASE}")
    message("CMAKE_C_FLAGS_DEBUG: ${GMXC_CFLAGS_DEBUG}")
    message("CMAKE_CXX_FLAGS: ${ACCELERATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_CXX_FLAGS} ${GMXC_CXXFLAGS}")
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
set(PKG_CFLAGS "${PKG_CFLAGS} ${OpenMP_C_FLAGS}")

########################################################################
# Specify install locations
########################################################################
set(GMX_LIB_INSTALL_DIR lib CACHE STRING
    "Library installation directory (default: lib)")
set(GMX_DATA_INSTALL_DIR gromacs CACHE STRING
    "Data installation directory under share/ (default: gromacs)")
mark_as_advanced(GMX_LIB_INSTALL_DIR GMX_DATA_INSTALL_DIR)

set(LIB_INSTALL_DIR  ${GMX_INSTALL_PREFIX}${GMX_LIB_INSTALL_DIR})
set(BIN_INSTALL_DIR  ${GMX_INSTALL_PREFIX}bin)
set(DATA_INSTALL_DIR ${GMX_INSTALL_PREFIX}share/${GMX_DATA_INSTALL_DIR})
set(MAN_INSTALL_DIR  ${GMX_INSTALL_PREFIX}share/man)
set(INCL_INSTALL_DIR ${GMX_INSTALL_PREFIX}include)

set(GMXLIB_SEARCH_DIR share/${GMX_DATA_INSTALL_DIR}/top)
set(GMXLIB_FALLBACK   ${CMAKE_INSTALL_PREFIX}/${DATA_INSTALL_DIR}/top)

# Binary and library suffix options
include(gmxManageSuffixes)

##################################################################
# Shared library settings
##################################################################
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "\\\$ORIGIN/../${GMXLIB}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    if(CMAKE_SYSTEM_VERSION VERSION_GREATER 8.0) #rpath supported for >10.4
        set(CMAKE_INSTALL_NAME_DIR "@rpath")
        set(GMX_EXE_LINKER_FLAGS ${GMX_EXE_LINKER_FLAGS} "-Wl,-rpath,@executable_path/../lib")
    else()
        set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
    endif()
endif()

#COPYING file: Only necessary for binary distributions.
#Simpler to always install.
install(FILES COPYING DESTINATION ${DATA_INSTALL_DIR} COMPONENT data)

if(GMX_EXTERNAL_BOOST)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    include_directories(${CMAKE_SOURCE_DIR}/src/external/boost)
    # typeid not supported for minimal internal version
    # (would add significant amount of code)
    add_definitions(-DBOOST_NO_TYPEID)
    # TODO: Propagate the above settings to the installed CMakeFiles.txt template
    # (from share/template/)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DBOOST_NO_TYPEID -I${CMAKE_INSTALL_PREFIX}/${INCL_INSTALL_DIR}/gromacs/external/boost")
    if (NOT GMX_BUILD_MDRUN_ONLY)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/external/boost/boost
                DESTINATION ${INCL_INSTALL_DIR}/gromacs/external/boost
                COMPONENT development)
    endif ()
endif()

# TODO The Windows builds are still broken (we are avoiding some
# unknown bug, perhaps with checkpointing?), because defaulting to
# .tng output leads to calling null-op implementation functions, so
# regressiontests will fail on Windows at the post-mdrun stages. We
# can fix this during the 5.0 beta phase.
gmx_add_cache_dependency(GMX_USE_TNG BOOL "NOT GMX_NATIVE_WINDOWS" OFF)

if(GMX_USE_TNG)
    find_package(ZLIB)
    set(TNG_BUILD_FORTRAN OFF CACHE BOOL "Build Fortran compatible TNG library and examples for testing")
    set(TNG_BUILD_EXAMPLES OFF CACHE BOOL "Build examples showing usage of the TNG API")
    set(TNG_BUILD_COMPRESSION_TESTS OFF CACHE BOOL "Build tests of the TNG compression library")
    # TODO reconsider use of OpenMP in TNG if/when trajectory writing is not serial in mdrun
    set(TNG_USE_OPENMP ${GMX_OPENMP} CACHE BOOL "Try to use the OpenMP library (if available) with TNG")
    set(TNG_BUILD_DOCUMENTATION OFF CACHE BOOL "Use Doxygen to create the HTML based TNG API documentation")
    set(TNG_BUILD_TEST OFF CACHE BOOL "Build TNG testing binary.")
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/external/tng_io)
    set(GMX_TNG_LIBRARIES tng_io)
endif()
mark_as_advanced(TNG_BUILD_FORTRAN)
mark_as_advanced(TNG_BUILD_EXAMPLES)
mark_as_advanced(TNG_BUILD_COMPRESSION_TESTS)
mark_as_advanced(TNG_USE_OPENMP)
mark_as_advanced(TNG_BUILD_DOCUMENTATION)
mark_as_advanced(TNG_BUILD_TEST)
mark_as_advanced(TNG_EXAMPLE_FILES_DIR)

if (GMX_BUILD_FOR_COVERAGE)
    # Code heavy with asserts makes conditional coverage close to useless metric,
    # as by design most of the false branches are impossible to trigger in
    # correctly functioning code.  And the benefit of testing those that could
    # be triggered by using an API against its specification isn't usually
    # worth the effort.
    add_definitions(-DNDEBUG -DBOOST_DISABLE_ASSERTS -DGMX_DISABLE_ASSERTS)
endif ()

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
    add_subdirectory(doxygen)
    add_subdirectory(share)
    add_subdirectory(scripts)
endif ()
add_subdirectory(src)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Issue a warning if NVIDIA GPUs were detected, but CUDA was not found.
# Don't bother the user after the first configure pass.
if ((CUDA_NOTFOUND_AUTO AND GMX_DETECT_GPU_AVAILABLE) AND NOT GMX_GPU_DETECTION_DONE)
    message(WARNING "${CUDA_NOTFOUND_MESSAGE}")
endif()
set(GMX_GPU_DETECTION_DONE TRUE CACHE INTERNAL "Whether GPU detection has already been done")

#######################
## uninstall target
#######################
    CONFIGURE_FILE(
                   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
                   IMMEDIATE @ONLY)
###########################
ADD_CUSTOM_TARGET(uninstall
                  "${CMAKE_COMMAND}" -P
                  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")
###########################

########################################################################
# Manual                                                               #
########################################################################

option(GMX_BUILD_MANUAL "Whether to try to configure to build the PDF manual" OFF)
mark_as_advanced(GMX_BUILD_MANUAL)
if(GMX_BUILD_MANUAL)
    # Make sure we only do detection of manual-building dependencies
    # when the user opted in for that.
    add_subdirectory(manual)
endif()
