#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#
cmake_minimum_required(VERSION 2.8)
# Keep CMake suitably quiet on Cygwin
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

# Allows CPack to act differently for normal tools and mdrun (e.g. because of MPI)
set(CPACK_COMPONENT_GROUP_TOOLS_DESCRIPTION "All GROMACS executable tools")
set(CPACK_COMPONENT_GROUP_MDRUN_DESCRIPTION "GROMACS executable for running simulations")

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform)

project(Gromacs C)
include(Dart)
mark_as_advanced(DART_ROOT)

# PROJECT_VERSION should have the following structure: 
# VERSION-dev[-SUFFIX] where the VERSION should have the for: vMajor.vMinor.vPatch
#
# The "-dev" suffix is important to keep because it makes possible to distinguish 
# between a build from official release and a build from git release branch on a 
# machine with no git. 
#
# NOTE: when releasing the "-dev" suffix needs to be stripped off!
set(PROJECT_VERSION "4.6.4-dev")
# The version number of the regressiontest tarball against which this
# git branch can be tested. Normally, this will be the version of the
# last patch release. Comment the next line out for branches leading
# to a major/minor release.
set(REGRESSIONTEST_VERSION "4.6.2")
set(CUSTOM_VERSION_STRING ""
    CACHE STRING "Custom version string (if empty, use hard-coded default)")
mark_as_advanced(CUSTOM_VERSION_STRING)
if (CUSTOM_VERSION_STRING)
    set(PROJECT_VERSION ${CUSTOM_VERSION_STRING})
endif (CUSTOM_VERSION_STRING)
set(SOVERSION 8)
# It is a bit irritating, but this has to be set separately for now!
SET(CPACK_PACKAGE_VERSION_MAJOR "4")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
SET(CPACK_PACKAGE_VERSION_PATCH "1")

# The numerical gromacs version. It is 40600 for 4.6.0.
# The #define GMX_VERSION in gmx_header_config_h is set to this value.
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

include(gmxBuildTypeReference)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Reference." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(C)

set(GMX_USE_RELATIVE_INSTALL_PATH OFF CACHE STRING "Use relative paths not absolute paths for cmake install. Has only an effect on cpack.")
mark_as_advanced(GMX_USE_RELATIVE_INSTALL_PATH)

set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_VENDOR "gromacs.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gromacs - a toolkit for high-performance molecular simulation")
if (NOT GMX_USE_RELATIVE_INSTALL_PATH)
    set(CPACK_SET_DESTDIR "ON")
endif()
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/admin/InstallWelcome.txt")
# Its GPL/LGPL, so they do not have to agree to a license for mere usage, but some installers require this...
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/admin/InstallInfo.txt")
set(CPACK_SOURCE_IGNORE_FILES "\\\\.isreposource$;\\\\.git/;\\\\.gitignore$")
set(CPACK_PROJECT_CONFIG_FILE "${CMAKE_SOURCE_DIR}/CPackInit.cmake")
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_SOURCE_DIR};/;${CMAKE_BINARY_DIR}/man;man")
set(CPACK_PACKAGE_CONTACT "gmx-users@gromacs.org")

#must come after all cpack settings!
include(CPack)

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
# User input options - enable C++ - before any CXX flags are changed   #
########################################################################

# decide on GPU settings based on user-settings and GPU/CUDA detection
include(gmxManageGPU)

option(GMX_FORCE_CXX "Enable C++ compilation even if not necessary" OFF)
mark_as_advanced(GMX_FORCE_CXX)

option(GMX_COOL_QUOTES "Enable Gromacs cool quotes" ON)
mark_as_advanced(GMX_COOL_QUOTES)

if(GMX_GPU OR GMX_FORCE_CXX OR GMX_OPENMM)
    enable_language(CXX)
endif()
set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")

########################################################################
# User input options                                                   #
########################################################################
option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" OFF)
option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREAD_MPI  "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" ON)
option(GMX_SOFTWARE_INVSQRT "Use GROMACS software 1/sqrt" ON)
mark_as_advanced(GMX_SOFTWARE_INVSQRT)
option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
mark_as_advanced(GMX_FAHCORE)

include(gmxDetectAcceleration)
if(NOT DEFINED GMX_CPU_ACCELERATION)
    if(CMAKE_CROSSCOMPILING)
        if("${CMAKE_SYSTEM_NAME}" MATCHES "BlueGeneQ")
            set(GMX_SUGGESTED_CPU_ACCELERATION "IBM_QPX")
        else()
            set(GMX_SUGGESTED_CPU_ACCELERATION "None")
        endif()
    else(CMAKE_CROSSCOMPILING)
        gmx_detect_acceleration(GMX_SUGGESTED_CPU_ACCELERATION)
    endif(CMAKE_CROSSCOMPILING)
endif(NOT DEFINED GMX_CPU_ACCELERATION)

set(GMX_CPU_ACCELERATION "@GMX_SUGGESTED_CPU_ACCELERATION@"
    CACHE STRING "Accelerated CPU kernels. Pick one of: None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, IBM_QPX, Sparc64_HPC_ACE")

set(GMX_FFT_LIBRARY "fftw3" 
    CACHE STRING "FFT library choices: fftw3,mkl,fftpack[built-in]")
option(GMX_BUILD_OWN_FFTW "Download and build FFTW 3 during the GROMACS build process, rather than fall back on the really slow fftpack." OFF)
mark_as_advanced(GMX_BUILD_OWN_FFTW)
option(GMX_DISABLE_FFTW_MEASURE 
       "Do not optimize FFTW setups (not needed with SSE)" OFF)
mark_as_advanced(GMX_DISABLE_FFTW_MEASURE)
set(GMX_QMMM_PROGRAM "none" 
    CACHE STRING "QM package choices: none,gaussian,mopac,gamess,orca")
option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
mark_as_advanced(GMX_BROKEN_CALLOC)
option(GMX_MPI_IN_PLACE "Enable MPI_IN_PLACE for MPIs that have it defined" ON)
mark_as_advanced(GMX_MPI_IN_PLACE)
option(GMX_LOAD_PLUGINS "Compile with plugin support, needed to read VMD supported file formats" ON)
mark_as_advanced(GMX_LOAD_PLUGINS)

option(GMX_OPENMP "Enable OpenMP-based multithreading" ON)

option(USE_VERSION_H "Generate development version string/information" ON)
mark_as_advanced(USE_VERSION_H)

option(GMX_DEFAULT_SUFFIX "Use default suffixes for GROMACS binaries and libs (_d for double, _mpi for MPI; rerun cmake after changing to see relevant options)" ON)

if(UNIX)
    option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer static archives (not available on non-*nix platforms and it will only work if static versions of external dependencies are available and found)!" OFF)
    mark_as_advanced(GMX_PREFER_STATIC_LIBS)
endif()

option(GMX_CYCLE_SUBCOUNTERS "Enable cycle subcounters to get a more detailed cycle timings" OFF)
mark_as_advanced(GMX_CYCLE_SUBCOUNTERS)

option(GMX_SKIP_DEFAULT_CFLAGS "Don't automatically add suggested/required Compiler flags." OFF)
mark_as_advanced(GMX_SKIP_DEFAULT_CFLAGS)

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
if(CMAKE_CXX_COMPILER_LOADED)
    include(CheckCXXCompilerFlag)
endif()

# Get compiler version information, needs to be done early as check that depend
# on compiler verison follow below.
include(gmxGetCompilerInfo)
get_compiler_version()

# First exclude compilers known to not work with OpenMP although claim to support it:
# gcc 4.2.1 and gcc-llvm 4.2.1 (also claims to be 4.2.1) on Mac OS X
# This fixes redmine 900 and needs to run before OpenMP flags are set below.
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND
    (CMAKE_COMPILER_IS_GNUCC AND C_COMPILER_VERSION AND C_COMPILER_VERSION VERSION_LESS 4.3))
    message(STATUS "OpenMP multithreading not supported with gcc/llvm-gcc 4.2 on Mac OS X, disabled")
    set(GMX_OPENMP OFF CACHE BOOL
        "OpenMP multithreading not not supported with gcc/llvm-gcc 4.2 on Mac OS X, disabled!" FORCE)
endif()

# OpenMP check must come before other CFLAGS!
if(GMX_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        # CMake on Windows doesn't support linker flags passed to target_link_libraries
        # (i.e. it treats /openmp as \openmp library file). Also, no OpenMP linker flags are needed.
        if(NOT (WIN32 AND NOT CYGWIN))
            if(CMAKE_COMPILER_IS_GNUCC AND GMX_PREFER_STATIC_OPENMP)
                set(OpenMP_LINKER_FLAGS "-Wl,-static -lgomp -lrt -Wl,-Bdynamic -lpthread")
                set(OpenMP_SHARED_LINKER_FLAGS "")
            else()
                # Only set a linker flag if the user didn't set them manually
                if(NOT DEFINED OpenMP_LINKER_FLAGS)
                    set(OpenMP_LINKER_FLAGS "${OpenMP_C_FLAGS}")
                endif()
                if(NOT DEFINED OpenMP_SHARED_LINKER_FLAGS)
                    set(OpenMP_SHARED_LINKER_FLAGS "${OpenMP_C_FLAGS}")
                endif()
            endif()
        endif()
    else(OPENMP_FOUND)
        message(WARNING
                "The compiler you are using does not support OpenMP parallelism. This might hurt your performance a lot, in particular with GPUs. Try using a more recent version, or a different compiler. For now, we are proceeding by turning off OpenMP.")
        set(GMX_OPENMP OFF CACHE STRING "Whether GROMACS will use OpenMP parallelism." FORCE)
    endif(OPENMP_FOUND)
endif()


include(gmxCFlags)
gmx_c_flags()

# This variable should be used for additional compiler flags which are not
# generated in gmxCFlags nor are acceleration or MPI related.
set(EXTRA_C_FLAGS "")

# gcc 4.4.x is buggy and crashes when compiling some files with O3 and OpenMP on.
# Detect here whether applying a workaround is needed and will apply it later
# on the affected files.
include(gmxGCC44O3BugWorkaround)
gmx_check_gcc44_bug_workaround_needed(GMX_USE_GCC44_BUG_WORKAROUND)

# clang 3.0 is buggy for some unknown reason detected during adding
# the SSE2 group kernels for GROMACS 4.6. If we ever work out what
# that is, we should replace these tests with a compiler feature test,
# update GROMACS Redmine task #1039 and perhaps report a clang bug.
#
# In the meantime, until we require CMake 2.8.10 we cannot rely on it to detect
# the compiler version for us. So we need a manual check for clang 3.0.
include(gmxDetectClang30)
gmx_detect_clang_3_0(COMPILER_IS_CLANG_3_0)
if(COMPILER_IS_CLANG_3_0)
    message(FATAL_ERROR "Your compiler is clang version 3.0, which is known to be buggy for GROMACS. Use a different compiler.")
endif()

# clang <=3.2 contains a bug that causes incorrect code to be generated for the
# vfmaddps instruction and therefore the bug is triggered with AVX_128_FMA.
# (see: http://llvm.org/bugs/show_bug.cgi?id=15040).
# We can work around this by not using the integrated assembler (except on OS X
# which has an outdated assembler that does not support AVX instructions).
if (${CMAKE_C_COMPILER_ID} MATCHES "Clang" AND C_COMPILER_VERSION VERSION_LESS "3.3")
    set(GMX_USE_CLANG_FMA_BUG_WORKAROUND TRUE)
endif()

if (CMAKE_C_COMPILER_ID STREQUAL "PGI")
    message(WARNING "All tested PGI compiler versions (up to 12.9.0) generate binaries which produce incorrect results, or even fail to compile Gromacs. Highly recommended to use a different compiler. If you choose to use PGI, make sure to run the regressiontests.")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND C_COMPILER_VERSION VERSION_LESS "12.0.0")
    message(WARNING "Intel compilers before 12.0.0 are not routinely tested, so there may be problems. Version 11.1 with SSE4.1 is known to produce incorrect results. It is highly recommended to use a more up-to-date compiler. If you choose to use this version, make sure you run the regressiontests.")
endif()

set(PKG_CFLAGS "")
if(GMX_DOUBLE)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_DOUBLE")
endif(GMX_DOUBLE)
if(GMX_SOFTWARE_INVSQRT)
  set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_SOFTWARE_INVSQRT")
endif(GMX_SOFTWARE_INVSQRT)



########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
check_include_files(string.h     HAVE_STRING_H)
check_include_files(math.h       HAVE_MATH_H)
check_include_files(limits.h     HAVE_LIMITS_H)
check_include_files(memory.h     HAVE_MEMORY_H)
check_include_files(unistd.h	 HAVE_UNISTD_H)
check_include_files(direct.h	 HAVE_DIRECT_H)
check_include_files(pwd.h        HAVE_PWD_H)
check_include_files(stdint.h	 HAVE_STDINT_H)
check_include_files(stdlib.h	 HAVE_STDLIB_H)
check_include_files(pthread.h    HAVE_PTHREAD_H)
check_include_files(dirent.h     HAVE_DIRENT_H)
check_include_files(inttypes.h   HAVE_INTTYPES_H)
check_include_files(regex.h      HAVE_REGEX_H)
check_include_files(sys/types.h  HAVE_SYS_TYPES_H)
check_include_files(sys/stat.h   HAVE_SYS_STAT_H)
check_include_files(sys/time.h   HAVE_SYS_TIME_H)
check_include_files(rpc/rpc.h    HAVE_RPC_RPC_H)
check_include_files("rpc/rpc.h;rpc/xdr.h"    HAVE_RPC_XDR_H)
check_include_files(io.h  		 HAVE_IO_H)
check_include_files(sched.h      HAVE_SCHED_H)

include(CheckFunctionExists)
check_function_exists(strcasecmp        HAVE_STRCASECMP)
check_function_exists(strdup            HAVE_STRDUP)
check_function_exists(vprintf           HAVE_VPRINTF)
check_function_exists(memcmp            HAVE_MEMCMP)
check_function_exists(posix_memalign    HAVE_POSIX_MEMALIGN)
check_function_exists(memalign          HAVE_MEMALIGN)
check_function_exists(_aligned_malloc   HAVE__ALIGNED_MALLOC)
check_function_exists(gettimeofday      HAVE_GETTIMEOFDAY)
check_function_exists(isnan             HAVE_ISNAN)
check_function_exists(_isnan            HAVE__ISNAN)
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
check_library_exists(m cbrt "" HAVE_CBRT)


include(CheckTypeSize)

check_type_size("bool"          SIZEOF_BOOL) # will also set HAVE_BOOL
check_type_size("int"           SIZEOF_INT) 
check_type_size("long int"      SIZEOF_LONG_INT) 
check_type_size("long long int" SIZEOF_LONG_LONG_INT) 
check_type_size("off_t"         SIZEOF_OFF_T)
check_type_size("void *"        SIZEOF_VOIDP)

if (CMAKE_C_SIZEOF_DATA_PTR EQUAL 8)
    set(GMX_64_BIT TRUE)
else (CMAKE_C_SIZEOF_DATA_PTR EQUAL 8)
    set(GMX_64_BIT FALSE)
endif (CMAKE_C_SIZEOF_DATA_PTR EQUAL 8)

# Check for some basic types that we *need*, so set these to int if they are not present 
check_type_size(uid_t uid_t)
if(NOT uid_t)
  set(uid_t int)
else(NOT uid_t)
  set(uid_t 0)
endif(NOT uid_t)

check_type_size(gid_t gid_t)
if(NOT gid_t)
  set(gid_t 1)
else(NOT gid_t)
  set(gid_t 0)
endif(NOT gid_t)

check_type_size(size_t size_t)
if(NOT size_t)
  set(size_t int)
else(NOT size_t)
  set(size_t 0)
endif(NOT size_t)

check_type_size(off_t off_t)
if(NOT off_t)
  set(off_t int)
else(NOT off_t)
  set(off_t 0)
endif(NOT off_t)

include(TestBigEndian)
test_big_endian(GMX_INTEGER_BIG_ENDIAN)


if(APPLE OR CYGWIN OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux|.*BSD")
    # Maybe Solaris should be here? Patch this if you know!
    SET(SHARED_LIBS_DEFAULT ON)
elseif(WIN32 OR ${CMAKE_SYSTEM_NAME} MATCHES "BlueGene")
    # Support for shared libs on native Windows is a bit new. Its
    # default might change later if/when we sort things out. Also,
    # Cray should go here. What variable value can detect it?
    SET(SHARED_LIBS_DEFAULT OFF)
else()
    message(STATUS "Defaulting to building static libraries")
    SET(SHARED_LIBS_DEFAULT OFF)
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

if(UNIX AND GMX_PREFER_STATIC_LIBS AND SHARED_LIBS_DEFAULT)
    if(BUILD_SHARED_LIBS)
        # Warn the user about the combination. But don't overwrite the request.
        message(WARNING "Searching for static libraries requested, and building shared Gromacs libraries requested. This might cause problems linking later.")
    elseif(NOT DEFINED BUILD_SHARED_LIBS)
        # Change default to OFF. Don't warn if it's already off.
        message(WARNING "Searching for static libraries requested, so the GROMACS libraries will also be built statically (BUILD_SHARED_LIBS=OFF)")
        set(SHARED_LIBS_DEFAULT OFF)
    endif()
endif()

# By now, all tool chains should have spoken up if they care about
# the setting of SHARED_LIBS_DEFAULT.
option(BUILD_SHARED_LIBS "Enable shared libraries (can be problematic e.g. with MPI, or on some HPC systems)" ${SHARED_LIBS_DEFAULT})

########################################################################
#Process MPI settings
########################################################################
include(gmxManageMPI)

########################################################################
# Find external packages                                               #
########################################################################
if(UNIX AND GMX_PREFER_STATIC_LIBS)
    # On Linux .a is the static library suffix, on Mac OS X .lib can also
    # be used, so we'll add both to the preference list.
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

IF( WIN32 AND NOT CYGWIN)
  if (NOT BUILD_SHARED_LIBS)
      option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer static system libraries (MT instead of MD)!" ON)
      if(NOT GMX_PREFER_STATIC_LIBS)
          message(WARNING "Shared system libraries requested, and static Gromacs libraries requested.")
      endif()
  else()
      option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer static system libraries (MT instead of MD)!" OFF)
      if(GMX_PREFER_STATIC_LIBS)
          #this combination segfaults (illigal passing of file handles)
          message(FATAL_ERROR "Static system libraries requested, and shared Gromacs libraries requested.")
      endif()
      add_definitions(-DUSE_VISIBILITY -DTMPI_USE_VISIBILITY)
      set(PKG_CFLAGS "$PKG_CFLAGS -DUSE_VISIBILITY -DTMPI_USE_VISIBILITY")
  endif()
  mark_as_advanced(GMX_PREFER_STATIC_LIBS)

  IF (GMX_PREFER_STATIC_LIBS)
      #Only setting Debug and Release flags. Others configurations are current not used.
      STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
      STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
      if(CMAKE_CXX_COMPILER_LOADED)
          STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
          STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
      endif()
  ENDIF()
  IF( CMAKE_C_COMPILER_ID MATCHES "Intel" )
    if(BUILD_SHARED_LIBS) #not sure why incremental building with shared libs doesn't work
        STRING(REPLACE "/INCREMENTAL:YES" "" CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS})
    endif()
  ENDIF()
ENDIF()

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
        message(WARNING "X11 include files and/or libraries were not found. Will not build the GROMACS X11-binaries, such as ngmx")
    endif()
endif(GMX_X11)

include(ThreadMPI)
set(THREAD_MPI_LIB thread_mpi)
# Enable core threading facilities
tmpi_enable_core("${CMAKE_SOURCE_DIR}/include")
# Enable tMPI C++ support
# tmpi_enable_cxx()
if(GMX_THREAD_MPI)
    # enable MPI functions
    tmpi_enable()
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_THREAD_MPI")
    set(GMX_MPI 1)
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
endif()

# only bother with finding git and using version.h if the source is a git repo
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    if(USE_VERSION_H)
        # We need at least git v1.5.3 be able to parse git's date output. If not
        # found or the version is too small, we can't generate version information.
        find_package(Git)

	# Find out the git version
	if(GIT_FOUND AND NOT GIT_VERSION)
	  execute_process(COMMAND ${GIT_EXECUTABLE} "--version"
            OUTPUT_VARIABLE _exec_out
            OUTPUT_STRIP_TRAILING_WHITESPACE)
	  string(REGEX REPLACE "git version (.*)" "\\1" GIT_VERSION ${_exec_out})
	  set(GIT_VERSION ${GIT_VERSION} CACHE STRING "Git version")
	  mark_as_advanced(GIT_VERSION)
	endif()

        if(NOT GIT_FOUND OR GIT_VERSION VERSION_LESS "1.5.3")
          message("No compatible git version found, won't be able to generate proper development version information.")
          set(USE_VERSION_H OFF)
        endif()
    endif()
else()
    set(USE_VERSION_H OFF)
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
include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

include(gmxSetBuildInformation)
gmx_set_build_information()
if(BUILD_CPU_FEATURES MATCHES "rdtscp" AND NOT GMX_DISTRIBUTABLE_BUILD)
    # The timestep counter headers do not include config.h
    add_definitions(-DHAVE_RDTSCP)
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

include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
if(NOT GMX_SYSTEM_XDR)
    set(GMX_INTERNAL_XDR 1)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_INTERNAL_XDR")
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
string(TOUPPER ${GMX_CPU_ACCELERATION} GMX_CPU_ACCELERATION)
if(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    # nothing to do
elseif(${GMX_CPU_ACCELERATION} STREQUAL "SSE2")

    GMX_TEST_CFLAG(GNU_SSE2_CFLAG "-msse2" ACCELERATION_C_FLAGS)
    if(NOT GNU_SSE2_CFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CFLAG(MSVC_SSE2_CFLAG "/arch:SSE2" ACCELERATION_C_FLAGS)
    endif(NOT GNU_SSE2_CFLAG AND GMX_NATIVE_WINDOWS)

    if (CMAKE_CXX_COMPILER_LOADED)
        GMX_TEST_CXXFLAG(GNU_SSE2_CXXFLAG "-msse2" ACCELERATION_CXX_FLAGS)
        if(NOT GNU_SSE2_CXXFLAG AND GMX_NATIVE_WINDOWS)
            GMX_TEST_CXXFLAG(MSVC_SSE2_CXXFLAG "/arch:SSE2" ACCELERATION_CXX_FLAGS)
        endif(NOT GNU_SSE2_CXXFLAG AND GMX_NATIVE_WINDOWS)
    endif()

    # We dont warn for lacking SSE2 flag support, since that is probably standard today.

    # Only test the include after we have tried to add the correct flag for SSE2 support
    check_include_file(emmintrin.h  HAVE_EMMINTRIN_H ${ACCELERATION_C_FLAGS})

    if(NOT HAVE_EMMINTRIN_H)
        message(FATAL_ERROR "Cannot find emmintrin.h, which is required for SSE2 intrinsics support.")
    endif(NOT HAVE_EMMINTRIN_H)

    set(GMX_CPU_ACCELERATION_X86_SSE2 1)
    # The user should not be able to set this orthogonally to the acceleration
    set(GMX_X86_SSE2 1)
    if (NOT ACCELERATION_QUIETLY)
      message(STATUS "Enabling SSE2 Gromacs acceleration, and it will help compiler optimization.")
    endif()

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

    if (CMAKE_CXX_COMPILER_LOADED)
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
    endif()

    # This must come after we have added the -msse4.1 flag on some platforms.
    check_include_file(smmintrin.h  HAVE_SMMINTRIN_H ${ACCELERATION_C_FLAGS})

    if(NOT HAVE_SMMINTRIN_H)
        message(FATAL_ERROR "Cannot find smmintrin.h, which is required for SSE4.1 intrinsics support.")
    endif(NOT HAVE_SMMINTRIN_H)

    set(GMX_CPU_ACCELERATION_X86_SSE4_1 1)
    # The user should not be able to set this orthogonally to the acceleration
    set(GMX_X86_SSE4_1 1)
    set(GMX_X86_SSE2   1)
    if (NOT ACCELERATION_QUIETLY)
      message(STATUS "Enabling SSE4.1 Gromacs acceleration, and it will help compiler optimization.")
    endif()

    if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND C_COMPILER_VERSION VERSION_EQUAL "11.1")
        message(FATAL_ERROR "You are using Intel compiler version 11.1, and that compiler is known to produce incorrect results with SSE4.1 acceleration. You need to use another compiler (e.g. icc 12 or newer) or different acceleration (probably slower simulations).")
    endif()
elseif(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" OR ${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")

    # Set the AVX compiler flag for both these choices!

    GMX_TEST_CFLAG(GNU_AVX_CFLAG "-mavx" ACCELERATION_C_FLAGS)
    if (NOT GNU_AVX_CFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CFLAG(MSVC_AVX_CFLAG "/arch:AVX" ACCELERATION_C_FLAGS)
    endif (NOT GNU_AVX_CFLAG AND GMX_NATIVE_WINDOWS)
    if (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)
        message(WARNING "No C AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance) giving the -DGMX_CPU_ACCELERATION=SSE4.1 to cmake.")
    endif (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)

    if (CMAKE_CXX_COMPILER_LOADED)
        GMX_TEST_CXXFLAG(GNU_AVX_CXXFLAG "-mavx" ACCELERATION_CXX_FLAGS)
        if (NOT GNU_AVX_CXXFLAG AND GMX_NATIVE_WINDOWS)
            GMX_TEST_CXXFLAG(MSVC_AVX_CXXFLAG "/arch:AVX" ACCELERATION_CXX_FLAGS)
        endif (NOT GNU_AVX_CXXFLAG AND GMX_NATIVE_WINDOWS)
        if (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)
            message(WARNING "No C++ AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance) giving the -DGMX_CPU_ACCELERATION=SSE4.1 to cmake.")
        endif (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)
    endif()

    # Set the FMA4 flags (MSVC doesn't require any)
    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" AND NOT MSVC)
        GMX_TEST_CFLAG(GNU_FMA_CFLAG "-mfma4" ACCELERATION_C_FLAGS)
        if (NOT GNU_FMA_CFLAG)
            message(WARNING "No C FMA4 flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
        endif(NOT GNU_FMA_CFLAG)
        GMX_TEST_CFLAG(GNU_XOP_CFLAG "-mxop" ACCELERATION_C_FLAGS)
        # No big deal if we do not have xop, so no point yelling warnings about it.
        if (CMAKE_CXX_COMPILER_LOADED)
            GMX_TEST_CXXFLAG(GNU_FMA_CXXFLAG "-mfma4" ACCELERATION_CXX_FLAGS)
            if (NOT GNU_FMA_CXXFLAG)
                message(WARNING "No C++ FMA flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
            endif (NOT GNU_FMA_CXXFLAG)
            GMX_TEST_CXXFLAG(GNU_XOP_CXXFLAG "-mxop" ACCELERATION_CXX_FLAGS)
            # No big deal if we do not have xop, so no point yelling warnings about it.
        endif()
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
        set(GMX_CPU_ACCELERATION_X86_AVX_128_FMA 1)
        set(GMX_X86_AVX_128_FMA 1)
        if (NOT ACCELERATION_QUIETLY)
          message(STATUS "Enabling 128-bit AVX Gromacs acceleration (with fused-multiply add), and it will help compiler optimization.")
        endif()

        # We don't have the full compiler version string yet (BUILD_C_COMPILER),
        # so we can't distinguish vanilla and Apple clang, but catering for AMD
        # hackintoshes is not worth the effort.
        if (APPLE AND ${CMAKE_C_COMPILER_ID} STREQUAL "Clang")
            message(WARNING "Due to a known compiler bug, Clang up to version 3.2 (and Apple Clang up to version 4.1) produces incorrect code with AVX_128_FMA acceleration. As we can not work around this bug on OS X, you will have to select a different compiler or CPU acceleration.")
        endif()

        if (GMX_USE_CLANG_FMA_BUG_WORKAROUND)
            # we assume that we have an external assembler that supports AVX
            message(STATUS "Clang ${C_COMPILER_VERSION} detected, enabling FMA bug workaround")
            set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -no-integrated-as")
        endif()

    else()
        # If we are not doing AVX_128, it must be AVX_256...
        set(GMX_CPU_ACCELERATION_X86_AVX_256 1)
        set(GMX_X86_AVX_256 1)
        if (NOT ACCELERATION_QUIETLY)
          message(STATUS "Enabling 256-bit AVX Gromacs acceleration, and it will help compiler optimization.")
        endif()
    endif()

    # Unfortunately gcc-4.5.2 and gcc-4.6.0 has a bug where they use the wrong datatype for the formal
    # parameter of the mask for maskload/maskstore arguments. Check if this is present, since we can work around it.
    gmx_test_avx_gcc_maskload_bug(${ACCELERATION_C_FLAGS} GMX_X86_AVX_GCC_MASKLOAD_BUG)

elseif(${GMX_CPU_ACCELERATION} STREQUAL "IBM_QPX")
    # Used on BlueGene/Q
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        GMX_TEST_CFLAG(XLC_BLUEGENEQ_CFLAG "-qarch=qp -qtune=qp" ACCELERATION_C_FLAGS)
        try_compile(TEST_QPX ${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/cmake/TestQPX.c"
            COMPILE_DEFINITIONS "${ACCELERATION_C_FLAGS}")
        if(NOT TEST_QPX)
            message(FATAL_ERROR "Cannot compile the requested IBM QPX intrinsics.")
        endif()
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "XL" AND CMAKE_CXX_COMPILER_LOADED)
        GMX_TEST_CXXFLAG(XLC_BLUEGENEQ_CXXFLAG "-qarch=qp -qtune=qp" ACCELERATION_CXX_FLAGS)
        try_compile(TEST_QPX ${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/cmake/TestQPX.c"
            COMPILE_DEFINITIONS "${ACCELERATION_CXX_FLAGS}")
        if(NOT TEST_QPX)
            message(FATAL_ERROR "Cannot compile the requested IBM QPX intrinsics.")
        endif()
    endif()

    if (TEST_QPX)
        message(WARNING "IBM QPX acceleration was selected and could be compiled, but the accelerated kernels are not yet available.")
        set(GMX_CPU_ACCELERATION_IBM_QPX 1)
    else()
        message(FATAL_ERROR "Cannot compile IBM QPX intrinsics without the XL compiler. If you are compiling for BlueGene/Q, use 'cmake .. -DCMAKE_TOOLCHAIN_FILE=BlueGeneQ-static-XL-C' to set up the tool chain.")
    endif()
elseif(${GMX_CPU_ACCELERATION} STREQUAL "SPARC64_HPC_ACE")
    set(GMX_CPU_ACCELERATION_SPARC64_HPC_ACE 1)
else(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    MESSAGE(FATAL_ERROR "Unrecognized option for accelerated kernels: ${GMX_CPU_ACCELERATION}. Pick one of None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, IBM_QPX, Sparc64_HPC_ACE")
endif(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
set(ACCELERATION_QUIETLY TRUE CACHE INTERNAL "")

# Process QM/MM Settings
string(TOUPPER ${GMX_QMMM_PROGRAM} GMX_QMMM_PROGRAM)
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
else(${GMX_QMMM_PROGRAM} STREQUAL "GAUSSIAN")
    MESSAGE(FATAL_ERROR "Invalid QM/MM program option: ${GMX_QMMM_PROGRAM}. Choose one of: Gaussian, Mopac, Gamess, Orca, None")
endif(${GMX_QMMM_PROGRAM} STREQUAL "GAUSSIAN")

# Process FFT library settings
string(TOUPPER ${GMX_FFT_LIBRARY} GMX_FFT_LIBRARY)
set(PKG_FFT "")
set(PKG_FFT_LIBS "")
set(MKL_LIBRARIES_FORMAT_DESCRIPTION "Use full paths to library files, in the right order, and separated by semicolons.")
set(MKL_LIBRARIES "" CACHE STRING "List of libraries for linking to MKL. Only used with GMX_FFT_LIBRARY=mkl. ${MKL_LIBRARIES_FORMAT_DESCRIPTION}")
set(MKL_INCLUDE_DIR "" CACHE PATH "Path to mkl.h (non-inclusive). Only used with GMX_FFT_LIBRARY=mkl.")
mark_as_advanced(FORCE MKL_LIBRARIES)
mark_as_advanced(FORCE MKL_INCLUDE_DIR)
if(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    if(GMX_DOUBLE)
        set(FFTW fftw)
    else()
        set(FFTW fftwf)
    endif()

    if(GMX_BUILD_OWN_FFTW)
      add_subdirectory(src/contrib/fftw)
    else()
      find_package(FFTW COMPONENTS ${FFTW})
    endif()

    string(TOUPPER "${FFTW}" FFTW)
    if(NOT ${FFTW}_FOUND)
      MESSAGE(FATAL_ERROR "Cannot find FFTW 3 (with correct precision - libfftw3f for single-precision GROMACS or libfftw3 for double-precision GROMACS). Either choose the right precision, choose another FFT(W) library (-DGMX_FFT_LIBRARY), enable the advanced option to let GROMACS build FFTW 3 for you (-DGMX_BUILD_OWN_FFTW=ON) , or use the really slow GROMACS built-in fftpack library (-DGMX_FFT_LIBRARY=fftpack).")
    endif()

    set(PKG_FFT "${${FFTW}_PKG}")
    include_directories(${${FFTW}_INCLUDE_DIRS})
    set(FFT_LIBRARIES ${${FFTW}_LIBRARIES})
    set(GMX_FFT_FFTW3 1)

    if ((${GMX_CPU_ACCELERATION} MATCHES "SSE" OR ${GMX_CPU_ACCELERATION} MATCHES "AVX") AND NOT ${FFTW}_HAVE_SIMD)
      message(WARNING "The fftw library found is compiled without SIMD support, which makes it slow. Consider recompiling it or contact your admin")
    endif()

    if((${GMX_CPU_ACCELERATION} MATCHES "SSE" OR ${GMX_CPU_ACCELERATION} MATCHES "AVX") AND ${FFTW}_HAVE_AVX)
        # If we're not doing CPU acceleration, we don't care about FFTW performance on x86 either
        message(WARNING "The FFTW library was compiled with --enable-avx to enable AVX SIMD instructions. That might sound like a good idea for your processor, but for FFTW versions up to 3.3.3, these are slower than the SSE/SSE2 SIMD instructions for the way GROMACS uses FFTs. Limitations in the way FFTW allows GROMACS to measure performance make it awkward for either GROMACS or FFTW to make the decision for you based on runtime performance. You should compile a different FFTW library with --enable-sse or --enable-sse2. If you have a more recent FFTW, you may like to compare the performance of GROMACS with FFTW libraries compiled with and without --enable-avx. However, the GROMACS developers do not really expect the FFTW AVX optimization to help, because the performance is limited by memory access, not computation.")
    endif()

    if(NOT "${GMX_FFT_LIBRARY}" STREQUAL "${GMX_FFT_LIBRARY_PREVIOUS_VALUE}")
        MESSAGE(STATUS "Using external FFT library - FFTW3")
    endif()
elseif(${GMX_FFT_LIBRARY} STREQUAL "MKL")
    # Intel 11 and up makes life somewhat easy if you just want to use
    # all their stuff. It's not easy if you only want some of their
    # stuff...
    if (CMAKE_C_COMPILER_ID MATCHES "Intel" AND C_COMPILER_VERSION VERSION_GREATER "11")
        # The next line takes care of everything for MKL
        if (WIN32)
            set(FFT_LINKER_FLAGS "/Qmkl=sequential")
        else()
            set(FFT_LINKER_FLAGS "-mkl=sequential")
        endif()
        # Some versions of icc require this in order that mkl.h can be
        # found at compile time.
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} ${FFT_LINKER_FLAGS}")

        set(MKL_ERROR_MESSAGE "Make sure you have configured your compiler so that ${FFT_LINKER_FLAGS} will work.")
    else()
        # The user will have to provide the set of magic libraries in
        # MKL_LIBRARIES, which we cache (non-advanced), so that they
        # don't have to keep specifying it, and can easily see that
        # CMake is still using that information.
        set(MKL_LIBRARIES "${MKL_LIBRARIES}" CACHE STRING "User-specified libraries for linking to MKL")
        mark_as_advanced(CLEAR MKL_LIBRARIES)
        include_directories(${MKL_INCLUDE_DIR})
        mark_as_advanced(CLEAR MKL_INCLUDE_DIR)
        set(FFT_LIBRARIES "${MKL_LIBRARIES}")
        set(MKL_ERROR_MESSAGE "The include path to mkl.h in MKL_INCLUDE_DIR, and the link libraries in MKL_LIBRARIES=${MKL_LIBRARIES} need to match what the MKL documentation says you need for your system. ${MKL_LIBRARIES_FORMAT_DESCRIPTION}")
        # Convert the semi-colon separated list to a list of
        # command-line linker arguments so that code using our
        # pkgconfig setup can use it.
        string(REGEX REPLACE ";" " " PKG_FFT_LIBS "${MKL_LIBRARIES}")
    endif()

    # Check MKL works. If we were in a non-global scope, we wouldn't
    # have to play nicely.
    set(old_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FFT_LINKER_FLAGS}")
    set(old_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES "${FFT_LIBRARIES}")

    check_function_exists(DftiCreateDescriptor TEST_MKL)

    set(CMAKE_REQUIRED_FLAGS "${old_CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_LIBRARIES "${old_CMAKE_REQUIRED_LIBRARIES}")

    if(NOT TEST_MKL)
        # Hack to help the user vary MKL settings until they work.
        # TODO Make this logic more useful.
        unset(TEST_MKL CACHE)
        message(FATAL_ERROR "Linking with MKL was requested, but was not successful. ${MKL_ERROR_MESSAGE}")
    endif()

    # Set variables to signal that we have MKL available and should use it for FFTs.
    set(GMX_FFT_MKL 1)
    set(HAVE_LIBMKL 1)

    if(NOT "${GMX_FFT_LIBRARY}" STREQUAL "${GMX_FFT_LIBRARY_PREVIOUS_VALUE}")
        MESSAGE(STATUS "Using external FFT library - Intel MKL")
    endif()

#elseif(${GMX_FFT_LIBRARY} STREQUAL "ACML")
#    MESSAGE(STATUS "Using external FFT library - AMD core math library")
#    set(GMX_FFT_ACML 1)
elseif(${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
    set(GMX_FFT_FFTPACK 1)
    if(NOT "${GMX_FFT_LIBRARY}" STREQUAL "${GMX_FFT_LIBRARY_PREVIOUS_VALUE}")
        MESSAGE(STATUS "Using internal FFT library - fftpack")
    endif()
else(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    MESSAGE(FATAL_ERROR "Invalid FFT library setting: ${GMX_FFT_LIBRARY}. Choose one of: fftw3, mkl, fftpack")
endif(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
set(GMX_FFT_LIBRARY_PREVIOUS_VALUE ${GMX_FFT_LIBRARY} CACHE INTERNAL "Previous value for GMX_FFT_LIBRARY, so that we can detect when it changes.")

# enable threaded fftw3 if we've found it 
if(FFTW3_THREADS OR FFTW3F_THREADS)
    add_definitions(-DFFT5D_FFTW_THREADS)
endif()

include(gmxManageLinearAlgebraLibraries)

set(GMX_USE_PLUGINS OFF CACHE INTERNAL "Whether GROMACS will really try to compile support for VMD plugins")

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

if(EXISTS "${CMAKE_SOURCE_DIR}/admin/.isreposource")
    if(NOT CMAKE_CROSSCOMPILING)
        option(GMX_BUILD_MANPAGES "Build man pages" ON)
    else()
        message(STATUS "Building the man pages is not available when "
            "cross-compiling the developer version from git")
    endif()
else()
    #make sure source package contains all man pages
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/man/man1/ngmx.1")
        message(FATAL_ERROR "Man pages are missing from source package.")
    endif()
endif()
mark_as_advanced(GMX_BUILD_MANPAGES)

# Math and thread libraries must often come after all others when linking...
if(HAVE_LIBM)
    list(APPEND	GMX_EXTRA_LIBRARIES m)
endif(HAVE_LIBM)

if(GMX_FAHCORE)
  set(COREWRAP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../corewrap" CACHE STRING 
      "Path to swindirect.h")
  include_directories(${COREWRAP_INCLUDE_DIR})
  set_property(CACHE GMX_COOL_QUOTES PROPERTY VALUE OFF)
endif(GMX_FAHCORE)

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (NOT GMX_SKIP_DEFAULT_CFLAGS)
    set(CMAKE_C_FLAGS "${ACCELERATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${ACCELERATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${FFT_LINKER_FLAGS} ${MPI_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${MPI_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
else()
    message("Recommended flags which are not added because GMX_SKIP_DEFAULT_CFLAGS=yes:")
    message("CMAKE_C_FLAGS: ${ACCELERATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${GMXC_CFLAGS}")
    message("CMAKE_C_FLAGS_RELEASE: ${GMXC_CFLAGS_RELEASE}")
    message("CMAKE_C_FLAGS_DEBUG: ${GMXC_CFLAGS_DEBUG}")
    if(CMAKE_CXX_COMPILER_LOADED)
        message("CMAKE_CXX_FLAGS: ${ACCELERATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${GMXC_CXXFLAGS}")
        message("CMAKE_CXX_FLAGS_RELEASE: ${GMXC_CXXFLAGS_RELEASE}")
        message("CMAKE_CXX_FLAGS_DEBUG: ${GMXC_CXXFLAGS_DEBUG}")
    endif()
    message("CMAKE_EXE_LINKER_FLAGS: ${FFT_LINKER_FLAGS} ${MPI_LINKER_FLAGS}")
    message("CMAKE_SHARED_LINKER_FLAGS: ${MPI_LINKER_FLAGS}")
endif()

if(NOT GMX_OPENMP)
    #Unset all OpenMP flags in case OpenMP was disabled either by the user
    #or because it was only partially detected (e.g. only for C but not C++ compiler)
    unset(OpenMP_C_FLAGS CACHE) 
    unset(OpenMP_CXX_FLAGS CACHE)
    unset(OpenMP_LINKER_FLAGS CACHE)
    unset(OpenMP_SHARED_LINKER_FLAGS)
endif()
set(PKG_CFLAGS "${PKG_CFLAGS} ${OpenMP_C_FLAGS}")

######################################
# Output compiler and CFLAGS used
######################################
get_compiler_info(C BUILD_C_COMPILER BUILD_CFLAGS)
if (CMAKE_CXX_COMPILER_LOADED)
    get_compiler_info(CXX BUILD_CXX_COMPILER BUILD_CXXFLAGS)
endif ()
if(GMX_GPU)
    get_cuda_compiler_info(CUDA_NVCC_COMPILER_INFO CUDA_NVCC_COMPILER_FLAGS)
endif(GMX_GPU)


########################################################################
# Specify install locations and which subdirectories to process        #
########################################################################
if (GMX_USE_RELATIVE_INSTALL_PATH)
    set(GMX_INSTALL_PREFIX "" CACHE STRING "Prefix gets appended to CMAKE_INSTALL_PREFIX. For cpack it sets the root folder of the archive.")
    mark_as_advanced(GMX_INSTALL_PREFIX)
else()
    set(GMX_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/")
endif()

if ( NOT DEFINED GMXLIB )
    set(GMXLIB lib)
endif()
set(LIB_INSTALL_DIR "${GMX_INSTALL_PREFIX}${GMXLIB}")
set(BIN_INSTALL_DIR  ${GMX_INSTALL_PREFIX}bin)
set(DATA_INSTALL_DIR ${GMX_INSTALL_PREFIX}share/gromacs)
set(MAN_INSTALL_DIR  ${GMX_INSTALL_PREFIX}share/man)
set(INCL_INSTALL_DIR ${GMX_INSTALL_PREFIX}include)

set(GMXLIBDIR        ${DATA_INSTALL_DIR}/top)

########################################################################
# Set up binary and library suffixing
########################################################################
set(GMX_BINARY_SUFFIX "" CACHE STRING "Suffix for GROMACS binaries (default: _d for double, _mpi for MPI, _mpi_d for MPI and double).")
set(GMX_LIBS_SUFFIX ""
  CACHE STRING "Suffix for GROMACS libs (default: _d for double, _mpi for MPI, _mpi_d for MPI and double).")
if (GMX_DEFAULT_SUFFIX)
  set(GMX_BINARY_SUFFIX "")
  set(GMX_LIBS_SUFFIX "")
  if (GMX_LIB_MPI)
    set(GMX_BINARY_SUFFIX "_mpi")
    set(GMX_LIBS_SUFFIX "_mpi")
  endif()
  if (GMX_DOUBLE)
    set (GMX_BINARY_SUFFIX "${GMX_BINARY_SUFFIX}_d")
    set (GMX_LIBS_SUFFIX "${GMX_LIBS_SUFFIX}_d")
  endif(GMX_DOUBLE)
  mark_as_advanced(FORCE GMX_BINARY_SUFFIX GMX_LIBS_SUFFIX)
  if (NOT SUFFIX_QUIETLY)
    message(STATUS "Using default binary suffix: \"${GMX_BINARY_SUFFIX}\"")
    message(STATUS "Using default library suffix: \"${GMX_LIBS_SUFFIX}\"")
  endif (NOT SUFFIX_QUIETLY)
else(GMX_DEFAULT_SUFFIX)
  mark_as_advanced(CLEAR GMX_BINARY_SUFFIX GMX_LIBS_SUFFIX)
  if (NOT SUFFIX_QUIETLY)
    message(STATUS "Using manually set binary suffix: \"${GMX_BINARY_SUFFIX}\"")
    message(STATUS "Using manually set library suffix: \"${GMX_LIBS_SUFFIX}\"")
  endif (NOT SUFFIX_QUIETLY)
endif(GMX_DEFAULT_SUFFIX)
set(SUFFIX_QUIETLY TRUE CACHE INTERNAL "")

##################################################################
# Shared library settings - Darwin uses INSTALL_NAME_DIR instead!
##################################################################
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "\\\$ORIGIN/../${GMXLIB}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

#COPYING file: Only necessary for binary distributions.
#Simpler to always install.
install(FILES COPYING DESTINATION ${DATA_INSTALL_DIR} COMPONENT data)

add_subdirectory(share)
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(scripts)

# Issue a warning if NVIDIA GPUs were detected, but CUDA was not found.
# Don't bother the user after the first configure pass.
if ((CUDA_NOTFOUND_AUTO AND GMX_DETECT_GPU_AVAILABLE) AND NOT GMX_GPU_DETECTION_DONE)
    message(WARNING "${CUDA_NOTFOUND_MESSAGE}")
    unset(CUDA_NOTFOUND_AUTO)
    unset(CUDA_NOTFOUND_MESSAGE)
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
# Tests                                                                #
########################################################################

include(CTest)
mark_as_advanced(BUILD_TESTING)
#gmxtests target builds all binaries required for running gmxtest
add_custom_target(gmxtests DEPENDS grompp mdrun pdb2gmx gmxcheck editconf)
IF(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
    if(REGRESSIONTEST_PATH)
        #check target builds all to run tests and the runs tests
        add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
        add_dependencies(check gmxtests)
    else()
        add_custom_target(check COMMAND ${CMAKE_COMMAND} -E echo "WARNING: No tests are run. Running the tests requires either of the cmake variables REGRESSIONTEST_PATH or REGRESSIONTEST_DOWNLOAD to be set.")
    endif()
ENDIF()

