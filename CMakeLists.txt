cmake_minimum_required(VERSION 2.8)
# Keep CMake suitably quiet on Cygwin
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform)

project(Gromacs)
include(Dart)
mark_as_advanced(DART_ROOT)

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
set(PROJECT_VERSION "5.0-dev")
# The version number of the regressiontest tarball against which this
# git branch can be tested. Normally, this will be the version of the
# last patch release. Comment the next line out for branches leading
# to a major/minor release.
set(REGRESSIONTEST_VERSION "4.6.1")
set(CUSTOM_VERSION_STRING ""
    CACHE STRING "Custom version string (if empty, use hard-coded default)")
mark_as_advanced(CUSTOM_VERSION_STRING)
if (CUSTOM_VERSION_STRING)
    set(PROJECT_VERSION ${CUSTOM_VERSION_STRING})
endif (CUSTOM_VERSION_STRING)
set(SOVERSION 8)
# It is a bit irritating, but this has to be set separately for now!
SET(CPACK_PACKAGE_VERSION_MAJOR "5")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
#SET(CPACK_PACKAGE_VERSION_PATCH "0")

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

set(GMX_INSTALL_PREFIX "" CACHE STRING "Prefix gets appended to CMAKE_INSTALL_PREFIX. For cpack it sets the root folder of the archive.")
mark_as_advanced(GMX_INSTALL_PREFIX)

include(gmxBuildTypeReference)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Reference." FORCE)
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
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_SOURCE_DIR};/;${CMAKE_BINARY_DIR}/man;man")
set(CPACK_PACKAGE_CONTACT "gmx-users@gromacs.org")

#must come after all cpack settings!
include(CPack)

set(SOURCE_IS_GIT_REPOSITORY OFF)
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(SOURCE_IS_GIT_REPOSITORY ON)
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
option(GMX_COOL_QUOTES "Enable Gromacs cool quotes" ON)
mark_as_advanced(GMX_COOL_QUOTES)

set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")
# User input options                                                   #
########################################################################
option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" OFF)
option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREAD_MPI  "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" ON)
option(GMX_SOFTWARE_INVSQRT "Use GROMACS software 1/sqrt" ON)
mark_as_advanced(GMX_SOFTWARE_INVSQRT)
option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
mark_as_advanced(GMX_FAHCORE)

# decide on GPU settings based on user-settings and GPU/CUDA detection
include(gmxManageGPU)

# TODO: move OpenMM to contrib
option(GMX_OPENMM "Accelerated execution on GPUs through the OpenMM library (rerun cmake after changing to see relevant options)" OFF)
mark_as_advanced(GMX_OPENMM)

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
    CACHE STRING "Accelerated CPU kernels. Pick one of: None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, IBM_QPX")

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

option(GMX_GPU  "Enable GPU acceleration" ON)
option(GMX_OPENMP "Enable OpenMP-based multithreading" ON)

option(GMX_GIT_VERSION_INFO "Generate git version information" ${SOURCE_IS_GIT_REPOSITORY})
mark_as_advanced(GMX_GIT_VERSION_INFO)

option(GMX_DEFAULT_SUFFIX "Use default suffixes for GROMACS binaries and libs (_d for double, _mpi for MPI; rerun cmake after changing to see relevant options)" ON)
if(UNIX)
    option(GMX_SYMLINK_OLD_BINARY_NAMES "Create symbolic links for pre-5.0 binary names" ON)
endif()

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
include(CheckCXXCompilerFlag)

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
            if(CMAKE_COMPILER_IS_GNUCC AND GMX_PREFER_STATIC_OPENMP AND NOT APPLE)
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

if (CMAKE_C_COMPILER_ID STREQUAL "PGI")
    message(WARNING "All tested PGI compiler versions (up to 12.9.0) generate binaries which produce incorrect results, or even fail to compile Gromacs. Highly recommended to use a different compiler. If you choose to use PGI, make sure to run the regressiontests.")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND C_COMPILER_VERSION VERSION_LESS "12.0.0")
    message(WARNING "Intel compilers before 12.0.0 are not routinely tested, so there may be problems. Version 11.1 with SSE4.1 is known to produce incorrect results. It is highly recommended to use a more up-to-date compiler. If you choose to use this version, make sure you run the regressiontests.")
endif()

########################################################################
# Set up binary and library suffixing 
########################################################################
set(GMX_BINARY_SUFFIX "" CACHE STRING "Suffix for GROMACS binaries (default: _d for double, _mpi for MPI, _mpi_d for MPI and double).")
set(GMX_LIBS_SUFFIX "" 
  CACHE STRING "Suffix for GROMACS libs (default: _d for double, _mpi for MPI, _mpi_d for MPI and double).")
if (GMX_DEFAULT_SUFFIX)
  set(GMX_BINARY_SUFFIX "")
  set(GMX_LIBS_SUFFIX "")
  if (GMX_MPI)
    set(GMX_BINARY_SUFFIX "_mpi")
    set(GMX_LIBS_SUFFIX "_mpi")
  endif(GMX_MPI)
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
include(CheckIncludeFileCXX)
check_include_files(unistd.h	 HAVE_UNISTD_H)
check_include_files(pwd.h        HAVE_PWD_H)
check_include_files(pthread.h    HAVE_PTHREAD_H)
check_include_files(dirent.h     HAVE_DIRENT_H)
check_include_files(sys/time.h   HAVE_SYS_TIME_H)
check_include_files(io.h  		 HAVE_IO_H)
check_include_files(sched.h      HAVE_SCHED_H)

check_include_files(regex.h      HAVE_POSIX_REGEX)
check_include_file_cxx(regex     HAVE_CXX11_REGEX)
# TODO: It could be nice to inform the user if no regex support is found,
# as selections won't be fully functional.

include(CheckFunctionExists)
check_function_exists(strdup            HAVE_STRDUP)
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
check_library_exists(m cbrt "" HAVE_CBRT)


include(CheckTypeSize)

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
  # This makes windows.h not declare min/max as macros that would break
  # C++ code using std::min/std::max.
  add_definitions(-DNOMINMAX)

  if (NOT BUILD_SHARED_LIBS)
      option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer static system libraries (MT instead of MD)!" ON)
      if(NOT GMX_PREFER_STATIC_LIBS)
          message(WARNING "Shared system libraries requested, and static Gromacs libraries requested.")
      endif()
  else()
      message(FATAL_ERROR "BUILD_SHARED_LIBS not yet working for Windows in the master branch")
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

# Unconditionally find the package, as it is also required for unit tests
find_package(LibXml2)
option(GMX_XML "Use libxml2 to parse xml files (currently has no effect)" ${LIBXML2_FOUND})
set(PKG_XML "")
mark_as_advanced(GMX_XML)
# Don't actually do anything, since libxml2 is currently not used by libgromacs
#if(GMX_XML AND NOT LIBXML2_FOUND)
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
    list(APPEND	GMX_EXTRA_LIBRARIES gsl)
  endif(GSL_FOUND)
endif (GMX_GSL)

option(GMX_CLN "Add support for cln" ON)
if (GMX_CLN)
  find_package(cln)
  set(PKG_CLN "")
  if(CLN_FOUND)
    include_directories(${CLN_INCLUDE_DIR})
    set(PKG_CLN cln)
    set(HAVE_LIBCLN 1)
    list(APPEND	GMX_EXTRA_LIBRARIES ${CLN_LIBRARIES})
  endif(CLN_FOUND)
endif (GMX_CLN)

option(GMX_SQLITE3 "Add support for sqlite3" ON)
if (GMX_SQLITE3)
  find_package(sqlite3)
  set(PKG_SQLITE3 "")
  if(SQLITE3_FOUND)
    include_directories(${SQLITE3_INCLUDE_DIR})
    set(PKG_SQLITE3sqlite3)
    set(HAVE_LIBSQLITE3 1)
    list(APPEND	GMX_EXTRA_LIBRARIES ${SQLITE3_LIBRARIES})
  endif(SQLITE3_FOUND)
endif (GMX_SQLITE3)

option(GMX_OPENBABEL2 "Add support for OpenBabel2 functionality" OFF)
if (GMX_OPENBABEL2)
  set(PKG_OPENBABEL2 "")
  find_package(OpenBabel2 REQUIRED)
  if (OPENBABEL2_FOUND)
    include_directories(${OPENBABEL2_INCLUDE_DIR})
    set(PKG_OPENBABEL2 openbabel2)
    set(HAVE_LIBOPENBABEL2 1)
    list(APPEND	GMX_EXTRA_LIBRARIES ${OPENBABEL2_LIBRARIES})
  endif(OPENBABEL2_FOUND)
endif (GMX_OPENBABEL2)

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
if(GMX_THREAD_MPI)
    tmpi_get_source_list(THREAD_MPI_SRC CXX)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_THREAD_MPI")
    set(GMX_MPI 1)
else(GMX_THREAD_MPI)
    tmpi_get_source_list(THREAD_MPI_SRC CXX NOMPI)
endif(GMX_THREAD_MPI)

if(GMX_GPU)
    # now that we have detected the dependencies, do the second configure pass
    gmx_gpu_setup()
endif(GMX_GPU)

if(APPLE)
   find_library(ACCELERATE_FRAMEWORK Accelerate)
   list(APPEND GMX_EXTRA_LIBRARIES ${ACCELERATE_FRAMEWORK})
endif(APPLE)

if(CYGWIN)
    set(GMX_CYGWIN 1)
endif(CYGWIN)

if(WIN32 AND NOT CYGWIN)
    set(GMX_NATIVE_WINDOWS 1)
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

option(GMX_BUILD_UNITTESTS "Build unit tests with BUILD_TESTING (uses Google C++ Testing and Mocking Frameworks, requires libxml2)" ${LIBXML2_FOUND})
mark_as_advanced(GMX_BUILD_UNITTESTS)
if (GMX_BUILD_UNITTESTS AND NOT LIBXML2_FOUND)
    message(FATAL_ERROR
        "Cannot build unit tests without libxml2. "
        "Either set GMX_BUILD_UNITTESTS=OFF or tell CMake how to find libxml2.")
endif()
set(MEMORYCHECK_SUPPRESSIONS_FILE ${CMAKE_SOURCE_DIR}/cmake/legacy_and_external.supp)

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

include(gmxTestCXX11)
gmx_test_cxx11(GMX_CXX11 GMX_CXX11_FLAGS)
if(CXX11_FLAG AND GMX_GPU)
    #FIXME: add proper solution for progate all but cxx11 flag
    set(CUDA_PROPAGATE_HOST_FLAGS no)
    message(WARNING "Please manually add compiler flags to CUDA_NVCC_FLAGS. Automatic propogation temporary not working.")
endif()

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
        message(WARNING "No C AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
    endif (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)

    GMX_TEST_CXXFLAG(GNU_AVX_CXXFLAG "-mavx" ACCELERATION_CXX_FLAGS)
    if (NOT GNU_AVX_CXXFLAG AND GMX_NATIVE_WINDOWS)
        GMX_TEST_CXXFLAG(MSVC_AVX_CXXFLAG "/arch:AVX" ACCELERATION_CXX_FLAGS)
    endif (NOT GNU_AVX_CXXFLAG AND GMX_NATIVE_WINDOWS)
    if (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)
       message(WARNING "No C++ AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
    endif (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)

    # Set the FMA4 flags (MSVC doesn't require any)
    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" AND NOT MSVC)
        if (${CMAKE_COMPILER_ID} MATCHES "Clang")
            message(FATAL_ERROR "Clang up to at least version 3.2 produces incorrect code for AVX_128_FMA. Sorry, but you will have to select a different compiler or acceleration.")
        endif()
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
            "cmake/TestQPX.c"
            COMPILE_DEFINITIONS "${ACCELERATION_CXX_FLAGS")
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

else(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    MESSAGE(FATAL_ERROR "Unrecognized option for accelerated kernels: ${GMX_CPU_ACCELERATION}. Pick one of None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, IBM_QPX")
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
    set(PKG_FFT "${${FFTW}_PKG}")
    include_directories(${${FFTW}_INCLUDE_DIRS})
    set(FFT_LIBRARIES ${${FFTW}_LIBRARIES})
    if(NOT ${FFTW}_FOUND)
      MESSAGE(FATAL_ERROR "Cannot find FFTW 3 (with correct precision - libfftw3f for single-precision GROMACS or libfftw3 for double-precision GROMACS). Either choose the right precision, choose another FFT(W) library, enable the advanced option to let GROMACS build FFTW 3 for you, or use the really slow GROMACS built-in fftpack library.")
    endif()

    set(GMX_FFT_FFTW3 1)

    if (NOT ${GMX_CPU_ACCELERATION} STREQUAL "NONE" AND NOT ${FFTW}_HAVE_SIMD) 
      message(WARNING "The fftw library found is compiled without SIMD support, which makes it slow. Consider recompiling it or contact your admin")
    endif()

    if(NOT ${GMX_CPU_ACCELERATION} STREQUAL "NONE" AND ${FFTW}_HAVE_AVX)
        # If we're not doing CPU acceleration, we don't care about FFTW performance on x86 either
        message(WARNING "The FFTW library was compiled with --enable-avx to enable AVX SIMD instructions. That might sound like a good idea for your processor, but for FFTW versions up to 3.3.3, these are slower than the SSE/SSE2 SIMD instructions for the way GROMACS uses FFTs. Limitations in the way FFTW allows GROMACS to measure performance make it awkward for either GROMACS or FFTW to make the decision for you based on runtime performance. You should compile a different FFTW library with --enable-sse or --enable-sse2. If you have a more recent FFTW, you may like to compare the performance of GROMACS with FFTW libraries compiled with and without --enable-avx. However, the GROMACS developers do not really expect the FFTW AVX optimization to help, because the performance is limited by memory access, not computation.")
    endif()

elseif(${GMX_FFT_LIBRARY} STREQUAL "MKL")
#    MESSAGE(STATUS "Using external FFT library - Intel MKL")
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIR})
    set(FFT_LIBRARIES ${MKL_LIBRARIES})
    set(PKG_FFT_LIBS ${MKL_LIBRARIES})

    set(GMX_FFT_MKL 1)
    set(HAVE_MKL 1)

#elseif(${GMX_FFT_LIBRARY} STREQUAL "ACML")
#    MESSAGE(STATUS "Using external FFT library - AMD core math library")
#    set(GMX_FFT_ACML 1)
elseif(${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
    MESSAGE(STATUS "Using internal FFT library - fftpack")
    set(GMX_FFT_FFTPACK 1)
else(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    MESSAGE(FATAL_ERROR "Invalid FFT library setting: ${GMX_FFT_LIBRARY}. Choose one of: fftw3, mkl, fftpack")
endif(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")

# enable threaded fftw3 if we've found it 
if(FFTW3_THREADS OR FFTW3F_THREADS)
    add_definitions(-DFFT5D_FFTW_THREADS)
endif()

set(GMX_EXTERNAL_BLAS TRUE CACHE BOOL "Use external BLAS instead of built-in")
set(GMX_EXTERNAL_LAPACK TRUE CACHE BOOL "Use external LAPACK instead of built-in")
# MKL has BLAS/LAPACK routines
if(NOT HAVE_MKL AND NOT ACCELERATE_FRAMEWORK)
  if(GMX_EXTERNAL_BLAS)
    if (GMX_BLAS_USER)
        list(APPEND GMX_EXTRA_LIBRARIES ${GMX_BLAS_USER})
    else(GMX_BLAS_USER)
        set(BLAS_FIND_QUIETLY ON)
        find_package(BLAS)
        if (BLAS_FOUND)
          list(APPEND GMX_EXTRA_LIBRARIES ${BLAS_LIBRARIES})
        else()
          MESSAGE(STATUS "Using internal BLAS library")
          set(GMX_EXTERNAL_BLAS FALSE CACHE BOOL "Use external BLAS instead of built-in" FORCE)
        endif()
    endif(GMX_BLAS_USER)
  endif(GMX_EXTERNAL_BLAS)
  if(GMX_EXTERNAL_LAPACK)
    if (GMX_LAPACK_USER)
        list(APPEND GMX_EXTRA_LIBRARIES ${GMX_LAPACK_USER})
    else(GMX_LAPACK_USER)
        set(LAPACK_FIND_QUIETLY ON)
        find_package(LAPACK)
        if (LAPACK_FOUND)
          list(APPEND GMX_EXTRA_LIBRARIES ${LAPACK_LIBRARIES})
        else()
          MESSAGE(STATUS "Using internal LAPACK library")
          set(GMX_EXTERNAL_LAPACK FALSE CACHE BOOL "Use external LAPACK instead of built-in" FORCE)
        endif()
    endif(GMX_LAPACK_USER)
  endif(GMX_EXTERNAL_LAPACK)
endif()
mark_as_advanced(GMX_EXTERNAL_LAPACK)
mark_as_advanced(GMX_EXTERNAL_BLAS)

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
  set_property(CACHE GMX_COOL_QUOTES VALUE OFF)
endif(GMX_FAHCORE)

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (NOT GMX_SKIP_DEFAULT_CFLAGS)
    set(CMAKE_C_FLAGS "${ACCELERATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${ACCELERATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${GMX_CXX11_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${MPI_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${MPI_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
else()
    message("Recommended flags which are not added because GMX_SKIP_DEFAULT_CFLAGS=yes:")
    message("CMAKE_C_FLAGS: ${ACCELERATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${GMXC_CFLAGS}")
    message("CMAKE_C_FLAGS_RELEASE: ${GMXC_CFLAGS_RELEASE}")
    message("CMAKE_C_FLAGS_DEBUG: ${GMXC_CFLAGS_DEBUG}")
    message("CMAKE_CXX_FLAGS: ${ACCELERATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${GMX_CXX11_FLAGS} ${GMXC_CXXFLAGS}")
    message("CMAKE_CXX_FLAGS_RELEASE: ${GMXC_CXXFLAGS_RELEASE}")
    message("CMAKE_CXX_FLAGS_DEBUG: ${GMXC_CXXFLAGS_DEBUG}")
    message("CMAKE_EXE_LINKER_FLAGS: ${MPI_LINKER_FLAGS}")
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

######################################
# Output compiler and CFLAGS used
######################################
get_compiler_info(C BUILD_C_COMPILER BUILD_CFLAGS)
if (CMAKE_CXX_COMPILER_LOADED)
    get_compiler_info(CXX BUILD_CXX_COMPILER BUILD_CXXFLAGS)
endif ()

########################################################################
# Specify install locations
########################################################################
if ( NOT DEFINED GMXLIB )
    set(GMXLIB lib)
endif()
set(LIB_INSTALL_DIR "${GMX_INSTALL_PREFIX}${GMXLIB}")
set(BIN_INSTALL_DIR  ${GMX_INSTALL_PREFIX}bin)
set(DATA_INSTALL_DIR ${GMX_INSTALL_PREFIX}share/gromacs)
set(MAN_INSTALL_DIR  ${GMX_INSTALL_PREFIX}share/man)
set(INCL_INSTALL_DIR ${GMX_INSTALL_PREFIX}include)

set(GMXLIBDIR        ${DATA_INSTALL_DIR}/top)

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
    set(PKG_CFLAGS "${PKG_CFLAGS} -DBOOST_NO_TYPEID -I${INCL_INSTALL_DIR}/gromacs/external/boost")
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/external/boost/boost
            DESTINATION ${INCL_INSTALL_DIR}/gromacs/external/boost
            COMPONENT development)
endif()

add_subdirectory(doxygen)
add_subdirectory(share)
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
add_custom_target(gmxtests DEPENDS grompp mdrun pdb2gmx gmxcheck gmx links)
IF(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
    #"check" target builds and runs all tests
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
    add_dependencies(check gmxtests tests)
    #TODO: Add warning if NOT REGRESSIONTEST_PATH OR NOT GMX_XML that regression/unit tests are not run.
ENDIF()

