cmake_minimum_required(VERSION 2.8)
# Keep CMake suitably quiet on Cygwin
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

# Allows CPack to act differently for normal tools and mdrun (e.g. because of MPI)
set(CPACK_COMPONENT_GROUP_TOOLS_DESCRIPTION "All GROMACS executable tools")
set(CPACK_COMPONENT_GROUP_MDRUN_DESCRIPTION "GROMACS executable for running simulations")

# override bugs on OS X where Cmake picks gcc (GNU) for C instead of system default cc (Clang).
if(APPLE)
    set(CMAKE_C_COMPILER_INIT "cc")
endif(APPLE)

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
set(PROJECT_VERSION "4.6-dev")
set(CUSTOM_VERSION_STRING ""
    CACHE STRING "Custom version string (if empty, use hard-coded default)")
mark_as_advanced(CUSTOM_VERSION_STRING)
if (CUSTOM_VERSION_STRING)
    set(PROJECT_VERSION ${CUSTOM_VERSION_STRING})
endif (CUSTOM_VERSION_STRING)
set(SOVERSION 6)
# It is a bit irritating, but this has to be set separately for now!
SET(CPACK_PACKAGE_VERSION_MAJOR "4")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
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

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND UNIX)
set(CMAKE_INSTALL_PREFIX "/usr/local/gromacs" CACHE STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
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
set(CPACK_SOURCE_IGNORE_FILES ".isreposource;.git;.gitignore")
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
option(GMX_GPU  "Enable GPU acceleration" ON)
option(GMX_OPENMM "Accelerated execution on GPUs through the OpenMM library (rerun cmake after changing to see relevant options)" OFF)
option(GMX_FORCE_CXX "Enable C++ compilation even if not necessary" OFF)
mark_as_advanced(GMX_FORCE_CXX)
if(GMX_GPU OR GMX_OPENMM OR GMX_FORCE_CXX)
    enable_language(CXX)
endif()
set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")

########################################################################
# Fix stupid flags on Windows
########################################################################
SET(SHARED_LIBS_DEFAULT ON) 
IF( WIN32 AND NOT CYGWIN)
  option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer static system libraries (MT instead of MD)!" ON)
  mark_as_advanced(GMX_PREFER_STATIC_LIBS)
  SET(SHARED_LIBS_DEFAULT OFF)  #is currently not working on Windows

  IF (GMX_PREFER_STATIC_LIBS)
    #Only setting Debug and Release flags. Others configurations current not used.
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)
    if(CMAKE_CXX_COMPILER_LOADED)
        STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
        SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} CACHE STRING "" FORCE)
        STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
        SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} CACHE STRING "" FORCE)
    endif()
  ENDIF()

  #Workaround for cmake bug 13174. Replace deprecated options.
  IF( CMAKE_C_COMPILER_ID MATCHES "Intel" )
    STRING(REPLACE /GZ /RTC1 CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)
  ENDIF()
  IF( CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_LOADED)
    STRING(REPLACE /GZ /RTC1 CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    STRING(REPLACE /GX /EHsc CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} CACHE STRING "" FORCE)

    STRING(REPLACE /GX /EHsc CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
    SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} CACHE STRING "" FORCE)
  ENDIF()
ENDIF()


########################################################################
# User input options                                                   #
########################################################################
option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" OFF)
option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREAD_MPI  "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" ON)
option(GMX_SOFTWARE_INVSQRT "Use GROMACS software 1/sqrt" ON)
mark_as_advanced(GMX_SOFTWARE_INVSQRT)
option(GMX_POWERPC_INVSQRT "Use PowerPC hardware 1/sqrt" OFF)
mark_as_advanced(GMX_POWERPC_INVSQRT)
option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
mark_as_advanced(GMX_FAHCORE)

include(gmxDetectAcceleration)
if(NOT DEFINED GMX_CPU_ACCELERATION)
    if(CMAKE_CROSSCOMPILING)
        set(GMX_SUGGESTED_CPU_ACCELERATION "None")
    else(CMAKE_CROSSCOMPILING)
        gmx_detect_acceleration(GMX_SUGGESTED_CPU_ACCELERATION)
    endif(CMAKE_CROSSCOMPILING)
endif(NOT DEFINED GMX_CPU_ACCELERATION)

set(GMX_CPU_ACCELERATION "@GMX_SUGGESTED_CPU_ACCELERATION@"
    CACHE STRING "Accelerated CPU kernels. Pick one of: None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, BlueGene, Power6, Fortran")

set(GMX_FFT_LIBRARY "fftw3" 
    CACHE STRING "FFT library choices: fftw3,mkl,fftpack[built-in]")
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

include(gmxGetCompilerVersion)
get_compiler_version()

# gcc 4.4.x is buggy and crashes when compiling some files with O3 and OpenMP on.
# Detect here whether applying a workaround is needed and will apply it later
# on the affected files.
include(gmxGCC44O3BugWorkaround)
gmx_check_gcc44_bug_workaround_needed(GMX_USE_GCC44_BUG_WORKAROUND)

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
  if (GMX_OPENMM)
    set (GMX_BINARY_SUFFIX "-openmm")
    set (GMX_LIBS_SUFFIX "_openmm")
  endif(GMX_OPENMM)
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
if(GMX_POWERPC_INVSQRT)
  set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_POWERPC_INVSQRT")
endif(GMX_POWERPC_INVSQRT)

########################################################################
#Process MPI settings
########################################################################
include(gmxManageMPI)

#######################################################################
# Check for options incompatible with OpenMM build                    #
#######################################################################
if(GMX_OPENMM)
    # we'll use the built-in fft to avoid unnecessary dependencies
    string(TOUPPER ${GMX_FFT_LIBRARY} GMX_FFT_LIBRARY)
    if(NOT ${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
        message(STATUS "No external FFT libraries needed for the OpenMM build, switching to fftpack!")
        set(GMX_FFT_LIBRARY "fftpack" CACHE STRING 
		"No external FFT libraries needed for the OpenMM build, switching to  fftpack!" FORCE)
    endif()
    if(GMX_MPI)
        message(FATAL_ERROR "The OpenMM build is not compatible with MPI!")
    endif(GMX_MPI)
    if(GMX_THREAD_MPI)
        message(STATUS "Thread-MPI not compatible with OpenMM, disabled!")
        set(GMX_THREAD_MPI OFF CACHE BOOL
		"Thread-MPI not compatible with OpenMM build, disabled!" FORCE)
    endif(GMX_THREAD_MPI)
    if(GMX_OPENMP)
        message(STATUS "OpenMP multithreading not compatible with OpenMM, disabled")
        set(GMX_OPENMP OFF CACHE BOOL
            "OpenMP multithreading not compatible with OpenMM, disabled!" FORCE)
    endif()
    if(GMX_SOFTWARE_INVSQRT)
        set(GMX_SOFTWARE_INVSQRT OFF CACHE STRING 
                "The OpenMM build does not need GROMACS software 1/sqrt!" FORCE)
    endif(GMX_SOFTWARE_INVSQRT)
    string(TOUPPER ${GMX_CPU_ACCELERATION} GMX_CPU_ACCELERATION)
    if(NOT GMX_CPU_ACCELERATION STREQUAL "NONE")
        message(STATUS "Switching off CPU-based acceleration, the OpenMM build does not support/need any!")
        set(GMX_CPU_ACCELERATION "None" CACHE STRING
            "Switching off CPU-based acceleration, the OpenMM build does not support/need any!" FORCE)
    endif()
    if(GMX_FAHCORE)
        message(FATAL_ERROR "The OpenMM build does not support FAH build!")
    endif(GMX_FAHCORE)
    if(GMX_DOUBLE)
        message(FATAL_ERROR  "The OpenMM-build does not support double precision calculations!")
    endif()
    # mark as advanced the unused variables
    mark_as_advanced(FORCE GMX_CPU_ACCELERATION GMX_MPI GMX_FFT_LIBRARY 
        GMX_QMMM_PROGRAM GMX_THREAD_MPI GMX_DOUBLE)
else(GMX_OPENMM)
     mark_as_advanced(CLEAR GMX_CPU_ACCELERATION GMX_MPI GMX_FFT_LIBRARY 
        GMX_QMMM_PROGRAM GMX_THREAD_MPI GMX_DOUBLE)
endif(GMX_OPENMM)


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
check_function_exists(lstat             HAVE_LSTAT)
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




########################################################################
# Find external packages                                               #
########################################################################
if(UNIX)
    if(GMX_PREFER_STATIC_LIBS)
        # On Linux .a is the static library suffix, on Mac OS X .lib can also
        # be used, so we'll add both to the preference list.
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
        if(SHARED_LIBS_DEFAULT)
            if(BUILD_SHARED_LIBS) #Warn the user about the combination. But don't overwrite the request.
                message(WARNING "Static libraries requested, and shared Gromacs libraries requested.")
            elseif(NOT DEFINED BUILD_SHARED_LIBS) #Change default to OFF. Don't warn if it's already off.
                message(WARNING "Static libraries requested, the GROMACS libraries will also be build static (BUILD_SHARED_LIBS=OFF)")
                set(SHARED_LIBS_DEFAULT OFF)
            endif()
        endif()
    endif()
endif()
option(BUILD_SHARED_LIBS "Enable shared libraries (can be problematic with MPI, Windows)" ${SHARED_LIBS_DEFAULT})

option(GMX_GSL "Add support for gsl" OFF)
if (GMX_GSL)
  find_package(gsl)
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
	if(X11_FOUND)
    	set(HAVE_X11 1)
	endif(X11_FOUND)
endif(GMX_X11)

include(ThreadMPI)
set(THREAD_MPI_LIB thread_mpi)
if(GMX_THREAD_MPI)
    tmpi_get_source_list(THREAD_MPI_SRC)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_THREAD_MPI")
    set(GMX_MPI 1)
else(GMX_THREAD_MPI)
    tmpi_get_source_list(THREAD_MPI_SRC NOMPI)
endif(GMX_THREAD_MPI)

if(GMX_OPENMM)
    set(CUDA_BUILD_EMULATION OFF)
    find_package(CUDA 3.1 REQUIRED)
    add_definitions(-DGMX_OPENMM)
    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")    
        set(CUDA_VERBOSE_BUILD ON)
    endif()
    find_package(OpenMM) 
endif(GMX_OPENMM)


if(GMX_GPU)
    if(GMX_DOUBLE)
        message(WARNING "GPU acceleration is not available in double precision, disabled!")
        set(GMX_GPU OFF CACHE BOOL "Enable GPU acceleration" FORCE)
    endif()
endif(GMX_GPU)

if(GMX_GPU)

    # We support CUDA >=v3.2 on *nix, but <= v4.1 doesn't work with MSVC
    if(MSVC)
        find_package(CUDA 4.1)
    else()
        find_package(CUDA 3.2)
    endif()

    if (NOT EXISTS ${CUDA_TOOLKIT_ROOT_DIR})
        message(FATAL_ERROR "
    mdrun supports native GPU acceleration on NVIDIA hardware with compute
    capability >=2.0. This requires the NVIDIA CUDA library, which was not
    found; the location can be hinted by setting CUDA_TOOLKIT_ROOT_DIR.

    CPU or GPU acceleration can be selected at runtime, but if you are
    sure you can not make use of GPU acceleration, disable it by setting
    the CMake variable GMX_GPU=OFF.")
    endif()

    if(NOT GMX_OPENMP)
        message(WARNING "
    In order to use GPU acceleration efficiently, mdrun requires OpenMP multithreding.
    Without OpenMP only a single CPU core per GPU can be used which is suboptimal.
    Note that with MPI multiple processes can be forced to use a single GPU, but this
    typically inefficient.")
    endif()

    include(gmxManageNvccConfig)

    # Check whether we can use atomic operations needed for polling wait for GPU
    # (to avoid the cudaStreamSynchronize + ECC bug).
    # With thread-MPI testing atomics has already been carried out, but without
    # thread-MPI we need to invoke the atomics test independently.
    if (NOT GMX_THREAD_MPI)
        set(TEST_TMPI_ATOMICS_ONLY ON CACHE INTERNAL
            "Test only the atomic operations of thread-MPI.")
        include(ThreadMPI)
    endif()

    # we need this linker flag in case if we have ld >= 2.22 (typically with gcc 4.5+),
    # but it's too cumbersome to check the ld version and the flag should not hurt
    if(CMAKE_COMPILER_IS_GNUCC)
        set(GROMACS_LINKER_FLAGS "-Wl,--add-needed ${GROMACS_LINKER_FLAGS}")
    endif()

    # annoyingly enough, FindCUDA leaves a few variables behind as non-advanced
    mark_as_advanced(CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD)
endif()

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

# Process nonbonded accelerated kernels settings
string(TOUPPER ${GMX_CPU_ACCELERATION} ${GMX_CPU_ACCELERATION})
if(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    # nothing to do
elseif(${GMX_CPU_ACCELERATION} STREQUAL "SSE2")

    GMX_TEST_CFLAG(GNU_SSE2_CFLAG "-msse2" GROMACS_C_FLAGS)
    if(NOT GNU_SSE2_CFLAG)
        GMX_TEST_CFLAG(MSVC_SSE2_CFLAG "/arch:SSE2" GROMACS_C_FLAGS)
    endif(NOT GNU_SSE2_CFLAG)

    if (CMAKE_CXX_COMPILER_LOADED)
        GMX_TEST_CXXFLAG(GNU_SSE2_CXXFLAG "-msse2" GROMACS_CXX_FLAGS)
        if(NOT GNU_SSE2_CXXFLAG)
            GMX_TEST_CXXFLAG(MSVC_SSE2_CXXFLAG "/arch:SSE2" GROMACS_CXX_FLAGS)
        endif(NOT GNU_SSE2_CXXFLAG)
    endif()

    # We dont warn for lacking SSE2 flag support, since that is probably standard today.

    # Only test the include after we have tried to add the correct flag for SSE2 support
    check_include_file(emmintrin.h  HAVE_EMMINTRIN_H ${GROMACS_C_FLAGS})

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

    GMX_TEST_CFLAG(GNU_SSE4_CFLAG "-msse4.1" GROMACS_C_FLAGS)
    if (NOT GNU_SSE4_CFLAG)
        GMX_TEST_CFLAG(MSVC_SSE4_CFLAG "/arch:SSE4.1" GROMACS_C_FLAGS)
    endif(NOT GNU_SSE4_CFLAG)
    if (NOT GNU_SSE4_CFLAG AND NOT MSVC_SSE4_CFLAG)
        message(WARNING "No C SSE4.1 flag found. Consider a newer compiler, or disable SSE4.1 for slightly lower performance.")
        # Not surprising if we end up here! MSVC current does not support the SSE4.1 flag. However, it appears to accept SSE4.1
        # intrinsics when SSE2 support is enabled, so we try that instead.
        GMX_TEST_CFLAG(MSVC_SSE2_CFLAG "/arch:SSE2" GROMACS_C_FLAGS)
    endif(NOT GNU_SSE4_CFLAG AND NOT MSVC_SSE4_CFLAG)

    if (CMAKE_CXX_COMPILER_LOADED)
        GMX_TEST_CXXFLAG(GNU_SSE4_CXXFLAG "-msse4.1" GROMACS_CXX_FLAG)
        if (NOT GNU_SSE4_CXXFLAG)
            GMX_TEST_CXXFLAG(MSVC_SSE4_CXXFLAG "/arch:SSE4.1" GROMACS_CXX_FLAGS)
        endif(NOT GNU_SSE4_CXXFLAG)
        if (NOT GNU_SSE4_CXXFLAG AND NOT MSVC_SSE4_CXXFLAG) 
            message(WARNING "No C++ SSE4.1 flag found. Consider a newer compiler, or disable SSE4.1 for slightly lower performance.")
            # Not surprising if we end up here! MSVC current does not support the SSE4.1 flag. However, it appears to accept SSE4.1
            # intrinsics when SSE2 support is enabled, so we try that instead.
            GMX_TEST_CXXFLAG(MSVC_SSE2_CXXFLAG "/arch:SSE2" GROMACS_CXX_FLAGS)
        endif(NOT GNU_SSE4_CXXFLAG AND NOT MSVC_SSE4_CXXFLAG)
    endif()

    # This must come after we have added the -msse4.1 flag on some platforms.
    check_include_file(smmintrin.h  HAVE_SMMINTRIN_H ${GROMACS_C_FLAGS})

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

elseif(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" OR ${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")

    # Set the AVX compiler flag for both these choices!

    GMX_TEST_CFLAG(GNU_AVX_CFLAG "-mavx" GROMACS_C_FLAGS)
    if (NOT GNU_AVX_CFLAG)
        GMX_TEST_CFLAG(MSVC_AVX_CFLAG "/arch:AVX" GROMACS_C_FLAGS)
    endif (NOT GNU_AVX_CFLAG)
    if (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)
        message(WARNING "No C AVX flag found. Consider a newer compiler, or disable AVX for much lower performance.")
    endif (NOT GNU_AVX_CFLAG AND NOT MSVC_AVX_CFLAG)

    if (CMAKE_CXX_COMPILER_LOADED)
        GMX_TEST_CXXFLAG(GNU_AVX_CXXFLAG "-mavx" GROMACS_CXX_FLAGS)
        if (NOT GNU_AVX_CXXFLAG)
            GMX_TEST_CXXFLAG(MSVC_AVX_CXXFLAG "/arch:AVX" GROMACS_CXX_FLAGS)
        endif (NOT GNU_AVX_CXXFLAG)
        if (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)
            message(WARNING "No C++ AVX flag found. Consider a newer compiler, or disable AVX for much lower performance.")
        endif (NOT GNU_AVX_CXXFLAG AND NOT MSVC_AVX_CXXFLAG)
    endif()

    # Set the FMA4 flags (MSVC doesn't require any)
    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" AND NOT MSVC)
        GMX_TEST_CFLAG(GNU_FMA_CFLAG "-mfma4" GROMACS_C_FLAGS)
        if (NOT GNU_FMA_CFLAG)
            message(WARNING "No C FMA4 flag found. Consider a newer compiler, or disable AVX_128_FMA for much lower performance.")
        endif(NOT GNU_FMA_CFLAG)
        if (CMAKE_CXX_COMPILER_LOADED)
            GMX_TEST_CXXFLAG(GNU_FMA_CXXFLAG "-mfma4" GROMACS_CXX_FLAGS)
            if (NOT GNU_FMA_CXXFLAG)
                message(WARNING "No C++ FMA flag found. Consider a newer compiler, or disable AVX_128_FMA for much lower performance.")
            endif (NOT GNU_FMA_CXXFLAG)
        endif()
    endif()

    # Only test the header after we have tried to add the flag for AVX support
    check_include_file(immintrin.h  HAVE_IMMINTRIN_H ${GROMACS_C_FLAGS})

    if(NOT HAVE_IMMINTRIN_H)
        message(FATAL_ERROR "Cannot find immintrin.h, which is required for AVX intrinsics support. Consider switching compiler.")
    endif(NOT HAVE_IMMINTRIN_H)

    if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")
        try_compile(TEST_AVX ${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/cmake/TestAVX.c"
            COMPILE_DEFINITIONS "${GROMACS_C_FLAGS}")
        if(NOT TEST_AVX)
            message(FATAL_ERROR "Cannot compile AVX intrinsics. Consider switching compiler.")
        endif()
    endif()

    # GCC requires x86intrin.h for FMA support. MSVC 2010 requires intrin.h for FMA support.
    check_include_file(x86intrin.h HAVE_X86INTRIN_H ${GROMACS_C_FLAGS})
    check_include_file(intrin.h HAVE_INTRIN_H ${GROMACS_C_FLAGS})

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

elseif(${GMX_CPU_ACCELERATION} STREQUAL "FORTRAN")

#    Fortran is temporarily disabled while we push in nbNxN kernels.
#    We need to fake it a bit here to avoid jenkins build errors!
#    add_definitions(-DGMX_FORTRAN)

elseif(${GMX_CPU_ACCELERATION} STREQUAL "BLUEGENE")
# GMX_CPU_ACCELERATION=BlueGene should be set in the Toolchain-BlueGene?-???.cmake file
    if (NOT ACCELERATION_QUIETLY)
      message(STATUS "Configuring for BlueGene")
    endif()
    set(GMX_BLUEGENE 1)
    if (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
        set(SHARED_LIBS_DEFAULT OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
    endif (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
    set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on BlueGene" FORCE)
    set(GMX_POWERPC_INVSQRT ON CACHE BOOL "Use hardware reciprocal square root on BlueGene" FORCE)
    set(GMX_X11 OFF CACHE BOOL "X11 not compatible with BlueGene, disabled!" FORCE)
    set(GMX_THREAD_MPI OFF CACHE BOOL "Thread-MPI not compatible with BlueGene, disabled!" FORCE)
    set(GMX_MPI ON CACHE BOOL "Use MPI on BlueGene" FORCE)
# Access to /etc/passwd is not available on the back end of BlueGene,
# despite being detected by CMake. This can cause linker warnings
# about harmless things in src/gmxlib/string2.h.
    set(HAVE_PWD_H OFF)
# The automatic testing for endianness does not work for the BlueGene cross-compiler
    set(GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER 1 CACHE INTERNAL "BlueGene has big endian FP byte order (by default)" FORCE)
    set(GMX_IEEE754_BIG_ENDIAN_WORD_ORDER 1 CACHE INTERNAL "BlueGene has big endian FP word order (by default)" FORCE)
elseif(${GMX_CPU_ACCELERATION} STREQUAL "POWER6")
    set(GMX_POWER6 1)
    set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on Power6" FORCE)
    set(GMX_POWERPC_INVSQRT ON CACHE BOOL "Use hardware reciprocal square root on Power6" FORCE)
else(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    MESSAGE(FATAL_ERROR "Unrecognized option for accelerated kernels: ${GMX_CPU_ACCELERATION}. Pick one of None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, Fortran, BlueGene, Power6")
endif(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
set(ACCELERATION_QUIETLY TRUE CACHE INTERNAL "")

if(GMX_FORTRAN OR GMX_POWER6)
    if (GMX_THREAD_MPI)
        message(FATAL_ERROR "FORTRAN/POWER6 is incompatible with thread-MPI and only provides a speed-up on certain IBM compilers. Disable FORTRAN (or threads if you really want to use FORTRAN kernels).")
    endif(GMX_THREAD_MPI)
    enable_language(Fortran)
    include(FortranCInterface)
    discover_fortran_mangling(prefix isupper suffix extra_under_score found)
    if(extra_under_score)
        set(extrasuffix "_")
    endif(extra_under_score)
    if(prefix)
      set(prefix "${prefix} ##")
    endif(prefix)
    if(suffix)
      set(suffix "## ${suffix}")
      if(extrasuffix)
	set(extrasuffix "${suffix}${extrasuffix}")
      endif(extrasuffix)
    else(suffix)
      if(extrasuffix)
	# Don't know if this is needed, but it can't hurt
	set(extrasuffix "## ${extrasuffix}")
      endif(extrasuffix)
    endif(suffix)

    if(isupper)
        set(F77_FUNCDEF   "${prefix} NAME ${suffix}")
        set(F77_FUNCDEF_  "${prefix} NAME ${extrasuffix}")
    else(isupper)
        set(F77_FUNCDEF   "${prefix} name ${suffix}")
        set(F77_FUNCDEF_  "${prefix} name ${extrasuffix}")
    endif(isupper)
else(GMX_FORTRAN OR GMX_POWER6)
        set(F77_FUNCDEF   "name ## _")
        set(F77_FUNCDEF_  "name ## _")
endif(GMX_FORTRAN OR GMX_POWER6)

# Process QM/MM Settings
string(TOUPPER ${GMX_QMMM_PROGRAM} ${GMX_QMMM_PROGRAM})
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

# Process FFT library settings - if not OpenMM build 
string(TOUPPER ${GMX_FFT_LIBRARY} ${GMX_FFT_LIBRARY})
set(PKG_FFT "")
set(PKG_FFT_LIBS "")
if(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    if(GMX_DOUBLE)
        set(FFTW fftw)
    else(GMX_DOUBLE)
        set(FFTW fftwf)
    endif(GMX_DOUBLE)
    find_package(FFTW COMPONENTS ${FFTW})
    string(TOUPPER "${FFTW}" FFTW)
    set(PKG_FFT "${${FFTW}_PKG}")
    include_directories(${${FFTW}_INCLUDE_DIRS})
    set(FFT_LIBRARIES ${${FFTW}_LIBRARIES})
    if(NOT ${FFTW}_FOUND)
      MESSAGE(FATAL_ERROR "Cannot find FFTW3 (with correct precision - libfftw3f for single precision GROMACS or libfftw3 for double precision GROMACS). Fix it, choose another FFT library, or use the Gromacs built-in fftpack (slower)!")
    endif(NOT ${FFTW}_FOUND)

    set(GMX_FFT_FFTW3 1)

    if (NOT ${GMX_CPU_ACCELERATION} STREQUAL "NONE" AND NOT ${FFTW}_HAVE_SIMD) 
      message(WARNING "The fftw library found is compiled without SIMD support, which makes it slow. Consider recompiling it or contact your admin")
    endif (NOT ${GMX_CPU_ACCELERATION} STREQUAL "NONE" AND NOT ${FFTW}_HAVE_SIMD) 

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
endif(GMX_FAHCORE)

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (NOT DEFINED GROMACS_C_FLAGS_SET)
    set(GROMACS_C_FLAGS_SET true CACHE INTERNAL "Whether to reset the C flags" 
        FORCE)
    set(CMAKE_C_FLAGS "${GROMACS_C_FLAGS} ${CMAKE_C_FLAGS}" CACHE STRING 
        "Flags used by the compiler during all build types" FORCE)
    if (CMAKE_CXX_COMPILER_LOADED)
        set(CMAKE_CXX_FLAGS "${GROMACS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}" CACHE STRING 
            "Flags used by the compiler during all build types" FORCE)
    endif()
    set(CMAKE_EXE_LINKER_FLAGS 
        "${GROMACS_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}" 
        CACHE STRING "Linker flags for creating executables" FORCE) 
    set(CMAKE_SHARED_LINKER_FLAGS 
        "${GROMACS_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}" 
        CACHE STRING "Linker flags for creating shared libraries" FORCE) 
endif (NOT DEFINED GROMACS_C_FLAGS_SET)

if(NOT GMX_OPENMP)
    #Unset all OpenMP flags in case OpenMP was disabled either by the user
    #or because it was only partially detected (e.g. only for C but not C++ compiler)
    unset(OpenMP_C_FLAGS CACHE) 
    unset(OpenMP_CXX_FLAGS CACHE)
    unset(OpenMP_LINKER_FLAGS CACHE)
    unset(OpenMP_SHARED_LINKER_FLAGS)
endif()
######################################
# Output compiler and CFLAGS used
######################################
execute_process(COMMAND ${CMAKE_C_COMPILER} --version RESULT_VARIABLE TMP_RESULT OUTPUT_VARIABLE CC_VERSION ERROR_VARIABLE CC_VERSION)
#try executing just the compiler command if that failed...
if(TMP_RESULT)
    execute_process(COMMAND ${CMAKE_C_COMPILER} RESULT_VARIABLE TMP_RESULT OUTPUT_VARIABLE CC_VERSION ERROR_VARIABLE CC_VERSION)
endif()
string(LENGTH "${CC_VERSION}" len)
if(len)
    string(REGEX MATCH "[^\n]+" CC_VERSION "${CC_VERSION}")
endif()

set(BUILD_COMPILER "@CMAKE_C_COMPILER@ @CMAKE_C_COMPILER_ID@ @CC_VERSION@")
set(BUILD_CFLAGS "@CMAKE_C_FLAGS@")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_CFLAGS "@BUILD_CFLAGS@ @CMAKE_C_FLAGS_DEBUG@")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_CFLAGS "@BUILD_CFLAGS@ @CMAKE_C_FLAGS_RELEASE@")
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(BUILD_CFLAGS "@BUILD_CFLAGS@ @CMAKE_C_FLAGS_MINSIZEREL@")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_CFLAGS "@BUILD_CFLAGS@ @CMAKE_C_FLAGS_RELWITHDEBINFO@")
endif()


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
IF(BUILD_TESTING)
	enable_testing()
	add_subdirectory(tests)
ENDIF()

