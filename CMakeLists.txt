cmake_minimum_required(VERSION 2.6)
# Keep CMake suitably quiet on Cygwin
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

project(Gromacs)
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
set(PROJECT_VERSION "4.5.6-dev")
set(CUSTOM_VERSION_STRING ""
    CACHE STRING "Custom version string (if empty, use hard-coded default)")
mark_as_advanced(CUSTOM_VERSION_STRING)
if (CUSTOM_VERSION_STRING)
    set(PROJECT_VERSION ${CUSTOM_VERSION_STRING})
endif (CUSTOM_VERSION_STRING)
set(SOVERSION 6)
# It is a bit irritating, but this has to be set separately for now!
SET(CPACK_PACKAGE_VERSION_MAJOR "4")
SET(CPACK_PACKAGE_VERSION_MINOR "5")
SET(CPACK_PACKAGE_VERSION_PATCH "6")


# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(CMAKE_INSTALL_PREFIX "/usr/local/gromacs" CACHE STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(C)

set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_VENDOR "gromacs.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gromacs - a toolkit for high-performance molecular simulation")
set(CPACK_SET_DESTDIR "ON")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0)
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/admin/InstallWelcome.txt")
# Its GPL/LGPL, so they do not have to agree to a license for mere usage, but some installers require this...
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/admin/InstallLicense.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/admin/InstallInfo.txt")

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
# Fix stupid flags on Windows
########################################################################
SET(SHARED_LIBS_DEFAULT ON) 
IF( MSVC )
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE) 
    SET(SHARED_LIBS_DEFAULT OFF)	
ENDIF()

IF( WIN32 AND NOT CYGWIN AND CMAKE_C_COMPILER_ID MATCHES "Intel" )
    STRING(REPLACE /GZ /RTC1 CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /GZ /RTC1 CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE) 
    SET(SHARED_LIBS_DEFAULT OFF)	
ENDIF()

set(GMX_EXTRA_LIBRARIES)



######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

include(gmxCFlags)
gmx_c_flags()

########################################################################
# User input options                                                   #
########################################################################
option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" OFF)
option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREADS    "Build a parallel (thread-based) version of GROMACS (cannot be combined with MPI yet)" ON)
option(GMX_SOFTWARE_INVSQRT "Use GROMACS software 1/sqrt" ON)
mark_as_advanced(GMX_SOFTWARE_INVSQRT)
option(GMX_POWERPC_INVSQRT "Use PowerPC hardware 1/sqrt" OFF)
mark_as_advanced(GMX_POWERPC_INVSQRT)
option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
mark_as_advanced(GMX_FAHCORE)
option(GMX_OPENMM "Accelerated execution on GPUs through the OpenMM library (rerun cmake after changing to see relevant options)" OFF)
set(GMX_ACCELERATION "auto" 
    CACHE STRING "Accelerated kernels. Pick one of: auto, none, SSE, BlueGene, Power6, ia64, altivec, fortran")

set(GMX_FFT_LIBRARY "fftw3" 
    CACHE STRING "FFT library choices: fftw3,fftw2,mkl,fftpack[built-in]")
option(GMX_DISABLE_FFTW_MEASURE 
       "Do not optimize FFTW setups (not needed with SSE)" OFF)
mark_as_advanced(GMX_DISABLE_FFTW_MEASURE)
set(GMX_QMMM_PROGRAM "none" 
    CACHE STRING "QM package choices: none,gaussian,mopac,gamess,orca")
option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
mark_as_advanced(GMX_BROKEN_CALLOC)
option(BUILD_SHARED_LIBS "Enable shared libraries (can be problematic with MPI, Windows)" ${SHARED_LIBS_DEFAULT})
option(GMX_MPI_IN_PLACE "Enable MPI_IN_PLACE for MPIs that have it defined" ON)
mark_as_advanced(GMX_MPI_IN_PLACE)
option(GMX_DLOPEN "Compile with dlopen, needed to read VMD supported file formats" ON)
mark_as_advanced(GMX_DLOPEN)


option(GMX_IA32_ASM "Add SSE assembly files for IA32" OFF)
mark_as_advanced(GMX_IA32_ASM)
option(GMX_X86_64_ASM "Add SSE assembly files for X86_64" OFF)
mark_as_advanced(GMX_X86_64_ASM)

option(USE_VERSION_H "Generate development version string/information" ON)
mark_as_advanced(USE_VERSION_H)

option(GMX_DEFAULT_SUFFIX "Use default suffixes for GROMACS binaries and libs (_d for double, _mpi for MPI; rerun cmake after changing to see relevant options)" ON)

if(UNIX AND NOT APPLE)
    option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer \".a\" static archives (NOTE: this is enabled only for UNIX (excluding APPLE) platforms but it might not always work!" OFF)
    mark_as_advanced(GMX_PREFER_STATIC_LIBS)
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
  if (GMX_OPENMM)
    set (GMX_BINARY_SUFFIX "-gpu")
    set (GMX_LIBS_SUFFIX "_gpu")
  endif(GMX_OPENMM)
  mark_as_advanced(FORCE GMX_BINARY_SUFFIX GMX_LIBS_SUFFIX)
  message(STATUS "Using default binary suffix: \"${GMX_BINARY_SUFFIX}\"")    
  message(STATUS "Using default library suffix: \"${GMX_LIBS_SUFFIX}\"") 
else(GMX_DEFAULT_SUFFIX)
  mark_as_advanced(CLEAR GMX_BINARY_SUFFIX GMX_LIBS_SUFFIX)
  message(STATUS "Using manually set binary suffix: \"${GMX_BINARY_SUFFIX}\"")    
  message(STATUS "Using manually set library suffix: \"${GMX_LIBS_SUFFIX}\"")
endif(GMX_DEFAULT_SUFFIX)

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
include(CheckCSourceCompiles) # for gmxTestMPI_IN_PLACE
if(GMX_MPI)
    if(GMX_THREADS)
        #message(FATAL_ERROR "Thread-based parallelization conflicts with MPI.")
        set(GMX_THREADS OFF CACHE BOOL 
            "Thread-based parallelization conflicts with MPI." FORCE)
    endif(GMX_THREADS)
    find_package(MPI)
    if(MPI_FOUND)
        set(GROMACS_C_FLAGS ${GROMACS_FLAGS} ${MPI_COMPILE_FLAGS})
	set(GROMACS_LINKER_FLAGS ${GROMACS_LINKER_FLAGS} ${MPI_LINK_FLAGS})
        
        include_directories(${MPI_INCLUDE_PATH})
        list(APPEND GMX_EXTRA_LIBRARIES ${MPI_LIBRARIES})
        if(GMX_FAHCORE)
            add_definitions( -DMPI ) #for FAHCORE
        endif(GMX_FAHCORE)
        include(gmxTestMPI_IN_PLACE)
        if (GMX_MPI_IN_PLACE)
            gmx_test_mpi_in_place(MPI_IN_PLACE_EXISTS)
        endif (GMX_MPI_IN_PLACE)

	# test for unsuitable versions of OpenMPI
	exec_program(${MPIEXEC}
	  ARGS --version
	  OUTPUT_VARIABLE MPI_TYPE
	  RETURN_VALUE MPI_EXEC_RETURN)
	if(MPI_EXEC_RETURN EQUAL 0)
	  if(MPI_TYPE MATCHES "Open MPI|OpenRTE")
	    string(REGEX MATCH "[0-9]+\\.[0-9]*\\.?[0-9]*" MPI_VERSION ${MPI_TYPE})
	    if(MPI_VERSION VERSION_LESS "1.4.1")
	      MESSAGE(WARNING "
            There are known problems with OpenMPI version < 1.4.1.
            Please consider updating your OpenMPI.")
	    endif(MPI_VERSION VERSION_LESS "1.4.1")
	    unset(MPI_VERSION)
	  else(MPI_TYPE MATCHES "Open MPI|OpenRTE")
	    # This is not OpenMPI, so give the old generic warning message
	    MESSAGE(WARNING "
            There are known problems with some MPI implementations:
                     MVAPICH2 version <= 1.4.1
            Please consider updating your MPI if applicable.")
	  endif(MPI_TYPE MATCHES "Open MPI|OpenRTE")
	  unset(MPI_TYPE)
	endif(MPI_EXEC_RETURN EQUAL 0)
    else(MPI_FOUND)
        message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
    endif(MPI_FOUND)
    include(gmxTestCatamount)
    gmx_test_catamount(GMX_CRAY_XT3)
    if(GMX_CRAY_XT3)
        set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_CRAY_XT3")
        set(GMX_NO_SYSTEM 1)
        set(GMX_NO_NICE 1)
    endif(GMX_CRAY_XT3)
    set(GMX_LIB_MPI 1)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_LIB_MPI")
endif(GMX_MPI)


#######################################################################
# Check for options incompatible with OpenMM build                    #
#######################################################################
if(GMX_OPENMM)
    cmake_minimum_required(VERSION 2.6.4)
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
    if(GMX_THREADS)
        message(STATUS "Threads are  not compatible with OpenMM build, disabled!")
        set(GMX_THREADS OFF CACHE BOOL 
		"Threads are not compatible with OpenMM build, disabled!" FORCE)
    endif(GMX_THREADS)
    if(GMX_SOFTWARE_INVSQRT)
        set(GMX_SOFTWARE_INVSQRT OFF CACHE STRING 
                "The OpenMM build does not need GROMACS software 1/sqrt!" FORCE)
    endif(GMX_SOFTWARE_INVSQRT)
    string(TOUPPER ${GMX_ACCELERATION} GMX_ACCELERATION)
    if(NOT GMX_ACCELERATION STREQUAL "NONE")
        message(STATUS "Switching off CPU-based acceleration, the OpenMM build does not support/need any!")	
        set(GMX_ACCELERATION "none" CACHE STRING 
		"Switching off CPU-based acceleration, the OpenMM build does not support/need any!" FORCE)
    endif()
    if(GMX_FAHCORE)
        message(FATAL_ERROR "The OpenMM build does not support FAH build!")
    endif(GMX_FAHCORE)
    if(GMX_DOUBLE)
        message(FATAL_ERROR  "The OpenMM-build does not support double precision calculations!")
    endif()
    # mark as advanced the unused variables
    mark_as_advanced(FORCE GMX_ACCELERATION GMX_MPI GMX_FFT_LIBRARY 
        GMX_QMMM_PROGRAM GMX_THREADS GMX_DOUBLE)
else(GMX_OPENMM)
     mark_as_advanced(CLEAR GMX_ACCELERATION GMX_MPI GMX_FFT_LIBRARY 
        GMX_QMMM_PROGRAM GMX_THREADS GMX_DOUBLE)   
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
# SSE support
check_include_files(xmmintrin.h  HAVE_XMMINTRIN_H)
check_include_files(emmintrin.h  HAVE_EMMINTRIN_H)
check_include_files(pmmintrin.h  HAVE_PMMINTRIN_H)
check_include_files(smmintrin.h  HAVE_SMMINTRIN_H)
check_include_files(io.h  		 HAVE_IO_H)


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
if(UNIX AND NOT APPLE)
    if(GMX_PREFER_STATIC_LIBS)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
        if(BUILD_SHARED_LIBS)
            message(WARNING "Static libraries requested, the GROMACS libraries will also be build static (BUILD_SHARED_LIBS=OFF)")
            set(BUILD_SHARED_LIBS OFF CACHE BOOL "Enable shared libraries (can be problematic with MPI, Windows)" FORCE)
        endif()
    endif()
endif()

option(GMX_XML "Use libxml2 to parse xml files" ON)
if (GMX_XML)
  find_package(LibXml2)
  set(PKG_XML "")
  if(LIBXML2_FOUND)
    include_directories(${LIBXML2_INCLUDE_DIR})
    set(PKG_XML libxml-2.0)
    set(XML_LIBRARIES ${LIBXML2_LIBRARIES})
    set(HAVE_LIBXML2 1)
  endif(LIBXML2_FOUND)
endif(GMX_XML)

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

if(GMX_THREADS)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_THREADS")
    include(ThreadMPI)
    set(THREAD_MPI_LIB thread_mpi)
    set(GMX_MPI 1)
    string(TOUPPER ${GMX_FFT_LIBRARY} ${GMX_FFT_LIBRARY})
    if(${GMX_FFT_LIBRARY} STREQUAL "FFTW2")
        message(FATAL_ERROR "FFTW2 can't be used with threads. Try fftw3 or mkl.")
    endif()
endif(GMX_THREADS)

if(GMX_OPENMM)
    set(CUDA_BUILD_EMULATION OFF)
    find_package(CUDA 3.1 REQUIRED)
    add_definitions(-DGMX_OPENMM)
    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")    
        set(CUDA_VERBOSE_BUILD ON)
    endif()
    find_package(OpenMM) 
endif(GMX_OPENMM)

if(APPLE)
   find_library(ACCELERATE_FRAMEWORK Accelerate)
   list(APPEND GMX_EXTRA_LIBRARIES ${ACCELERATE_FRAMEWORK})
endif(APPLE)

if(CYGWIN)
    set(GMX_CYGWIN 1)
endif(CYGWIN)

# only bother with finding git and using version.h if the source is a git repo
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    if(USE_VERSION_H)
        # We need at least git v1.5.3 be able to parse git's date output. If not
        # fund or the version is too small, we can't generate version information.
        find_package(Git 1.5.3)
        # this should at some point become VERSION_LESS
        if(NOT Git_FOUND OR Git_VERSION STRLESS "1.5.3")
            message("No compatible git version found, won't be able to generate proper development version information.")
            set(USE_VERSION_H OFF)
        endif()
    endif()
else()
    set(USE_VERSION_H OFF)
endif()

if (GMX_DLOPEN)
    list(APPEND GMX_EXTRA_LIBRARIES ${CMAKE_DL_LIBS})
    set(PKG_DL_LIBS "-l${CMAKE_DL_LIBS}")
else(GMX_DLOPEN)
    set(PKG_DL_LIBS)
endif (GMX_DLOPEN)

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
include_directories(${CMAKE_SOURCE_DIR}/include)

include(gmxCheckBuildUserTime)
gmx_check_build_user_time(BUILD_TIME BUILD_USER BUILD_MACHINE)

include(gmxTestFloatFormat)
gmx_test_float_format(GMX_FLOAT_FORMAT_IEEE754 
                      GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER
                      GMX_IEEE754_BIG_ENDIAN_WORD_ORDER)

include(gmxTestLargeFiles)
gmx_test_large_files(GMX_LARGEFILES)

include(gmxTestSignal)
gmx_test_retsigtype(RETSIGTYPE)
gmx_test_sigusr1(HAVE_SIGUSR1)

include(gmxTestInline)
gmx_test_inline(INLINE_KEYWORD)

include(gmxTestRestrict)
gmx_test_inline(RESTRICT_KEYWORD)

include(gmxTestPipes)
gmx_test_pipes(HAVE_PIPES)

include(gmxTestIsfinite)
gmx_test_isfinite(HAVE_ISFINITE)
gmx_test__isfinite(HAVE__ISFINITE)
gmx_test__finite(HAVE__FINITE)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)
gmx_test_inline_asm_msvc_x86(GMX_X86_MSVC_INLINE_ASM)

# turn on SSE if supported with reasonable defaults.
if (${GMX_ACCELERATION} STREQUAL "auto" AND NOT GMX_OPENMM)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(i.86|x86|x64|x86_64|AMD64|amd64)" OR CYGWIN)

    set(GMX_ACCELERATION "SSE" CACHE STRING "Accelerated kernels. Pick one of: auto, none, SSE, BlueGene, Power6, ia64, altivec, fortran" FORCE)
    
    # Determine the assembler/compiler to use
  else()
    set(GMX_ACCELERATION "none" CACHE STRING "Accelerated kernels. Pick one of: auto, none, SSE, BlueGene, Power6, ia64, altivec, fortran" FORCE)
  endif()
endif (${GMX_ACCELERATION} STREQUAL "auto" AND NOT GMX_OPENMM)

include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
if(NOT GMX_SYSTEM_XDR)
    set(GMX_INTERNAL_XDR 1)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_INTERNAL_XDR")
endif(NOT GMX_SYSTEM_XDR)

# Process nonbonded accelerated kernels settings
string(TOUPPER ${GMX_ACCELERATION} ${GMX_ACCELERATION})
if(${GMX_ACCELERATION} STREQUAL "NONE")
  # nothing to do
elseif(${GMX_ACCELERATION} STREQUAL "SSE")
    
    if (GMX_64_BIT)
      set(GMX_X86_64_ASM ON CACHE BOOL "Add SSE assembly files for x86_64" FORCE)
    else (GMX_64_BIT)
      set(GMX_IA32_ASM ON CACHE BOOL "Add SSE assembly files for i386" FORCE)
    endif (GMX_64_BIT)

    if( WIN32 AND NOT CYGWIN )
      option(GMX_ASM_USEASM_NASM "Use Nasm for assembly instead of compiler (necessary on windows)" ON)
    else()
      option(GMX_ASM_USEASM_NASM "Use Nasm for assembly instead of compiler (necessary on windows)" OFF)
    endif()
	mark_as_advanced(GMX_ASM_USEASM_NASM)

    if (NOT GMX_64_BIT)
        # for 32-bit compiles, we might need to turn on sse 
        CHECK_C_COMPILER_FLAG("-msse2" XFLAGS_SSE)
        if (XFLAGS_SSE)
            set(GROMACS_C_FLAGS "-msse2 ${GROMACS_C_FLAGS}")
        endif (XFLAGS_SSE)
        CHECK_CXX_COMPILER_FLAG("-msse2" XXFLAGS_SSE)
        if (XXFLAGS_SSE)
            set(GROMACS_CXX_FLAGS "-msse2 ${GROMACS_CXX_FLAGS}")
        endif (XXFLAGS_SSE)
    endif (NOT GMX_64_BIT)

    if(GMX_IA32_ASM)
      if(GMX_DOUBLE)
        set(GMX_IA32_SSE2 1)
      else()
       set(GMX_IA32_SSE 1)
      endif()
    elseif(GMX_X86_64_ASM)
      if(GMX_DOUBLE)
       set(GMX_X86_64_SSE2 1)
      else()
        set(GMX_X86_64_SSE 1)
      endif()
    endif()

#    if(HAVE_PMMINTRIN_H)
#        set(GMX_SSE3 1)
#    endif(HAVE_PMMINTRIN_H)
#    if(HAVE_SMMINTRIN_H)
#        set(GMX_SSE4_1 1)
#    endif(HAVE_SMMINTRIN_H)

elseif(${GMX_ACCELERATION} STREQUAL "FORTRAN")
    set(GMX_FORTRAN 1)
    #these are switch on by default sometimes
    set(GMX_IA32_ASM 0)
    set(GMX_GMX_X86_64_ASM 0)
elseif(${GMX_ACCELERATION} STREQUAL "BLUEGENE")
# GMX_ACCELERATION=BlueGene should be set in the Toolchain-BlueGene?-???.cmake file
    message(STATUS "Configuring for BlueGene")
    set(GMX_BLUEGENE 1)
    if (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
        set(SHARED_LIBS_DEFAULT OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
    endif (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
    set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on BlueGene" FORCE)
    set(GMX_POWERPC_INVSQRT ON CACHE BOOL "Use hardware reciprocal square root on BlueGene" FORCE)
    set(GMX_X11 OFF CACHE BOOL "X11 not compatible with BlueGene, disabled!" FORCE)
    set(GMX_THREADS OFF CACHE BOOL "Threads not compatible with BlueGene, disabled!" FORCE)
    set(GMX_MPI ON CACHE BOOL "Use MPI on BlueGene" FORCE)
# Access to /etc/passwd is not available on the back end of BlueGene,
# despite being detected by CMake. This can cause linker warnings
# about harmless things in src/gmxlib/string2.h.
    set(HAVE_PWD_H OFF)
# The automatic testing for endianness does not work for the BlueGene cross-compiler
    set(GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER 1 CACHE INTERNAL "BlueGene has big endian FP byte order (by default)" FORCE)
    set(GMX_IEEE754_BIG_ENDIAN_WORD_ORDER 1 CACHE INTERNAL "BlueGene has big endian FP word order (by default)" FORCE)
elseif(${GMX_ACCELERATION} STREQUAL "POWER6")
    set(GMX_POWER6 1)
    set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on Power6" FORCE)
    set(GMX_POWERPC_INVSQRT ON CACHE BOOL "Use hardware reciprocal square root on Power6" FORCE)
elseif(${GMX_ACCELERATION} STREQUAL "IA64")
    set(GMX_IA64_ASM 1)
    set(DISABLE_WATERWATER_NLIST 1)
    set(DISABLE_WATER_NLIST 1)
elseif(${GMX_ACCELERATION} STREQUAL "ALTIVEC")
    check_include_files(altivec.h HAVE_ALTIVEC_H)
    if(HAVE_ALTIVEC_H)
        set(GMX_PPC_ALTIVEC 1)
    endif(HAVE_ALTIVEC_H)
else(${GMX_ACCELERATION} STREQUAL "NONE")
    MESSAGE(FATAL_ERROR "Unrecognized option for accelerated kernels: ${GMX_ACCELERATION}. Pick one of auto, none, SSE, Fortran, BlueGene, Power6, ia64, altivec")
endif(${GMX_ACCELERATION} STREQUAL "NONE")

if(GMX_FORTRAN OR GMX_POWER6)
    if (GMX_THREADS)
        message(FATAL_ERROR "FORTRAN/POWER6 is incompatible with threads and only provides a speed-up on certain IBM compilers. Disable FORTRAN (or threads if you really want to use FORTRAN kernels).")
    endif(GMX_THREADS)
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
#    MESSAGE(STATUS "Using external FFT library - fftw3")
    if(GMX_DOUBLE)
        find_package(FFTW3 REQUIRED)
		include_directories(${FFTW3_INCLUDE_DIR})
        set(FFT_LIBRARIES ${FFTW3_LIBRARIES})
        set(PKG_FFT "fftw3")
    else(GMX_DOUBLE)
        find_package(FFTW3F REQUIRED)
        include_directories(${FFTW3F_INCLUDE_DIR})
        set(FFT_LIBRARIES ${FFTW3F_LIBRARIES})
        set(PKG_FFT "fftw3f")
    endif(GMX_DOUBLE)

    if(NOT FFTW3_FOUND AND NOT FFTW3F_FOUND)
        MESSAGE(FATAL_ERROR "Cannot find fftw3 (with correct precision). Fix it, choose another FFT library, or use the Gromacs built-in fftpack (slower)!")
    endif(NOT FFTW3_FOUND AND NOT FFTW3F_FOUND)

    set(GMX_FFT_FFTW3 1)

elseif(${GMX_FFT_LIBRARY} STREQUAL "FFTW2")
#    MESSAGE(STATUS "Using external FFT library - fftw2")
    if(GMX_DOUBLE)
        find_package(FFTW2 REQUIRED)
    else(GMX_DOUBLE)
        find_package(FFTW2F REQUIRED)
    endif(GMX_DOUBLE)

    if(NOT FFTW2_FOUND)
      	   MESSAGE(FATAL_ERROR "Cannot find fftw2 (with correct precision). Fix it, choose another FFT library, or use the Gromacs built-in fftpack (slower)!")
    endif(NOT FFTW2_FOUND)
    include_directories(${FFTW2_INCLUDE_DIR})
    set(FFT_LIBRARIES ${FFTW2_LIBRARIES})
    set(PKG_FFT_LIBS ${FFTW2_LIBRARIES})

    if("${FFTW2_LIBRARIES}" MATCHES "dfftw")
        set(FFTW2_NAME_DFFTW 1)
    elseif("${FFTW2_LIBRARIES}" MATCHES "sfftw")
        set(FFTW2_NAME_SFFTW 1)
    else("${FFTW2_LIBRARIES}" MATCHES "dfftw")
        set(FFTW2_NAME_FFTW 1) 
    endif("${FFTW2_LIBRARIES}" MATCHES "dfftw")

    set(GMX_FFT_FFTW2 1)
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
    MESSAGE(FATAL_ERROR "Invalid FFT library setting: ${GMX_FFT_LIBRARY}. Choose one of: fftw3, fftw2, mkl, acml, fftpack")
endif(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")

# MKL has BLAS/LAPACK routines
if(HAVE_MKL OR ACCELERATE_FRAMEWORK)
  set(GMX_EXTERNAL_BLAS TRUE CACHE BOOL "Use external BLAS instead of built-in")
  set(GMX_EXTERNAL_LAPACK TRUE CACHE BOOL "Use external LAPACK instead of built-in")
else(HAVE_MKL OR ACCELERATE_FRAMEWORK)
  set(GMX_EXTERNAL_BLAS FALSE CACHE BOOL "Use external BLAS instead of built-in") 
  set(GMX_EXTERNAL_LAPACK FALSE CACHE BOOL "Use external LAPACK instead of built-in") 
  if(GMX_EXTERNAL_BLAS)
    #FindBLAS needs Fortran
    enable_language(Fortran)
    if (GMX_BLAS_USER)
        list(APPEND GMX_EXTRA_LIBRARIES ${GMX_BLAS_USER})
    else(GMX_BLAS_USER)
        find_package(BLAS REQUIRED)
        list(APPEND GMX_EXTRA_LIBRARIES ${BLAS_LIBRARIES})
    endif(GMX_BLAS_USER)
  endif(GMX_EXTERNAL_BLAS)
  if(GMX_EXTERNAL_LAPACK)
    #FindLAPACK needs Fortran
    enable_language(Fortran)
	if (GMX_LAPACK_USER)
        list(APPEND GMX_EXTRA_LIBRARIES ${GMX_LAPACK_USER})
	else(GMX_LAPACK_USER)
        find_package(LAPACK REQUIRED)
        list(APPEND GMX_EXTRA_LIBRARIES ${LAPACK_LIBRARIES})
    endif(GMX_LAPACK_USER)
  endif(GMX_EXTERNAL_LAPACK)
endif(HAVE_MKL OR ACCELERATE_FRAMEWORK)
mark_as_advanced(GMX_EXTERNAL_LAPACK)
mark_as_advanced(GMX_EXTERNAL_BLAS)

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
    set(CMAKE_CXX_FLAGS "${GROMACS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}" CACHE STRING 
        "Flags used by the compiler during all build types" FORCE)
    set(CMAKE_EXE_LINKER_FLAGS 
        "${GROMACS_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}" 
        CACHE STRING "Linker flags" FORCE) 
endif (NOT DEFINED GROMACS_C_FLAGS_SET)

########################################################################
# Specify install locations and which subdirectories to process        #
########################################################################
if ( DEFINED GMXLIB )
    set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${GMXLIB}")
else()
    set(LIB_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/lib)
endif()
set(BIN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/bin)
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/gromacs)
set(MAN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/share/man)
set(INCL_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)

set(GMXLIBDIR        ${DATA_INSTALL_DIR}/top)

##################################################################
# Shared library settings - Darwin uses INSTALL_NAME_DIR instead!
##################################################################
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "${LIB_INSTALL_DIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

add_subdirectory(share)
add_subdirectory(include)
add_subdirectory(man)
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

