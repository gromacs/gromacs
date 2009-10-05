cmake_minimum_required(VERSION 2.6)

project(Gromacs)
set(PROJECT_VERSION "4.0.99_development_20090605"
    CACHE STRING "Gromacs version string")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

########################################################################
# User input options                                                   #
########################################################################
option(GMX_DOUBLE "Use double precision" OFF)
option(GMX_MPI    "Build a parallel (message-passing) version of GROMACS" OFF)
option(GMX_THREADS    "Build a parallel (threaded-based) version of GROMACS: (doesn't work yet!)" OFF)
option(GMX_SOFTWARE_INVSQRT "Use GROMACS software 1/sqrt" ON)
option(GMX_FAHCORE "Build a library with mdrun functionality" OFF)
set(GMX_ACCELERATION "none" 
    CACHE STRING "Accelerated kernels. Pick one of: SSE, BlueGene, Power6, ia64, altivec")

set(GMX_FFT_LIBRARY "fftw3" 
    CACHE STRING "FFT library choices: fftw3,fftw2,mkl,fftpack[built-in]")
option(GMX_DISABLE_FFTW_MEASURE 
       "Do not optimize FFTW setups (not needed with SSE)" OFF)
set(GMX_QMMM_PROGRAM "none" 
    CACHE STRING "QM package choices: none,gaussian,mopac,gamess")
option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
option(BUILD_SHARED_LIBS "Enable shared libraries (can be problematic with MPI)" OFF)
option(GMX_MPI_IN_PLACE "Enable MPI_IN_PLACE for MPIs that have it defined" ON)
mark_as_advanced(GMX_MPI_IN_PLACE)

enable_language(C)

########################################################################
# Fix stupid flags on MSVC
########################################################################
IF(CMAKE_GENERATOR MATCHES "Visual Studio")
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)   
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio")

set(GMX_EXTRA_LIBRARIES)

########################################################################
#Process MPI settings
########################################################################
if(GMX_MPI)
    if(GMX_THREADS)
        message(FATAL_ERROR "Thread-based parallelization conflicts with MPI.")
    endif(GMX_THREADS)
    find_package(MPI)
    if(MPI_FOUND)
        set(GROMACS_C_FLAGS ${GROMACS_FLAGS} ${MPI_COMPILE_FLAGS})
	set(GROMACS_LINKER_FLAGS ${GROMACS_LINKER_FLAGS} ${MPI_LINK_FLAGS})
        
        include_directories(${MPI_INCLUDE_PATH})
        list(APPEND GMX_EXTRA_LIBRARIES ${MPI_LIBRARIES})
        add_definitions( -DMPI ) #for FAHCORE
        include(gmxTestMPI_IN_PLACE)
        if (GMX_MPI_IN_PLACE)
            gmx_test_mpi_in_place(MPI_IN_PLACE_EXISTS)
        endif (GMX_MPI_IN_PLACE)
    else(MPI_FOUND)
        message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
    endif(MPI_FOUND)
    include(gmxTestCatamount)
    gmx_test_catamount(GMX_CRAY_XT3)
    if(GMX_CRAY_XT3)
        set(GMX_NO_SYSTEM 1)
        set(GMX_NO_NICE 1)
    endif(GMX_CRAY_XT3)
    set(GMX_LIB_MPI 1)
endif(GMX_MPI)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
check_include_files(string.h     HAVE_STRING_H)
check_include_files(limits.h     HAVE_LIMITS_H)
check_include_files(memory.h     HAVE_MEMORY_H)
check_include_files(unistd.h	 HAVE_UNISTD_H)
check_include_files(pwd.h        HAVE_PWD_H)
check_include_files(stdint.h	 HAVE_STDINT_H)
check_include_files(stdlib.h	 HAVE_STDLIB_H)
check_include_files(pthread.h    HAVE_PTHREAD_H)
check_include_files(inttypes.h   HAVE_INTTYPES_H)
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

include(CheckFunctionExists)
check_function_exists(strcasecmp        HAVE_STRCASECMP)
check_function_exists(strdup            HAVE_STRDUP)
check_function_exists(vprintf           HAVE_VPRINTF)
check_function_exists(memcmp            HAVE_MEMCMP)
check_function_exists(gettimeofday      HAVE_GETTIMEOFDAY)
check_function_exists(isnan             HAVE_ISNAN)
check_function_exists(_isnan            HAVE__ISNAN)
# After large files? check_function_exists(fseeko     HAVE_FSEEKO)

include(CheckLibraryExists)
check_library_exists(m sqrt "" HAVE_LIBM)

include(CheckTypeSize)
check_type_size("bool"          SIZEOF_BOOL) # will also set HAVE_BOOL
check_type_size("int"           SIZEOF_INT) 
check_type_size("long int"      SIZEOF_LONG_INT) 
check_type_size("long long int" SIZEOF_LONG_LONG_INT) 
check_type_size("off_t"         SIZEOF_OFF_T)
check_type_size("void *"        SIZEOF_VOIDP)

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

find_package(LibXml2)
if(LIBXML2_FOUND)
    include_directories(${LIBXML2_INCLUDE_DIR})
    list(APPEND GMX_EXTRA_LIBRARIES ${LIBXML2_LIBRARIES})
    set(HAVE_LIBXML2 1)
endif(LIBXML2_FOUND)

find_package(X11)
# X11 includes/libraries are only set in the ngmx subdirectory!
if(X11_FOUND)
    set(HAVE_X11 1)
endif(X11_FOUND)

if(GMX_THREADS)
    include(ThreadMPI)
    set(THREAD_MPI_LIB thread_mpi)
    set(GMX_MPI 1)
    string(TOUPPER ${GMX_FFT_LIBRARY} ${GMX_FFT_LIBRARY})
    if(${GMX_FFT_LIBRARY} STREQUAL "FFTW2")
        message(FATAL_ERROR "FFTW2 can't be used with threads. Try fftw3 or mkl.")
    endif()
endif(GMX_THREADS)


if(APPLE)
   find_library(ACCELERATE_FRAMEWORK Accelerate)
   list(APPEND GMX_EXTRA_LIBRARIES ${ACCELERATE_FRAMEWORK})
endif(APPLE)

########################################################################
# Our own GROMACS tests
########################################################################

include(gmxCFlags)
gmx_c_flags()

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

include(gmxTestPipes)
gmx_test_pipes(HAVE_PIPES)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)
gmx_test_inline_asm_msvc_x86(GMX_X86_MSVC_INLINE_ASM)

include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
if(NOT GMX_SYSTEM_XDR)
    set(GMX_INTERNAL_XDR 1)
endif(NOT GMX_SYSTEM_XDR)

if(GMX_FORTRAN OR GMX_POWER6)
    enable_language(Fortran)
    include(FortranCInterface)
    discover_fortran_mangling(prefix isupper suffix extra_under_score found)
    if(extra_under_score)
        set(extrasuffix "_")
    endif(extra_under_score)

    if(isupper)
        set(F77_FUNCDEF   "${prefix} ## NAME ## ${suffix}")
        set(F77_FUNCDEF_  "${prefix} ## NAME ## ${suffix}${extrasuffix}")
    else(isupper)
        set(F77_FUNCDEF   "${prefix} ## name ## ${suffix}")
        set(F77_FUNCDEF_  "${prefix} ## name ## ${suffix}${extrasuffix}")
    endif(isupper)
else(GMX_FORTRAN OR GMX_POWER6)
        set(F77_FUNCDEF   "name ## _")
        set(F77_FUNCDEF_  "name ## _")
endif(GMX_FORTRAN OR GMX_POWER6)


# Process nonbonded accelerated kernels settings
string(TOUPPER ${GMX_ACCELERATION} ${GMX_ACCELERATION})
if(${GMX_ACCELERATION} STREQUAL "NONE")
  # nothing to do
elseif(${GMX_ACCELERATION} STREQUAL "SSE")
    if(HAVE_XMMINTRIN_H)
        set(GMX_SSE 1)
    endif(HAVE_XMMINTRIN_H)
    if(HAVE_EMMINTRIN_H)
        set(GMX_SSE2 1)
    endif(HAVE_EMMINTRIN_H)
#    if(HAVE_PMMINTRIN_H)
#        set(GMX_SSE3 1)
#    endif(HAVE_PMMINTRIN_H)
#    if(HAVE_SMMINTRIN_H)
#        set(GMX_SSE4_1 1)
#    endif(HAVE_SMMINTRIN_H)

elseif(${GMX_ACCELERATION} STREQUAL "FORTRAN")
    set(GMX_FORTRAN 1)
elseif(${GMX_ACCELERATION} STREQUAL "BLUEGENE")
    set(GMX_BLUEGENE 1)
elseif(${GMX_ACCELERATION} STREQUAL "POWER6")
    set(GMX_POWER6 1)
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
    MESSAGE(FATAL_ERROR "Unrecognized option for accelerated kernels: ${GMX_ACCELERATION}. Pick one of none, SSE, Fortran, BlueGene, Power6, ia64, altivec")
endif(${GMX_ACCELERATION} STREQUAL "NONE")







# Process QM/MM Settings
string(TOUPPER ${GMX_QMMM_PROGRAM} ${GMX_QMMM_PROGRAM})
if(${GMX_QMMM_PROGRAM} STREQUAL "GAUSSIAN")
    set(GMX_QMMM_GAUSSIAN 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "MOPAC")
    set(GMX_QMMM_MOPAC 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "GAMESS")
    set(GMX_QMMM_GAMESS 1)
elseif(${GMX_QMMM_PROGRAM} STREQUAL "NONE")
    # nothing to do
else(${GMX_QMMM_PROGRAM} STREQUAL "GAUSSIAN")
    MESSAGE(FATAL_ERROR "Invalid QM/MM program option: ${GMX_QMMM_PROGRAM}. Choose one of: Gaussian, Mopac, Gamess, None")
endif(${GMX_QMMM_PROGRAM} STREQUAL "GAUSSIAN")

# Process FFT library settings
string(TOUPPER ${GMX_FFT_LIBRARY} ${GMX_FFT_LIBRARY})
if(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
#    MESSAGE(STATUS "Using external FFT library - fftw3")
    if(GMX_DOUBLE)
        find_package(FFTW3 REQUIRED)
    else(GMX_DOUBLE)
        find_package(FFTW3F REQUIRED)
    endif(GMX_DOUBLE)

    if(NOT FFTW3_FOUND)
        MESSAGE(FATAL_ERROR "Cannot find fftw3 (with correct precision). Fix it, choose another FFT library, or use the Gromacs built-in fftpack (slower)!")
    endif(NOT FFTW3_FOUND)

    set(GMX_FFT_FFTW3 1)
    include_directories(${FFTW3_INCLUDE_DIR})
    list(APPEND	GMX_EXTRA_LIBRARIES ${FFTW3_LIBRARIES})

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
    list(APPEND	GMX_EXTRA_LIBRARIES ${FFTW2_LIBRARIES})

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
    list(APPEND	GMX_EXTRA_LIBRARIES ${MKL_LIBRARIES})

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
  set(GMX_EXTERNAL_BLAS FALSE CACHE BOOL "Use external LAPACK instead of built-in") 
  set(GMX_EXTERNAL_LAPACK FALSE CACHE BOOL "Use external LAPACK instead of built-in") 
endif(HAVE_MKL OR ACCELERATE_FRAMEWORK)

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
set(CMAKE_C_FLAGS ${GROMACS_C_FLAGS} CACHE STRING "Compiler flags")
set(CMAKE_EXE_LINKER_FLAGS ${GROMACS_LINKER_FLAGS} CACHE STRING "Linker flags")

########################################################################
# Specify install locations and which subdirectories to process        #
########################################################################
set(LIB_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/lib)
set(BIN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/bin)
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/gromacs)
set(MAN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/share/man)

set(GMXLIBDIR        ${DATA_INSTALL_DIR}/top)


add_subdirectory(share)
add_subdirectory(include)
add_subdirectory(man)

add_subdirectory(src)


