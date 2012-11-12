# Download and build a suitable copy of FFTW.
# The GROMACS team won't distribute source or binaries linked to FFTW
# because we are choosing to be very clear about distributing only
# LGPL-licensed code, to suit requirements from our funding source.
#
# Input: FFTW variable contains the FFTW component to build,
#        either fftw or fftwf for double or single precision

set(GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION "" CACHE INTERNAL "Optimization flags for FFTW compilation")
if(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[[:space:]].*")
    message(FATAL_ERROR "An internal limitation of FFTW means GROMACS cannot build FFTW in a directory with whitespace in its name. Either use a system FFTW, build it yourself, or build GROMACS in a different location.")
endif()

if(NOT GMX_DOUBLE)
    set(GMX_BUILD_OWN_FFTW_PREC --enable-float)
endif()

# Always build a static lib, so it gets added to libmd and doesn't need to be installed
set(GMX_BUILD_OWN_FFTW_SHARED_FLAG --disable-shared --enable-static)
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS) # FFTW doesn't use -DPIC by default
    set(GMX_BUILD_OWN_FFTW_SHARED_FLAG ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} --with-pic)
endif()

# Testing shows FFTW configured with --enable-avx is never better than --enable-sse2, so we do the latter always.
if(${GMX_CPU_ACCELERATION} MATCHES "^(SSE|AVX)")
    set(GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION --enable-sse2 CACHE INTERNAL "Optimization flags for FFTW compilation")
endif()

# Machinery for running the external project
set(EXTERNAL_FFTW_VERSION 3.3.2)
set(EXTERNAL_FFTW_VERSION_MD5_SUM 6977ee770ed68c85698c7168ffa6e178)
include(ExternalProject)
if(CMAKE_VERSION VERSION_LESS 2.8.3)
    # Can't check MD5 sum with this CMake version
    message(WARNING "GROMACS is downloading FFTW ${EXTERNAL_FFTW_VERSION}, but your version of CMake is too old to allow GROMACS to check that the archive it receives is correct. GROMACS will build and link to FFTW anyway, but there is a security risk if you execute a GROMACS tool that calls this library. Either build your own FFTW, or update CMake to at least version 2.8.3.")
        ExternalProject_add(gmxfftw
            URL "http://www.fftw.org/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz"
            CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
            ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}
            ${GMX_BUILD_OWN_FFTW_PREC})
else()
    ExternalProject_add(gmxfftw
        URL "http://www.fftw.org/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz"
        URL_MD5 ${EXTERNAL_FFTW_VERSION_MD5_SUM}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
        ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}
        ${GMX_BUILD_OWN_FFTW_PREC})
endif()
externalproject_get_property(gmxfftw INSTALL_DIR)
# The dependency that triggers building the gmxfftw target gets made where libmd is constructed

string(TOUPPER "${FFTW}" UPPERFFTW)
string(TOLOWER "${FFTW}" LOWERFFTW)
string(REGEX REPLACE "fftw" "fftw3" FFTW_LIBNAME ${LOWERFFTW})
set(${UPPERFFTW}_LIBRARIES ${CMAKE_BINARY_PREFIX}/${INSTALL_DIR}/lib/lib${FFTW_LIBNAME}${CMAKE_STATIC_LIBRARY_SUFFIX} PARENT_SCOPE)
set(${UPPERFFTW}_INCLUDE_DIRS ${CMAKE_BINARY_PREFIX}/${INSTALL_DIR}/include PARENT_SCOPE)
set(${UPPERFFTW}_FOUND TRUE PARENT_SCOPE)
set(${UPPERFFTW}_HAVE_SIMD TRUE PARENT_SCOPE)

message(STATUS "The GROMACS-managed build of FFTW 3 will configure with the following optimizations: ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}")
