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

# Allow cross-compiles
if (TARGET_HOST)
    set(GMX_BUILD_OWN_FFTW_TARGET --host=${TARGET_HOST})
endif()

# Machinery for running the external project
set(EXTERNAL_FFTW_VERSION 3.3.2)
# cmake make eats slashes //// -> //
set(GMX_BUILD_OWN_FFTW_URL "http:////www.fftw.org/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz" CACHE PATH "URL from where to download fftw, (use an absolute path when offline)")
mark_as_advanced(GMX_BUILD_OWN_FFTW_URL)
include(ExternalProject)
# TODO in master branch - show this warning only on the first run
# by using gmx_check_if_changed_result from I21b791ab8e4f3 when
# that becomes available
message(WARNING "The GROMACS build will download FFTW ${EXTERNAL_FFTW_VERSION} as requested, but it will not know the file it receives is correct. GROMACS will build and link to FFTW anyway, but there is a possible security risk if you execute a GROMACS tool that calls this library. You can use\nmake gmxfftw\n to do the download and build, and then run\nmd5sum src/contrib/fftw/gmxfftw-prefix/src/fftw-3.3.2.tar.gz\nto see if it matches 6977ee770ed68c85698c7168ffa6e178. If so, everything is OK and you should proceed with the rest of the GROMACS build. Alternatively, you can follow the GROMACS installation instructions to build FFTW yourself.")
# TODO if/when CMake fixes http://www.cmake.org/Bug/view.php?id=14330
# (ie. at least version > 2.8.11.2), consider reverting to using an
# md5sum check to avoid needing the above warning
    ExternalProject_add(fftwBuild
        URL "${GMX_BUILD_OWN_FFTW_URL}"
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/lib --disable-fortran
        ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}
        ${GMX_BUILD_OWN_FFTW_PREC}
        ${GMX_BUILD_OWN_FFTW_TARGET})
externalproject_get_property(fftwBuild INSTALL_DIR)

string(TOUPPER "${FFTW}" UPPERFFTW)
string(TOLOWER "${FFTW}" LOWERFFTW)
string(REGEX REPLACE "fftw" "fftw3" FFTW_LIBNAME ${LOWERFFTW})
set(${UPPERFFTW}_LIBRARIES ${CMAKE_BINARY_PREFIX}/${INSTALL_DIR}/lib/lib${FFTW_LIBNAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(${UPPERFFTW}_INCLUDE_DIRS ${CMAKE_BINARY_PREFIX}/${INSTALL_DIR}/include PARENT_SCOPE)
set(${UPPERFFTW}_FOUND TRUE PARENT_SCOPE)
set(${UPPERFFTW}_HAVE_SIMD TRUE PARENT_SCOPE)

if(CMAKE_VERSION VERSION_GREATER 2.8.7) #add_library GLOBAL is broken in 2.8.7
    add_library(gmxfftw STATIC IMPORTED GLOBAL)
    set_target_properties(gmxfftw PROPERTIES IMPORTED_LOCATION ${${UPPERFFTW}_LIBRARIES})
    set(${UPPERFFTW}_LIBRARIES gmxfftw PARENT_SCOPE)
else()
    # The dependency that triggers building the gmxfftw target gets made where libmd is constructed
    add_custom_target(gmxfftw)
    set(${UPPERFFTW}_LIBRARIES ${${UPPERFFTW}_LIBRARIES} PARENT_SCOPE)
endif()
add_dependencies(gmxfftw fftwBuild)

message(STATUS "The GROMACS-managed build of FFTW 3 will configure with the following optimizations: ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}")
