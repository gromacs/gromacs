#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.

#  To help us fund GROMACS development, we humbly ask that you cite
#  the research papers on the package. Check out http://www.gromacs.org
#
#
# - Find FFTW 3
# Find the native FFTW headers and libraries.
#
# Usage:
#
#  find_package(FFTW COMPONENTS [fftw|fftwf])
#
#  fftw or fftwf selects double or single precision FFTW3, respectively
#
# Input:
#
#  FFTW_FIND_COMPONENTS - which of fftw or fftwf component is sought? (set by the find_package() command)
#  FFTW_PATH            - user-specified path prefix in which the FFTW include and lib directories reside
#
# Outputs:
#
#  ${FFTW}_INCLUDE_DIRS - where to find FFTW headers
#  ${FFTW}_LIBRARIES    - List of libraries when using FFTW.
#  ${FFTW}_PKG          - The name of the pkg-config package needed
#  ${FFTW}_HAVE_SIMD    - True if FFTW was build with SIMD support
#  ${FFTW}_FOUND        - True if FFTW was found
#  where ${FFTW} is FFTW or FFTWF
#
#  This file is part of Gromacs        Copyright (c) 2012

#
# Is the component valid for GROMACS?
#
list(LENGTH FFTW_FIND_COMPONENTS FFTW_NUM_COMPONENTS_WANTED)
if(${FFTW_NUM_COMPONENTS_WANTED} LESS 1)
  message(FATAL_ERROR "find_package(FFTW) requires you specify a component: fftw or fftwf.")
elseif(${FFTW_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "find_package(FFTW) can only find one FFTW component at a time.")
elseif(${FFTW_FIND_COMPONENTS} MATCHES "^fftw(f)?$")
  if (NOT FFTW_FIND_VERSION OR FFTW_FIND_VERSION EQUAL 3) # Find FFTW3 by default
    string(TOUPPER "${FFTW_FIND_COMPONENTS}" FFTW)
    string(REGEX REPLACE "fftw" "fftw3" FFTW_COMPONENT "${FFTW_FIND_COMPONENTS}")
    set(FFTW_FUNCTION_PREFIX "${FFTW_FIND_COMPONENTS}")
  else()
    message(FATAL_ERROR "find_package(FFTW) can only find FFTW version 3.")
  endif()
else()
  message(FATAL_ERROR "find_package(FFTW) does not support finding component ${FFTW_FIND_COMPONENTS}. Choose fftw or fftwf.")
endif()

#
# If we've got pkgconfig and we haven't already found FFTW, try pkgconfig for clues
#
find_package(PkgConfig)
if(NOT __pkg_config_checked_PKGCONFIG_FFTW OR NOT FFTW_LIBRARY)
  pkg_check_modules(PKGCONFIG_FFTW "${FFTW_COMPONENT}")
endif(NOT __pkg_config_checked_PKGCONFIG_FFTW OR NOT FFTW_LIBRARY)

#
# Find path and library
#
find_path(FFTW_INCLUDE_DIR "fftw3.h" HINTS "${FFTW_PATH}/include" ${PKGCONFIG_FFTW_INCLUDE_DIRS} DOC "Path to fftw3.h include file")
find_library(FFTW_LIBRARY NAMES "${FFTW_COMPONENT}" HINTS "${FFTW_PATH}/lib" ${PKGCONFIG_FFTW_LIBRARY_DIRS} DOC "Library for linking FFTW3")

#
# Set the variables we'll export to the parent scope (and we'd add any dependencies
# here if they existed).
#
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
set(FFTW_INCLUDE_DIRS "${FFTW_INCLUDE_DIR}")

# Better error message than find_package_handle_standard_args
if (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
  set(FFTW_FOUND TRUE)
else()
  if (NOT FFTW_LIBRARY)
    message("Could not find ${FFTW_COMPONENT} library named lib${FFTW_COMPONENT}, if it is installed, please specify its location in FFTW_LIBRARY by hand (e.g. cmake .. -DFFTW_LIBRARY='/path/to/lib${FFTW_COMPONENT}.so')")
  endif()
  if (NOT FFTW_INCLUDE_DIR)
    message("Could not the ${FFTW_COMPONENT} header fftw3.h, if it is installed, please specify its path in FFTW_INCLUDE_DIR by hand (e.g. cmake .. -DFFTW_INCLUDE_DIR='/path/to/include')")
  endif()
endif()

if (FFTW_FOUND)
  # The user could specify trash in FFTW_LIBRARY, so test if we can link it
  include(CheckLibraryExists)
  if (HAVE_LIBM)
    # Adding MATH_LIBRARIES here to allow static libs, which we'll link later for our own use, anyway
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif (HAVE_LIBM)
  check_library_exists("${FFTW_LIBRARIES}" "${FFTW_FUNCTION_PREFIX}_plan_r2r_1d" "" FOUND_FFTW_PLAN)

  if(NOT FOUND_FFTW_PLAN)
    message(FATAL_ERROR "Could not find ${FFTW_FUNCTION_PREFIX}_plan_r2r_1d in ${FFTW_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a) make sure you have specified all dependencies of ${FFTW_COMPONENT} in FFTW_LIBRARY by hand (e.g. cmake .. -DFFTW_LIBRARY='/path/to/lib${FFTW_COMPONENT}.so;/path/to/libm.so') !")
  endif(NOT FOUND_FFTW_PLAN)

  # In FFTW 3.3, the SSE function name has changed
  foreach(SIMD_FUNCTION ${FFTW_FUNCTION_PREFIX}_have_simd_sse2;${FFTW_FUNCTION_PREFIX}_have_simd_avx;${FFTW_FUNCTION_PREFIX}_have_simd_altivec;${FFTW_FUNCTION_PREFIX}_have_simd_neon;${FFTW_FUNCTION_PREFIX}_have_sse2;${FFTW_FUNCTION_PREFIX}_have_sse;${FFTW_FUNCTION_PREFIX}_have_altivec)
    check_library_exists("${FFTW_LIBRARIES}" "${SIMD_FUNCTION}" "" FFTW_HAVE_${SIMD_FUNCTION})
    if(FFTW_HAVE_${SIMD_FUNCTION})
      set(FFTW_HAVE_SIMD TRUE BOOL "If ${FFTW_COMPONENT} was built with SIMD support")
      break()
    endif(FFTW_HAVE_${SIMD_FUNCTION})
  endforeach()

  # Verify that FFTW is compiled with -fPIC (sometimes necessary for GROMACS to link shared libraries)
  if (CMAKE_OBJDUMP AND CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS)
      execute_process(COMMAND ${CMAKE_OBJDUMP} --reloc ${FFTW_LIBRARY} OUTPUT_VARIABLE FFTW_OBJDUMP)
      if (${FFTW_OBJDUMP} MATCHES "R_X86_64" # Should always be true for static libraries. Checks that objdump works properly and that the library isn't dynamic
              AND NOT ${FFTW_OBJDUMP} MATCHES "R_X86_64_PLT32")
          message(FATAL_ERROR "The FFTW library ${FFTW_LIBRARY} cannot be used to link shared libraries in GROMACS. Either provide a different FFTW library by setting FFTW_LIBRARY, recompile FFTW with \"--enable-shared\" or \"--with-pic\", or disable the building of shared libraries in GROMACS by setting BUILD_SHARED_LIBS to \"OFF\". Note: Disabling shared libraries can require up to ten times as much disk space when installing GROMACS.")
      endif()
  endif()
  set(CMAKE_REQUIRED_LIBRARIES)

endif(FFTW_FOUND)

set(FFTW_HAVE_SIMD FALSE BOOL "Whether FFTW was built with SIMD support")

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY FFTW_HAVE_SIMD)
