# - Find FFTW3F
# Find the native FFTW3 includes and library, single precision
#
#  FFTW3F_INCLUDE_DIR    - where to find fftw3.h
#  FFTW3F_LIBRARIES   - List of libraries when using FFTW.
#  FFTW3F_FOUND       - True if FFTW found.
#
# The FFTW3F root installation directory can be provided in the FFTW3F_ROOT_DIR

if (FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARIES)
  # Already in cache, be silent
  set (FFTW3F_FIND_QUIETLY TRUE)
endif (FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARIES)

file(TO_CMAKE_PATH "$ENV{FFTW3F_ROOT_DIR}" _env_FFTW3F_ROOT_DIR)

find_path (FFTW3F_INCLUDE_DIR fftw3.h 
            PATHS "${_env_FFTW3F_ROOT_DIR}/include"
	    CACHE STRING "Path to single precision FFTW3 headers")

find_library (FFTW3F_LIBRARIES 
		NAMES fftw3f
                PATHS "${_env_FFTW3F_ROOT_DIR}/lib"
                      "${FFTW3F_INCLUDE_DIR}/../lib" 
		CACHE STRING "Single precision FFTW3 libraries")

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
set(__MSG "Could not find FFTW3F. Provide the fftw3 install directory in the FFTW3F_ROOT_DIR environment variable.")
find_package_handle_standard_args (FFTW3F ${__MSG} FFTW3F_LIBRARIES FFTW3F_INCLUDE_DIR)

mark_as_advanced (FFTW3F_LIBRARIES FFTW3F_INCLUDE_DIR)
