# - Find FFTW3
# Find the native FFTW3 includes and library, double precision
#
#  FFTW3_INCLUDE_DIR    - where to find fftw3.h
#  FFTW3_LIBRARIES   - List of libraries when using FFTW.
#  FFTW3_FOUND       - True if FFTW found.
#
# The FFTW3 root installation directory can be provided in the FFTW3_ROOT_DIR

if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARIES)
  # Already in cache, be silent
  set (FFTW3_FIND_QUIETLY TRUE)
endif (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARIES)

file(TO_CMAKE_PATH "$ENV{FFTW3_ROOT_DIR}" _env_FFTW3_ROOT_DIR)

find_path (FFTW3_INCLUDE_DIR fftw3.h
                PATHS "${_env_FFTW3_ROOT_DIR}/include"
		CACHE STRING "Path to double precision FFTW3 headers")

find_library (FFTW3_LIBRARIES 
		NAMES fftw3
                PATHS "${_env_FFTW3_ROOT_DIR}/lib"
                      "${FFTW3_INCLUDE_DIR}/../lib" 
		CACHE STRING "Double precision FFTW3 libraries")

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
set(__MSG "Could not find FFTW3. Provide the fftw3 install directory in the FFTW3_ROOT_DIR environment variable.")
find_package_handle_standard_args (FFTW3 ${__MSG} FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

mark_as_advanced (FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)


