# - Find FFTW3
# Find the native FFTW3 includes and library, double precision
#
#  FFTW3_INCLUDE_DIR    - where to find fftw3.h
#  FFTW3_LIBRARIES   - List of libraries when using FFTW.
#  FFTW3_FOUND       - True if FFTW found.

if (FFTW3_INCLUDE_DIR)
  # Already in cache, be silent
  set (FFTW3_FIND_QUIETLY TRUE)
endif (FFTW3_INCLUDE_DIR)

find_path (FFTW3_INCLUDE_DIR fftw3.h)

find_library (FFTW3_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

mark_as_advanced (FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)


