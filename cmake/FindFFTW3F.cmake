# - Find FFTW3
# Find the native FFTW3 includes and library, single precision
#
#  FFTW3_INCLUDES    - where to find fftw3.h
#  FFTW3_LIBRARIES   - List of libraries when using FFTW.
#  FFTW3_FOUND       - True if FFTW found.

if (FFTW3_INCLUDES)
  # Already in cache, be silent
  set (FFTW3_FIND_QUIETLY TRUE)
endif (FFTW3_INCLUDES)

find_path (FFTW3_INCLUDES fftw3.h)

find_library (FFTW3_LIBRARIES NAMES fftw3f)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDES)

mark_as_advanced (FFTW3_LIBRARIES FFTW3_INCLUDES)
