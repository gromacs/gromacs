# - Find Intel MKL
# Find the Intel Math Kernel Library, version 6.0 or later
#
#  MKL_INCLUDE_DIR - where to find mkl_dfti.h
#  MKL_LIBRARIES   - List of libraries when using MKL
#  MKL_FOUND       - True if MKL was found

if (MKL_INCLUDE_DIR)
  # Already in cache, be silent
  set (MKL_FIND_QUIETLY TRUE)
endif (MKL_INCLUDE_DIR)

find_path (MKL_INCLUDE_DIR mkl_dfti.h)
find_library (MKL_LIBRARIES mkl_core)

# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIR)

# MKL Libraries change name ALL the time, so the user will have to set these...
# mark_as_advanced (MKL_LIBRARIES MKL_INCLUDE_DIR)


