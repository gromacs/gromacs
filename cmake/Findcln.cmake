# - Find CLN
# Find the Class Library for Numbers
#
#  CLN_INCLUDE_DIR - where to find cln/cln.h
#  CLN_LIBRARIES   - List of libraries when using CLN
#  CLN_FOUND       - True if CLN was found

if (CLN_INCLUDE_DIR)
  # Already in cache, be silent
  set (CLN_FIND_QUIETLY TRUE)
endif (CLN_INCLUDE_DIR)

find_path (CLN_INCLUDE_DIR cln/cln.h)
find_library (CLN_LIBRARIES cln)

# handle the QUIETLY and REQUIRED arguments and set CLN_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CLN DEFAULT_MSG CLN_LIBRARIES CLN_INCLUDE_DIR)

#mark_as_advanced (CLN_LIBRARIES CLN_INCLUDE_DIR)
