# - Find SQLITE3
# Find the SQLite version 3 library.
#
#  SQLITE3_INCLUDE_DIR - where to find sqlite3.h
#  SQLITE3_LIBRARIES   - List of libraries when using SQLITE3
#  SQLITE3_FOUND       - True if SQLITE3 was found

if (SQLITE3_INCLUDE_DIR)
  # Already in cache, be silent
  set (SQLITE3_FIND_QUIETLY TRUE)
endif (SQLITE3_INCLUDE_DIR)

find_path (SQLITE3_INCLUDE_DIR sqlite3.h)
find_library (SQLITE3_LIBRARIES sqlite3)

# handle the QUIETLY and REQUIRED arguments and set SQLITE3_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SQLITE3 DEFAULT_MSG SQLITE3_LIBRARIES SQLITE3_INCLUDE_DIR)

#mark_as_advanced (SQLITE3_LIBRARIES SQLITE3_INCLUDE_DIR)
