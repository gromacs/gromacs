# - Find OPENBABEL2
# Find the SQLite version 3 library.
#
#  OPENBABEL2_INCLUDE_DIR - where to find sqlite3.h
#  OPENBABEL2_LIBRARIES   - List of libraries when using OPENBABEL2
#  OPENBABEL2_FOUND       - True if OPENBABEL2 was found

if (OPENBABEL2_INCLUDE_DIR)
  # Already in cache, be silent
  set (OPENBABEL2_FIND_QUIETLY TRUE)
endif (OPENBABEL2_INCLUDE_DIR)

find_path (OPENBABEL2_INCLUDE_DIR openbabel/obconversion.h)
if(OPENBABEL2_INCLUDE_DIR)
  set(OPENBABEL2_INCLUDE_DIR ${OPENBABEL2_INCLUDE_DIR}/openbabel-2.0)
endif(OPENBABEL2_INCLUDE_DIR)
	
find_library (OPENBABEL2_LIBRARIES openbabel)

# handle the QUIETLY and REQUIRED arguments and set OPENBABEL2_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OPENBABEL2 DEFAULT_MSG OPENBABEL2_LIBRARIES OPENBABEL2_INCLUDE_DIR)

mark_as_advanced (OPENBABEL2_LIBRARIES OPENBABEL2_INCLUDE_DIR)
