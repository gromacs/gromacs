# - Try to find EXTRAE
# Once done this will define
#  EXTRAE_FOUND - System has EXTRAE
#  EXTRAE_INCLUDE_DIRS - The EXTRAE include directories
#  EXTRAE_LIBRARIES - The libraries needed to use EXTRAE

find_path(EXTRAE_INCLUDE_DIR extrae_user_events.h)

find_library(EXTRAE_LIBRARY NAMES seqtrace)

set(EXTRAE_LIBRARIES ${EXTRAE_LIBRARY} )
set(EXTRAE_INCLUDE_DIRS ${EXTRAE_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set EXTRAE_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXTRAE  DEFAULT_MSG
                                  EXTRAE_LIBRARY EXTRAE_INCLUDE_DIR)

mark_as_advanced(EXTRAE_INCLUDE_DIR EXTRAE_LIBRARY )

