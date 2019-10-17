# - Find OPENBABEL3
# Find the OpenBabel version 2 library.
#
#  OPENBABEL3_INCLUDE_DIR - where to find openbabel/obconversion.h
#  OPENBABEL3_LIBRARIES   - List of libraries when using OPENBABEL3
#  OPENBABEL3_FOUND       - True if OPENBABEL3 was found

if(OPENBABEL3_INCLUDE_DIR)
  # Already in cache, be silent
  set(OPENBABEL3_FIND_QUIETLY TRUE)
endif(OPENBABEL3_INCLUDE_DIR)

if(NOT OPENBABEL3_INCLUDE_DIR)
  find_path(OPENBABEL3_INCLUDE_DIR openbabel3/openbabel/obconversion.h
    PATHS
    ${_obIncDir}
    ${GNUWIN32_DIR}/include
    $ENV{OPENBABEL3_INCLUDE_DIR}
  )
  if(OPENBABEL3_INCLUDE_DIR)
    set(OPENBABEL3_INCLUDE_DIR ${OPENBABEL3_INCLUDE_DIR}/openbabel3 CACHE INTERNAL "OpenBabel3 include directory")
  endif(OPENBABEL3_INCLUDE_DIR)
endif(NOT OPENBABEL3_INCLUDE_DIR)

find_library(OPENBABEL3_LIBRARIES NAMES openbabel openbabel3
  PATHS
  ${_obLinkDir}
  ${GNUWIN32_DIR}/lib
  $ENV{OPENBABEL3_LIBRARIES}
)

# handle the QUIETLY and REQUIRED arguments and set OPENBABEL3_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OPENBABEL3 DEFAULT_MSG OPENBABEL3_LIBRARIES OPENBABEL3_INCLUDE_DIR)

mark_as_advanced (OPENBABEL3_LIBRARIES OPENBABEL3_INCLUDE_DIR)
