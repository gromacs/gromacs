# This module looks for Pandoc, and sets PANDOC_EXECUTABLE to the
# location of its binary.
#
# It respects the variable Pandoc_FIND_QUIETLY

include(FindPackageHandleStandardArgs)

if(Pandoc_FIND_QUIETLY OR PANDOC_EXECUTABLE)
  set(PANDOC_FIND_QUIETLY TRUE)
endif()

find_program(PANDOC_EXECUTABLE
  NAMES pandoc
  DOC "Pandoc - a universal document converter")

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pandoc REQUIRED_VARS PANDOC_EXECUTABLE)

mark_as_advanced(PANDOC_EXECUTABLE)
