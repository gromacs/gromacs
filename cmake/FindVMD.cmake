# The module defines the following variables:
#   VMD_EXECUTABLE - path to vmd command
#   VMD_FOUND - true if the command was found
# Example usage:
#   find_package(VMD)
#   if(VMD_FOUND)
#     message("vmd found: ${VMD_EXECUTABLE}")
#   endif()

# Look for 'vmd' or 'eg' (easy vmd)
#
set(vmd_names vmd)

# Prefer .cmd variants on Windows unless running in a Makefile
# in the MSYS shell.
#
if(WIN32)
  if(NOT CMAKE_GENERATOR MATCHES "MSYS")
    set(vmd_names vmd.cmd vmd)
  endif()
endif()

find_program(VMD_EXECUTABLE
  NAMES ${vmd_names}
  DOC "VMD command"
  )
mark_as_advanced(VMD_EXECUTABLE)

# Handle the QUIETLY and REQUIRED arguments and set VMD_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VMD DEFAULT_MSG VMD_EXECUTABLE)
