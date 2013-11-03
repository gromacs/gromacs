# This macro gets the version string for a compiler, and also extracts
# compiler flags.
#
# Parameters:
#   LANGUAGE       - C or CXX, the compiler to check for
#   BUILD_COMPILER - [output variable] string with compiler path, ID and
#                    some compiler-provided information
#   BUILD_FLAGS    - [output variable] flags for the compiler
#
macro(get_compiler_info LANGUAGE BUILD_COMPILER BUILD_FLAGS)
    set(${BUILD_COMPILER} "${CMAKE_${LANGUAGE}_COMPILER} ${CMAKE_${LANGUAGE}_COMPILER_ID} ${CMAKE_${LANGUAGE}_COMPILER_VERSION}")
    string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
    set(${BUILD_FLAGS} "${CMAKE_${LANGUAGE}_FLAGS} ${CMAKE_${LANGUAGE}_FLAGS_${_build_type}}")
endmacro()
