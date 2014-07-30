# - Check if the source code provided in the SOURCE argument compiles.
# CHECK_C_SOURCE_COMPILES(SOURCE VAR)
# - macro which checks if the source code compiles
#  SOURCE   - source code to try to compile
#  VAR      - variable to store whether the source code compiled
#
# The following variables may be set before calling this macro to
# modify the way the check is run:
#
#  CMAKE_REQUIRED_FLAGS = string of compile command line flags
#  CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#  CMAKE_REQUIRED_INCLUDES = list of include directories
#  CMAKE_REQUIRED_LIBRARIES = list of libraries to link

MACRO(MY_CHECK_C_SOURCE_COMPILES SOURCE VAR)
  IF("${VAR}" MATCHES "^${VAR}$")
    SET(MACRO_CHECK_FUNCTION_DEFINITIONS
      "-D${VAR} ${CMAKE_REQUIRED_FLAGS}")
    IF(CMAKE_REQUIRED_LIBRARIES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_LIBRARIES
        "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}")
    ELSE(CMAKE_REQUIRED_LIBRARIES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_LIBRARIES)
    ENDIF(CMAKE_REQUIRED_LIBRARIES)
    IF(CMAKE_REQUIRED_INCLUDES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_INCLUDES
        "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}")
    ELSE(CMAKE_REQUIRED_INCLUDES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_INCLUDES)
    ENDIF(CMAKE_REQUIRED_INCLUDES)
    FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.c"
      "${SOURCE}\n")

    MESSAGE(STATUS "Performing Test ${VAR}")
    TRY_COMPILE(${VAR}
      ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.c
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
      "${CHECK_C_SOURCE_COMPILES_ADD_LIBRARIES}"
      "${CHECK_C_SOURCE_COMPILES_ADD_INCLUDES}"
      OUTPUT_VARIABLE OUTPUT)
    message(STATUS "output from my-compile:\n${OUTPUT}\n")
    IF(${VAR})
      SET(${VAR} 1 CACHE INTERNAL "Test ${VAR}")
      MESSAGE(STATUS "Performing Test ${VAR} - Success")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Performing C SOURCE FILE Test ${VAR} succeded with the following output:\n"
        "${OUTPUT}\n"
        "Source file was:\n${SOURCE}\n")
    ELSE(${VAR})
      MESSAGE(STATUS "Performing Test ${VAR} - Failed")
      SET(${VAR} "" CACHE INTERNAL "Test ${VAR}")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Performing C SOURCE FILE Test ${VAR} failed with the following output:\n"
        "${OUTPUT}\n"
        "Source file was:\n${SOURCE}\n")
    ENDIF(${VAR})
  ENDIF("${VAR}" MATCHES "^${VAR}$")
ENDMACRO(MY_CHECK_C_SOURCE_COMPILES)


# - Check whether the C compiler supports a given flag.
# CHECK_C_COMPILER_FLAG(<flag> <var>)
#  <flag> - the compiler flag
#  <var>  - variable to store the result
# This internally calls the check_c_source_compiles macro.
# See help for CheckCSourceCompiles for a listing of variables
# that can modify the build.

#=============================================================================
# Copyright 2006-2011 Kitware, Inc.
# Copyright 2006 Alexander Neundorf <neundorf@kde.org>
# Copyright 2011 Matthias Kretz <kretz@kde.org>
#
# ORIGINAL Copyright notice (from Copyright.txt):
#
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2011 Kitware, Inc., Insight Software Consortium
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

INCLUDE(CheckCSourceCompiles)

MACRO (CHECK_C_COMPILER_FLAG _FLAG _RESULT)
   SET(SAFE_CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}")
   SET(CMAKE_REQUIRED_DEFINITIONS "${_FLAG}")
   MY_CHECK_C_SOURCE_COMPILES("int main(void) { return 0; }" ${_RESULT}
     # Some compilers do not fail with a bad flag
     FAIL_REGEX "command line option .* is valid for .* but not for C" # GNU
     FAIL_REGEX "unrecognized .*option"                     # GNU
     FAIL_REGEX "unknown .*option"                          # Clang
     FAIL_REGEX "ignoring unknown option"                   # MSVC
     FAIL_REGEX "warning D9002"                             # MSVC, any lang
     FAIL_REGEX "option.*not supported"                     # Intel
     FAIL_REGEX "invalid argument .*option"                 # Intel
     FAIL_REGEX "ignoring option .*argument required"       # Intel
     FAIL_REGEX "[Uu]nknown option"                         # HP
     FAIL_REGEX "[Ww]arning: [Oo]ption"                     # SunPro
     FAIL_REGEX "command option .* is not recognized"       # XL
     FAIL_REGEX "command option .* contains an incorrect subargument" # XL
     FAIL_REGEX "not supported in this configuration; ignored"       # AIX
     FAIL_REGEX "File with unknown suffix passed to linker" # PGI
     FAIL_REGEX "WARNING: unknown flag:"                    # Open64
     FAIL_REGEX "Incorrect command line option:"            # Borland
     FAIL_REGEX "Warning: illegal option"                   # SunStudio 12
     FAIL_REGEX "[Ww]arning: Invalid"                       # Fujitsu
     )
   SET (CMAKE_REQUIRED_DEFINITIONS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
ENDMACRO (CHECK_C_COMPILER_FLAG)
