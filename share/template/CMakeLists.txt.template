cmake_minimum_required(VERSION 2.6)

project(template)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(C)

option(GMX_DOUBLE "Use double precision" OFF)

########################################################################
# Fix stupid flags on MSVC
########################################################################
IF(CMAKE_GENERATOR MATCHES "Visual Studio")
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)   
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio")

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################

if (GMX_DOUBLE)
  set(LIBGROMACS "libgmx_d")
else(GMX_DOUBLE)
  set(LIBGROMACS "libgmx")
endif(GMX_DOUBLE)

FIND_PACKAGE(GROMACS COMPONENTS ${LIBGROMACS} REQUIRED)
message("GROMACS version ${GROMACS_VERSION_STRING} found")
if ("${GROMACS_VERSION_STRING}" VERSION_LESS "4.6" OR "${GROMACS_VERSION_STRING}" VERSION_GREATER "4.99")
  message(FATAL_ERROR "This template works only with GROMACS 4.6")
endif ("${GROMACS_VERSION_STRING}" VERSION_LESS "4.6" OR "${GROMACS_VERSION_STRING}" VERSION_GREATER "4.99")

add_definitions( ${GROMACS_DEFINITIONS} )
include_directories( ${GROMACS_INCLUDE_DIRS} )

########################################################################
# Specify install locations and which subdirectories to process        #
########################################################################
set(BIN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/bin)

add_executable(template template.c)
target_link_libraries(template ${GROMACS_LIBRARIES})
