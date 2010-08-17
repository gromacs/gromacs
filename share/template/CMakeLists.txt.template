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
  set(LIBMD "libmd_d")
else(GMX_DOUBLE)
  set(LIBMD "libmd")
endif(GMX_DOUBLE)

#only check for libmd, because it depends on libgmx
FIND_PACKAGE(PkgConfig)
#do not use REQUIRED here, we have a better error message then pkgconfig
PKG_CHECK_MODULES(LIBMD ${LIBMD})
if (NOT LIBMD_FOUND)
  message(FATAL_ERROR "${LIBMD} not found, source GMXRC.")
endif (NOT LIBMD_FOUND)
add_definitions( ${LIBMD_CFLAGS} )
include_directories( ${LIBMD_INCLUDE_DIRS} )
link_directories( ${LIBMD_LIBRARY_DIRS} )


########################################################################
# Specify install locations and which subdirectories to process        #
########################################################################
set(BIN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/bin)

add_executable(template template.c)
target_link_libraries(template ${LIBMD_LIBRARIES})
