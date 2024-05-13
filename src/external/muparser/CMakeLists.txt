# CMake based on work from @xantares
cmake_minimum_required (VERSION 3.16.3)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# By default, build in Release mode. Must appear before project() command
if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project(muParserProject)

# Bump versions on release
set(MUPARSER_VERSION_MAJOR 2)
set(MUPARSER_VERSION_MINOR 3)
set(MUPARSER_VERSION_PATCH 4)
set(MUPARSER_VERSION ${MUPARSER_VERSION_MAJOR}.${MUPARSER_VERSION_MINOR}.${MUPARSER_VERSION_PATCH})

# Build options
option(ENABLE_SAMPLES "Build the samples" ON)
option(ENABLE_OPENMP "Enable OpenMP for multithreading" ON)
option(ENABLE_WIDE_CHAR "Enable wide character support" OFF)
option(BUILD_SHARED_LIBS "Build shared/static libs" ON)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# Credit: https://stackoverflow.com/questions/2368811/how-to-set-warning-level-in-cmake/3818084
if(MSVC)
    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    # Update if necessary
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

FILE(GLOB_RECURSE MUPARSER_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp") #all .cpp
# The test file produces warnings with MSVC and was recommended for exclusion, see
# https://github.com/beltoforion/muparser/issues/142#issuecomment-2093561006
list(REMOVE_ITEM MUPARSER_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/muParserTest.cpp")
add_library(muparser ${MUPARSER_SOURCES})

# Use the headers in the build-tree or the installed ones
target_include_directories(muparser SYSTEM PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# This compiles the "DLL" interface (C API)
target_compile_definitions(muparser PRIVATE MUPARSER_DLL)

if (BUILD_SHARED_LIBS)
  target_compile_definitions(muparser PRIVATE MUPARSERLIB_EXPORTS)
  add_definitions( -DMUPARSERLIB_EXPORTS )
else ()
  target_compile_definitions(muparser PUBLIC MUPARSER_STATIC)
  add_definitions( -DMUPARSER_STATIC )
endif()

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  target_compile_definitions(muparser PRIVATE _DEBUG)
endif ()

if(ENABLE_OPENMP)
  target_compile_definitions(muparser PRIVATE MUP_USE_OPENMP)
  target_link_libraries(muparser PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_WIDE_CHAR)
  target_compile_definitions(muparser PUBLIC _UNICODE)
endif()

set_target_properties(muparser PROPERTIES
    VERSION ${MUPARSER_VERSION}
    SOVERSION ${MUPARSER_VERSION_MAJOR}
)

if(ENABLE_SAMPLES)
  add_executable(example1 samples/example1/example1.cpp)
  target_link_libraries(example1 muparser)

  add_executable(example2 samples/example2/example2.c)
  target_link_libraries(example2 muparser)
endif()

# The GNUInstallDirs defines ${CMAKE_INSTALL_DATAROOTDIR}
# See https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include (GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/muparser)

install(TARGETS muparser
    EXPORT muparser-export
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
)

FILE(GLOB_RECURSE MUPARSER_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h") #all .h
install(FILES
    ${MUPARSER_HEADERS}
    DESTINATION include
    COMPONENT Development
)

# Export the target under the build-tree (no need to install)
export(EXPORT muparser-export
    FILE "${CMAKE_BINARY_DIR}/muparser-targets.cmake"
    NAMESPACE muparser::
)
add_library(muparser::muparser ALIAS muparser)

# Export the installed target (typically for packaging)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/muparserConfigVersion.cmake"
    VERSION ${MUPARSER_VERSION}
    COMPATIBILITY AnyNewerVersion
)
configure_file(muparserConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/muparserConfig.cmake"
    COPYONLY
)
install(EXPORT muparser-export
    FILE muparser-targets.cmake
    NAMESPACE muparser::
    DESTINATION ${INSTALL_CONFIGDIR}
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/muparserConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/muparserConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
    COMPONENT Development
)

# Define variables for the pkg-config file
set(PACKAGE_NAME muparser)
if(ENABLE_WIDE_CHAR)
  set(PKG_CONFIG_FLAGS "-D_UNICODE")
endif(ENABLE_WIDE_CHAR)
configure_file(
    muparser.pc.in
    ${CMAKE_BINARY_DIR}/muparser.pc
    @ONLY
)
install(
    FILES ${CMAKE_BINARY_DIR}/muparser.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

include(CTest)

if (BUILD_TESTING)
    add_executable (t_ParserTest test/t_ParserTest.cpp)
    target_link_libraries(t_ParserTest muparser)
    add_test (NAME ParserTest COMMAND t_ParserTest)
endif()
