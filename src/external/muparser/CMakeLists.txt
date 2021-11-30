# CMake based on work from @xantares
cmake_minimum_required (VERSION 3.1.0)
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
set(MUPARSER_VERSION_PATCH 2)
set(MUPARSER_VERSION ${MUPARSER_VERSION_MAJOR}.${MUPARSER_VERSION_MINOR}.${MUPARSER_VERSION_PATCH})

# Build options
#option(ENABLE_SAMPLES "Build the samples" ON)
#option(ENABLE_OPENMP "Enable OpenMP for multithreading" ON)
set(ENABLE_OPENMP ${GMX_OPENMP})
option(BUILD_SHARED_LIBS "Build shared/static libs" ON)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_SHARED_LIBRARY_CXX_FLAGS}")
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

# include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
add_library(muparser
    src/muParserBase.cpp
    src/muParserBytecode.cpp
    src/muParserCallback.cpp
    src/muParser.cpp
    src/muParserDLL.cpp
    src/muParserError.cpp
    src/muParserInt.cpp
    src/muParserTest.cpp
    src/muParserTokenReader.cpp
)
target_include_directories(muparser SYSTEM PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# this compiles the "DLL" interface (C API)
target_compile_definitions(muparser PRIVATE MUPARSER_DLL)

if (BUILD_SHARED_LIBS)
  target_compile_definitions(muparser PRIVATE MUPARSERLIB_EXPORTS)
else ()
  target_compile_definitions(muparser PUBLIC MUPARSER_STATIC)
endif()

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  target_compile_definitions(muparser PRIVATE _DEBUG)
endif ()

if(ENABLE_OPENMP)
  target_compile_definitions(muparser PRIVATE MUP_USE_OPENMP)
endif()
set_target_properties(muparser PROPERTIES
    VERSION ${MUPARSER_VERSION}
    SOVERSION ${MUPARSER_VERSION_MAJOR}
)

# Install the export set for use with the install-tree
export(TARGETS muparser FILE "${CMAKE_BINARY_DIR}/muparser-targets.cmake")

if(ENABLE_SAMPLES)
  add_executable(example1 samples/example1/example1.cpp)
  target_link_libraries(example1 muparser)

  add_executable(example2 samples/example2/example2.c)
  target_link_libraries(example2 muparser)
endif()

# The GNUInstallDirs defines ${CMAKE_INSTALL_DATAROOTDIR}
# See https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include (GNUInstallDirs)

install(TARGETS muparser
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
)

install(FILES
    include/muParserBase.h
    include/muParserBytecode.h
    include/muParserCallback.h
    include/muParserDef.h
    include/muParserDLL.h
    include/muParserError.h
    include/muParserFixes.h
    include/muParser.h
    include/muParserInt.h
    include/muParserTemplateMagic.h
    include/muParserTest.h
    include/muParserToken.h
    include/muParserTokenReader.h
    DESTINATION include
    COMPONENT Development
)

# Define variables for the pkg-config file
set(PACKAGE_NAME muparser)
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
enable_testing()

add_executable (t_ParserTest test/t_ParserTest.cpp)
target_link_libraries(t_ParserTest muparser)
if (ENABLE_OPENMP)
    target_link_libraries(t_ParserTest OpenMP::OpenMP_CXX)
endif()
add_test (NAME ParserTest COMMAND t_ParserTest)
if (BUILD_TESTING)
    add_dependencies(tests t_ParserTest)
endif()
