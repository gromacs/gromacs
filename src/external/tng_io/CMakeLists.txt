cmake_minimum_required(VERSION 2.8)

if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W2")
else()
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()

project(TNG_IO)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(BUILD_SHARED_LIBS "Enable shared libraries" ON)

option(TNG_BUILD_FORTRAN "Build Fortran compatible library and examples for testing" OFF)

option(TNG_BUILD_EXAMPLES "Build examples showing usage of the TNG API" ON)
option(TNG_BUILD_TEST "Build TNG testing binary." ON)
option(TNG_BUILD_COMPRESSION_TESTS "Build tests of the TNG compression library" OFF)

option(TNG_USE_OPENMP "Try to use the OpenMP library (if available)" OFF)
if(TNG_USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

find_package(ZLIB)

include(CheckIncludeFile)
check_include_file(inttypes.h   HAVE_INTTYPES_H)

add_subdirectory(src)

install(FILES include/tng_io.h DESTINATION include/)

#-- Add an Option to toggle the generation of the API documentation
option(TNG_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(TNG_BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
  #-- Configure the Template Doxyfile for our specific project
  configure_file(Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (Docs ALL
                                        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
                                        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
endif()

