# ########################################################################
# Copyright 2013 Advanced Micro Devices, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ########################################################################

cmake_minimum_required( VERSION 3.1 )

if( CMAKE_GENERATOR MATCHES "NMake" )
	option( NMAKE_COMPILE_VERBOSE "Print compile and link strings to the console" OFF )
	if( NMAKE_COMPILE_VERBOSE )
		set( CMAKE_START_TEMP_FILE "" )
		set( CMAKE_END_TEMP_FILE "" )
		set( CMAKE_VERBOSE_MAKEFILE 1 )
	endif( )
endif( )

# This becomes the name of the solution file
project( clFFT )

# Define a version for the code
if( NOT DEFINED CLFFT_VERSION_MAJOR )
	set( CLFFT_VERSION_MAJOR 2 )
endif( )

if( NOT DEFINED CLFFT_VERSION_MINOR )
	set( CLFFT_VERSION_MINOR 14 )
endif( )

if( NOT DEFINED CLFFT_VERSION_PATCH )
	set( CLFFT_VERSION_PATCH 0 )
endif( )

set( CLFFT_VERSION "${CLFFT_VERSION_MAJOR}.${CLFFT_VERSION_MINOR}.${CLFFT_VERSION_PATCH}")

# This is incremented when the ABI to the library changes
set( CLFFT_SOVERSION 2 )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR} )

# On windows, it's convenient to change the default install prefix such that it does NOT point to 'program files'
# Need to check out CMAKE_RUNTIME_OUTPUT_DIRECTORY variable, and see if that eliminates the need to modify install path
if( WIN32 AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
	set( CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/package" CACHE PATH "Install path prefix, prepended onto install directories" FORCE )
endif( )

# Set the default of CMAKE_BUILD_TYPE to be release, unless user specifies with -D.  MSVC_IDE does not use CMAKE_BUILD_TYPE
if( NOT MSVC_IDE AND NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif()

# Options that the user or driver program can set to control various components of the build
option( BUILD_RUNTIME "Build the FFT runtime library" ON )
option( BUILD_CLIENT "Build a command line clFFT client program with a variety of configurable parameters (dependency on Boost)" OFF )
option( BUILD_TEST "Build the library testing suite (dependency on google test, Boost, and FFTW)" OFF )
option( BUILD_LOADLIBRARIES "Build the optional dynamic load libraries that the FFT runtime will search for" ON )
option( BUILD_SHARED_LIBS "Build shared libraries." ON )
option( BUILD_EXAMPLES "Build examples." ON )
option( BUILD_CALLBACK_CLIENT "Build a command line clFFT client program that tests callback functionality (dependency on Boost)" OFF )

# If BOOST_ROOT is defined as an environment value, use that value and cache it so it's visible in the cmake-gui.
# Otherwise, create a sensible default that the user can change
if( DEFINED ENV{BOOST_ROOT} )
	set( BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Environment variable defining the root of the Boost installation" )
endif( )

# Currently, linux has a problem outputing both narrow and wide characters,
# which happens in our client because openCL only supports narrow characters
if( WIN32 )
	option( UNICODE "Build with Unicode Support" ON )
	if( UNICODE )
		message( STATUS "UNICODE build" )
	endif( )
else( )
	set( CMAKE_SKIP_BUILD_RPATH ON )
	set( UNICODE OFF )
	message( STATUS "UNICODE feature disabled on linux" )
endif()

if( MSVC_IDE )
	set_property( GLOBAL PROPERTY USE_FOLDERS TRUE )
	set( BUILD64 ${CMAKE_CL_64} )
else()
	option( BUILD64 "Build a 64-bit product" ON )

	if( IS_DIRECTORY ${PROJECT_SOURCE_DIR}/tests )
		option( CODE_COVERAGE "Build makefiles with code coverage instrumentation" OFF )
		if( CODE_COVERAGE )
			message( STATUS "Code coverage instrumentation on" )
		endif()
	endif()
endif()

# These variables are meant to contain string which should be appended to the installation paths
# of library and executable binaries, respectively.  They are meant to be user configurable/overridable.
set( SUFFIX_LIB_DEFAULT "" )
set( SUFFIX_BIN_DEFAULT "" )

# Modify the global find property to help us find libraries like Boost in the correct paths for 64-bit
# Essentially, find_library calls will look for /lib64 instead of /lib; works for windows and linux
if( BUILD64 )
	set_property( GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE )
	message( STATUS "64bit build - FIND_LIBRARY_USE_LIB64_PATHS TRUE" )
	if( NOT APPLE )
		set( SUFFIX_LIB_DEFAULT "64" )
	endif()
else( )
	set_property( GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS FALSE )
	message( STATUS "32bit build - FIND_LIBRARY_USE_LIB64_PATHS FALSE" )
endif( )

set( SUFFIX_LIB ${SUFFIX_LIB_DEFAULT} CACHE STRING "String to append to 'lib' install path" )
set( SUFFIX_BIN ${SUFFIX_BIN_DEFAULT} CACHE STRING "String to append to 'bin' install path" )

# Useful variables to configure FindBoost.cake
# set( Boost_USE_MULTITHREADED ON )
# set( Boost_DETAILED_FAILURE_MSG   ON )
# set( Boost_DEBUG ON )
# set( Boost_NO_SYSTEM_PATHS ON )

# Client is built only if boost is found; on windows, we need vs10 or higher
# Find Boost on the system, and configure the type of boost build we want
if( NOT DEFINED Boost_USE_STATIC_LIBS )
	set( Boost_USE_STATIC_LIBS   ON )
endif( )

if( NOT DEFINED Boost_USE_STATIC_RUNTIME )
	set( Boost_USE_STATIC_RUNTIME OFF )
endif( )

# This will define Boost_FOUND
find_package( Boost 1.33.0 COMPONENTS program_options )
if( Boost_FOUND )
	message( STATUS "Boost_PROGRAM_OPTIONS_LIBRARY: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
else( )
	message( WARNING "Try setting Boost_DEBUG and Boost_DETAILED_FAILURE_MSG for more information" )
endif( )

# This will define OPENCL_FOUND
find_package( OpenCL )

# This will define FFTW_FOUND
find_package( FFTW )

# Enable building of the clACML client if both requested and all dependencies are found
if( BUILD_CLIENT AND Boost_FOUND )
	set( FFT_CLIENT ON )
else( )
	set( FFT_CLIENT OFF )
endif( )

# Enable building of the googletest unit test framework if requested and all dependencies are found
set( UNIT_TEST OFF )
if( BUILD_TEST )
	include(gtest.cmake)
	if( GTEST_FOUND AND Boost_FOUND AND FFTW_FOUND )
		set( UNIT_TEST ON )
	else( )
		message( "GoogleTest unit testing will NOT be built" )
	endif( )
endif( )

# Enable building of the callback client if requested and all dependencies are found
if( BUILD_CALLBACK_CLIENT AND Boost_FOUND AND FFTW_FOUND )
	set( FFT_CALLBACK_CLIENT ON )
else( )
	set( FFT_CALLBACK_CLIENT OFF )
endif( )

# FFLAGS depend on the compiler, grab the compiler name from the path
get_filename_component( C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME_WE )
# message( "C_COMPILER_NAME: " ${C_COMPILER_NAME} )
# message( "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER} )

if( APPLE )
	# avoid warning on newer cmake versions
	set( CMAKE_MACOSX_RPATH 0 )
endif()

# Set common compile and link options
if( MSVC )
	# Following options for nMake
	message( STATUS "Detected MSVS Ver: " ${MSVC_VERSION} )

	# CMake sets huge stack frames for windows, for whatever reason.  We go with compiler default.
	string( REGEX REPLACE "/STACK:[0-9]+" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )
	string( REGEX REPLACE "/STACK:[0-9]+" "" CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}" )
	string( REGEX REPLACE "/STACK:[0-9]+" "" CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}" )

elseif( CMAKE_COMPILER_IS_GNUCXX )
	message( STATUS "Detected GNU fortran compiler." )
	EXEC_PROGRAM( ${CMAKE_CXX_COMPILER} ARGS --version OUTPUT_VARIABLE vnum )
	STRING( REGEX REPLACE ".*([0-9])\\.([0-9])\\.([0-9]).*" "\\1\\2\\3" vnum ${vnum} )
	if( ${vnum} STREQUAL "452" )
		# we only want c++0x if we're using gcc 4.5.2
		set( CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}" )
	endif( )

	set( CMAKE_CXX_FLAGS "-pthread ${CMAKE_CXX_FLAGS}" )
	set( CMAKE_C_FLAGS "-pthread ${CMAKE_C_FLAGS}" )

	# For linux debug builds, define the same preprocessing symbols as win to keep it simple
	if( CMAKE_BUILD_TYPE MATCHES "Debug" )
		add_definitions( "/D_DEBUG" )
	endif( )

	if( BUILD64 )
		set( CMAKE_CXX_FLAGS "-m64 ${CMAKE_CXX_FLAGS}" )
		set( CMAKE_C_FLAGS "-m64 ${CMAKE_C_FLAGS}" )
	else( )
		set( CMAKE_CXX_FLAGS "-m32 ${CMAKE_CXX_FLAGS}" )
		set( CMAKE_C_FLAGS "-m32 ${CMAKE_C_FLAGS}" )
	endif( )

	if( CODE_COVERAGE )
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage" )
		set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage" )
	endif( )
endif( )

# If UNICODE is defined for microsoft compilers, pass extra definitions
if( MSVC AND UNICODE )
	add_definitions( "/DUNICODE /D_UNICODE" )
endif( )

# If UNICODE is defined for MinGW compilers, pass extra definitions
if( MINGW AND UNICODE )
	add_definitions( "-DUNICODE -D_UNICODE" )
endif( )

# Print out compiler flags for viewing/debug
message( STATUS "CMAKE_CXX_COMPILER flags: " ${CMAKE_CXX_FLAGS} )
message( STATUS "CMAKE_CXX_COMPILER debug flags: " ${CMAKE_CXX_FLAGS_DEBUG} )
message( STATUS "CMAKE_CXX_COMPILER release flags: " ${CMAKE_CXX_FLAGS_RELEASE} )
message( STATUS "CMAKE_CXX_COMPILER relwithdebinfo flags: " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO} )
message( STATUS "CMAKE_EXE_LINKER link flags: " ${CMAKE_EXE_LINKER_FLAGS} )

if( NOT BUILD_SHARED_LIBS )
	set( CLFFT_STATIC ON )
endif( )
# configure a header file to pass the CMake version settings to the source, and package the header files in the output archive
configure_file( "${PROJECT_SOURCE_DIR}/include/clFFT.version.h.in" "${PROJECT_BINARY_DIR}/include/clFFT.version.h" )
install( FILES
		"${PROJECT_BINARY_DIR}/include/clFFT.version.h"
		"include/clFFT.h"
		"include/clAmdFft.h"
		"include/clAmdFft.version.h"
	DESTINATION
		"./include" )

# Recurse into subdirectory and start building files there
if( BUILD_RUNTIME AND IS_DIRECTORY "${PROJECT_SOURCE_DIR}/library" )
	add_subdirectory( library )
else( )
	message( "Runtime library will NOT be built" )
endif( )

if( IS_DIRECTORY "${PROJECT_SOURCE_DIR}/scripts/perf" )
	add_subdirectory( scripts/perf )
endif( )

# We only want to build the following if the user options are set
if( FFT_CLIENT AND IS_DIRECTORY "${PROJECT_SOURCE_DIR}/client" )
	add_subdirectory( client )
else( )
	message( "FFT clients will NOT be built" )
endif( )

# Recurse into subdirectory and start building files there
if( BUILD_LOADLIBRARIES AND IS_DIRECTORY "${PROJECT_SOURCE_DIR}/statTimer" )
	add_subdirectory( statTimer )
else( )
	message( "LoadLibraries will NOT be built" )
endif( )

if( UNIT_TEST AND IS_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" )
	# enable_testing( )
	add_subdirectory( tests )
else( )
	message( "GoogleTest unit tests will NOT be built" )
endif( )

# We only want to build the following if the user options are set
if( FFT_CALLBACK_CLIENT AND IS_DIRECTORY "${PROJECT_SOURCE_DIR}/callback-client" )
	add_subdirectory( callback-client )
else( )
	message( "FFT callback client will NOT be built" )
endif( )

if( BUILD_EXAMPLES )
	add_subdirectory( examples )
endif( )

if( WIN32 )
	set( destdir CMake )
else( )
	set( destdir lib${SUFFIX_LIB}/cmake/clFFT )
endif( )
string( REGEX REPLACE "[^/]+" ".." reldir "${destdir}" )
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/clFFTConfigVersion.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/clFFTConfigVersion.cmake
	@ONLY )
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/clFFTConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/clFFTConfig.cmake
	@ONLY )
install( EXPORT Library DESTINATION ${destdir} FILE clFFTTargets.cmake )
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/clFFTConfigVersion.cmake
	${CMAKE_CURRENT_BINARY_DIR}/clFFTConfig.cmake
	DESTINATION ${destdir} )

# The following code is setting variables to control the behavior of CPack to generate our
if( WIN32 )
	set( CPACK_SOURCE_GENERATOR "ZIP" )
	set( CPACK_GENERATOR "ZIP" )
else( )
	set( CPACK_SOURCE_GENERATOR "TGZ" )
	set( CPACK_GENERATOR "TGZ" )
endif( )

if( BUILD64 )
	set( CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CLFFT_VERSION}-${CMAKE_HOST_SYSTEM_NAME}-x64" )
else( )
	set( CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CLFFT_VERSION}-${CMAKE_HOST_SYSTEM_NAME}-x32" )
endif( )

set( CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CLFFT_VERSION}-${CMAKE_HOST_SYSTEM_NAME}-Source" )

set( CPACK_PACKAGE_VERSION_MAJOR ${CLFFT_VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${CLFFT_VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH ${CLFFT_VERSION_PATCH} )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenCL implementation of an FFT library" )
set( CPACK_PACKAGE_VENDOR "Neutral" )
set( CPACK_SOURCE_IGNORE_FILES "/\\\\.hg/;/\\\\.svn/;/\\\\.git/" )

# Define all variables that influence CPack before including CPack, such as install targets
include( CPack )
