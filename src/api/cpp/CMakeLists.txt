# This list file provides the gmxapi::gmxapi cmake module. gmxapi::tests is the lone
# submodule, but subdirectories include the public headers.

# Define the gmxapi library target. Dependent targets can use properties of the target
# for include directories and header files.

# Configure / copy header files from the source to the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi)
configure_file(gmxapi/context.h gmxapi)
configure_file(gmxapi/exceptions.h gmxapi)
configure_file(gmxapi/gmxapi.h gmxapi)
configure_file(gmxapi/md.h gmxapi)
configure_file(gmxapi/runner.h gmxapi)
configure_file(gmxapi/system.h gmxapi)
configure_file(gmxapi/version.h gmxapi)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi/md)
configure_file(gmxapi/md/runnerstate.h gmxapi/md)

# Populate version constants from CMake.
configure_file(apiversion.in.h apiversion.h)

# Add to install target: copy the public API headers from the source directory
# into the gmxapi header subdirectory.
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

# Using $<TARGET_OBJECTS:mdrun_objlib> doesn't work in all environments when linking a library.
# The object files in programs/mdrun/ are compiled without -fPIC by default and so may not be
# linkable into a shared object. I had to tell cmake to build mdrun_objlib relocatable...
add_library(gmxapi SHARED $<TARGET_OBJECTS:mdrun_objlib>)

target_include_directories(gmxapi PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

target_include_directories(gmxapi PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           )

# libgmxapi source files.
target_sources(
        gmxapi PRIVATE
        apiversion.h
        context.cpp
        gmxapi.cpp
        md-impl.h
        md.cpp
        mdproxy.cpp
#        runnerimpl.cpp
        runnerstate.cpp
        system.cpp
        system-builder.h
        system-impl.h
        version.cpp
        )

# Ref. https://cmake.org/Wiki/CMake_RPATH_handling
# use, i.e. don't skip the full RPATH for the build tree
set_target_properties(gmxapi PROPERTIES SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set_target_properties(gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)

# duplicates link path when installed in LIB_INSTALL_DIR with libgromacs
#set_target_properties(gmxapi PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")

set_target_properties(gmxapi PROPERTIES
                      INSTALL_RPATH_USE_LINK_PATH TRUE
                      SOVERSION ${GMXAPI_MAJOR}
                      VERSION ${GMXAPI_RELEASE}
                      )

target_link_libraries(gmxapi PRIVATE libgromacs)

# Install the gmxapi target and simultaneously define the export target for
# which CMake will create a helper file. Specify the directory for clients to
# add to their include path to be able to `#include "gmxapi/some_header.h"`
install(TARGETS gmxapi
        EXPORT gmxapiExport
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
        RUNTIME DESTINATION ${LIB_INSTALL_DIR}
        INCLUDES DESTINATION include
        COMPONENT gmxapi)


# Create the CMake exports file to help other projects build against libgmxapi
# as a CMake import target Gromacs::gmxapi.
install(EXPORT gmxapiExport
        NAMESPACE Gromacs::
        DESTINATION share/cmake/gromacs/
        )
add_library(Gromacs::gmxapi ALIAS gmxapi )

include(CMakePackageConfigHelpers)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
# These include directories were not getting set right in the tests directory...
#target_include_directories( test_impl PRIVATE
#                            ${CMAKE_CURRENT_SOURCE_DIR}
#                            $<TARGET_PROPERTY:gmxapi,INCLUDE_DIRECTORIES>
#                            ${CMAKE_CURRENT_BINARY_DIR}
#                            )
