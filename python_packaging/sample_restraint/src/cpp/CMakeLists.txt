# Defines targets for the C++ restraints implemented here. These CMake targets are used by the
# unit tests and by the Python module target defined in ../pythonmodule/CMakeLists.txt

# Create a shared object library for our restrained ensemble plugin.
add_library(gmxapi_extension_ensemblepotential STATIC
            ensemblepotential.h
            ensemblepotential.cpp
            sessionresources.cpp)
set_target_properties(gmxapi_extension_ensemblepotential PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(gmxapi_extension_ensemblepotential PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(gmxapi_extension_ensemblepotential PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(gmxapi_extension_ensemblepotential PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(gmxapi_extension_ensemblepotential PRIVATE Gromacs::gmxapi)
