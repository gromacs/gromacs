# Defines targets for the C++ restraints implemented here. These CMake targets are used by the
# unit tests and by the Python module target defined in ../pythonmodule/CMakeLists.txt

# These targets depend both on Gromacs::libgromacs and on Gromacs::gmxapi.
# These dependencies are likely to evolve.
# In particular, Gromacs::libgromacs has been deprecated for some time.
# https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=API+restructuring

add_library(gmxapi_extension_resources STATIC
            sessionresources.cpp)
set_target_properties(gmxapi_extension_resources PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(gmxapi_extension_resources PUBLIC
                           ${CMAKE_CURRENT_SOURCE_DIR})
# Note that the sample_restraint gmxapi extension
# is generally built as its own standalone project,
# but can be built within the GROMACS project build tree when GMX_PYTHON_PACKAGE=ON.
# When build against an existing GROMACS installation, the target that provides the
# restraint API headers is `Gromacs::libgromacs`. When built within the GROMACS build
# tree, the target that provides the restraint API is `restraint`. This is subject to
# change in the context of #3152 and #3288.
if(GMXAPI_EXTENSION_MAIN_PROJECT)
    target_link_libraries(gmxapi_extension_resources PUBLIC
                          Gromacs::libgromacs
                          Gromacs::gmxapi)
else()
    target_link_libraries(gmxapi_extension_resources PUBLIC
                          restraint
                          Gromacs::gmxapi)
endif ()

# Create an archive library for our restrained ensemble plugin.
add_library(gmxapi_extension_ensemblepotential STATIC
            ensemblepotential.cpp)
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

target_link_libraries(gmxapi_extension_ensemblepotential PRIVATE
                      gmxapi_extension_resources
                      Gromacs::gmxapi)

# Note that the sample_restraint gmxapi extension
# is generally built as its own standalone project,
# but can be built within the GROMACS project build tree when GMX_PYTHON_PACKAGE=ON.
# When build against an existing GROMACS installation, the target that provides the
# restraint API headers is `Gromacs::libgromacs`. When built within the GROMACS build
# tree, the target that provides the restraint API is `restraint`. This is subject to
# change in the context of #3152 and #3288.
if(GMXAPI_EXTENSION_MAIN_PROJECT)
    target_link_libraries(gmxapi_extension_ensemblepotential PUBLIC
                          Gromacs::libgromacs)
else()
    target_link_libraries(gmxapi_extension_ensemblepotential PUBLIC
                          restraint)
endif()

# Create an archive library for a test plugin.
add_library(gmxapi_extension_test STATIC
            nullpotential.cpp
            )
set_target_properties(gmxapi_extension_test PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(gmxapi_extension_test PUBLIC
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           )
set_target_properties(gmxapi_extension_test PROPERTIES SKIP_BUILD_RPATH FALSE)
set_target_properties(gmxapi_extension_test PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
target_link_libraries(gmxapi_extension_test PRIVATE
                      gmxapi_extension_resources
                      Gromacs::gmxapi)
# Note that the sample_restraint gmxapi extension
# is generally built as its own standalone project,
# but can be built within the GROMACS project build tree when GMX_PYTHON_PACKAGE=ON.
# When build against an existing GROMACS installation, the target that provides the
# restraint API headers is `Gromacs::libgromacs`. When built within the GROMACS build
# tree, the target that provides the restraint API is `restraint`. This is subject to
# change in the context of #3152 and #3288.
if(GMXAPI_EXTENSION_MAIN_PROJECT)
    target_link_libraries(gmxapi_extension_test PRIVATE
                          Gromacs::libgromacs)
else()
    target_link_libraries(gmxapi_extension_test PRIVATE
                          restraint)
endif()
