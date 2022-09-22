# The first argument sets the name of the target for CMake bookkeeping purposes. It is the default name of the
# shared object library that is produced, but we can override that. There is no particular reason to
# change "gmxapi_extension" here unless you need different CMake target names to build several modules in
# a single project.
pybind11_add_module(gmxapi_extension MODULE export_plugin.cpp)

# Set the name of the shared object library (and the name of the Python module) to "myplugin".
# If you change "myplugin" you must also change the argument to the macro ``PYBIND11_MODULE(myplugin, m)`` in
# export_plugin.cpp
set_target_properties(gmxapi_extension PROPERTIES OUTPUT_NAME myplugin)

# We can't easily/reliably let a debug build of a Python module have a "d" suffix and still be importable with the same
# name.
set_target_properties(gmxapi_extension PROPERTIES DEBUG_POSTFIX "")

# We expect to be building against an installed GROMACS that we will continue to dynamically link against at runtime.
set_target_properties(gmxapi_extension PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
set_target_properties(gmxapi_extension PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# The Python module requires the new library we wrote as well as the gmxapi that we found in the top-level
# CMakeLists.txt
target_link_libraries(gmxapi_extension PRIVATE
                      Gromacs::gmxapi
                      gmxapi_extension_ensemblepotential
                      gmxapi_extension_test
                      )

if (GMXAPI_EXTENSION_MAIN_PROJECT)
    install(TARGETS gmxapi_extension
            LIBRARY DESTINATION ${GMXPLUGIN_INSTALL_PATH}
            ARCHIVE DESTINATION ${GMXPLUGIN_INSTALL_PATH}
            RUNTIME DESTINATION ${GMXPLUGIN_INSTALL_PATH}
            )
endif ()
