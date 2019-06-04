# Build a C++ library using libgmxapi to extend GROMACS.
# We process this subdirectory first to reaffirm that this code does not depend
# on other substantial infrastructure, such as a particular Python bindings
# framework or the Googletest framework.
add_subdirectory(cpp)


# Import the Googletest targets, unless the user has requested that the repository
# be downloaded instead (see tests/CMakeLists.txt)
if(NOT ${DOWNLOAD_GOOGLETEST})
    if(${BUILD_TESTING})
        add_subdirectory(external/googletest EXCLUDE_FROM_ALL)
    endif()
endif()


# Build a Python extension package from our new library.
add_subdirectory(pythonmodule)
