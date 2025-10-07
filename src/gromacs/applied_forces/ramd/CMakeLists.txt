gmx_add_libgromacs_sources(
    ramd.cpp
    ramdforceprovider.cpp
    ramdoptions.cpp
    ramdoutputprovider.cpp
    read_params.cpp
)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()