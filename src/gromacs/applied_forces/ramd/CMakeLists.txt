gmx_add_libgromacs_sources(
    ramd.cpp
    ramdforceprovider.cpp
    ramdoptions.cpp
)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()