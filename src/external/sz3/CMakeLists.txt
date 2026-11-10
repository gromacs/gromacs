cmake_minimum_required(VERSION 3.19)
project(SZ3 VERSION 3.4.0)

#data version defines the version of the compressed data format
#it is not always equal to the program version (e.g., SZ3 v3.1.0 and SZ3 v.3.1.1 may use the same data version of v.3.1.0)
#only update data version if the new version of the program changes compressed data format
set(SZ3_DATA_VERSION 3.4.0)

include(GNUInstallDirs)
# Off when SZ3 is added with add_subdirectory or FetchContent: the parent installs what it needs of SZ3.
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" SZ3_IS_TOP_LEVEL)
option(SZ3_INSTALL "install SZ3's headers, libraries, tools and CMake package" ${SZ3_IS_TOP_LEVEL})


# Cache defaults, so only at top level: set from a subproject they would apply to the parent's whole build.
if (SZ3_IS_TOP_LEVEL)
    option(BUILD_SHARED_LIBS "build shared libraries by default" ON)
    set(default_build_type "Release")
    if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
        set(CMAKE_BUILD_TYPE ${default_build_type} CACHE
                STRING "Choose the type to build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                "MinSizeRel" "RelWithDebInfo")
    endif ()
    option(BUILD_TESTING "Build tests" OFF)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif ()
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
else () #-Wunused-parameter
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreorder -Wsign-compare -Wsuggest-override -Wold-style-cast -Wcast-qual -Wextra-semi -Wnon-virtual-dtor -Wcast-align -Wunused-variable")
endif ()

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/include/SZ3/version.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/include/SZ3/version.hpp)

# SZ3::SZ3core is SZ3 with only the dependencies it cannot work without. SZ3::SZ3 adds the optional
# ones SZ3 was built with, such as OpenMP. SZ3's headers still use OpenMP in a consumer compiled with it.
add_library(SZ3core INTERFACE)
add_library(SZ3::SZ3core ALIAS SZ3core)
target_include_directories(
        SZ3core INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_compile_features(SZ3core
        INTERFACE cxx_std_17
)

add_library(${PROJECT_NAME} INTERFACE)
add_library(SZ3::SZ3 ALIAS ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} INTERFACE SZ3core)
find_package(OpenMP)
if (OpenMP_FOUND)
    # BUILD_INTERFACE only: SZ3Config.cmake names the OpenMP target for the language the consumer
    # actually enabled. Exported as it stands, OpenMP::OpenMP_CXX is a target a C-only consumer
    # cannot be given, and CMake ends the consumer's configure over it.
    target_link_libraries(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:OpenMP::OpenMP_CXX>)
endif ()

option(SZ3_DEBUG_TIMINGS "print debug timing information" OFF)
if (SZ3_DEBUG_TIMINGS)
    target_compile_definitions(SZ3core INTERFACE SZ3_DEBUG_TIMINGS=1)
else ()
    target_compile_definitions(SZ3core INTERFACE SZ3_DEBUG_TIMINGS=0)
endif ()

# "${MSVC}" reads as ON under MSVC and OFF everywhere else. Keep the quotes.
option(SZ3_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" "${MSVC}")
if (NOT SZ3_USE_BUNDLED_ZSTD)
    # The search SZ3Config.cmake repeats in the consumer's project.
    find_library(SZ3_ZSTD_LIBRARY NAMES zstd)
endif ()
if (SZ3_ZSTD_LIBRARY AND NOT SZ3_USE_BUNDLED_ZSTD)
    # BUILD_INTERFACE, and no include directory: both are paths on this machine only.
    target_link_libraries(SZ3core INTERFACE $<BUILD_INTERFACE:${SZ3_ZSTD_LIBRARY}>)
else ()
    # What this build uses, for SZ3Config.cmake: the option only states a preference.
    set(SZ3_USE_BUNDLED_ZSTD ON)
    add_subdirectory(tools/zstd)
    target_link_libraries(SZ3core INTERFACE sz3_zstd)
endif ()

# Only at top level: a parent's BUILD_TESTING is for its own tests, and SZ3's would fetch a second googletest.
if (BUILD_TESTING AND SZ3_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory(tools/test)
endif ()


option(BUILD_SZ3_BINARY "build SZ3 tools" ${SZ3_IS_TOP_LEVEL})
if (BUILD_SZ3_BINARY)
    add_subdirectory(tools/sz3)
    add_subdirectory(tools/sz3c)
endif ()

option(BUILD_H5Z_FILTER "build the H5Z-SZ3 filter" OFF)
if (BUILD_H5Z_FILTER)
    add_subdirectory(tools/H5Z-SZ3)
endif ()

option(BUILD_MDZ "build MDZ (MMD-SZ)" OFF)
if (BUILD_MDZ)
    add_subdirectory(tools/mdz)
endif ()

option(BUILD_PARAVIEW_PLUGIN "build the ParaView plugin" OFF)
if (BUILD_PARAVIEW_PLUGIN)
    add_subdirectory(tools/paraview)
endif ()

if (SZ3_INSTALL)
    install(TARGETS ${PROJECT_NAME} SZ3core
            EXPORT SZ3Targets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            PATTERN "version.hpp.in" EXCLUDE)
    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/include/SZ3/version.hpp"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/SZ3/"
    )
    install(EXPORT SZ3Targets NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
    include(CMakePackageConfigHelpers)
    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SZ3Config.cmake.in
            "${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake"
            INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
    )
    write_basic_package_version_file(
            "${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake"
            VERSION ${PROJECT_VERSION}
            COMPATIBILITY SameMajorVersion
    )
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
    )
endif ()
#export sz3 target for external use
export(TARGETS SZ3 SZ3core FILE SZ3.cmake)
