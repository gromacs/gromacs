cmake_minimum_required(VERSION 3.18)
project(SZ3 VERSION 3.1.7)
include(GNUInstallDirs)
include(CTest)

option(BUILD_SHARED_LIBS "build shared libraries by default" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
  set(CMAKE_BUILD_TYPE ${default_build_type} CACHE
          STRING "Choose the type to build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
          "MinSizeRel" "RelWithDebInfo")
endif()

find_package(PkgConfig)
find_package(OpenMP)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/SZ3/version.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/SZ3/version.hpp)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
target_compile_features(${PROJECT_NAME}
    INTERFACE cxx_std_17
  )

target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)

option(SZ3_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" OFF)
option(SZ3_DEBUG_TIMINGS "print debug timing information" ON)

if(SZ3_DEBUG_TIMINGS)
  target_compile_definitions(${PROJECT_NAME} INTERFACE SZ3_DEBUG_TIMINGS=1)
else()
  target_compile_definitions(${PROJECT_NAME} INTERFACE SZ3_DEBUG_TIMINGS=0)
endif()

pkg_search_module(ZSTD IMPORTED_TARGET libzstd)
if(ZSTD_FOUND AND NOT SZ3_USE_BUNDLED_ZSTD)
  target_link_libraries(${PROJECT_NAME} INTERFACE PkgConfig::ZSTD)
else()
  add_subdirectory(tools/zstd)
  target_link_libraries(${PROJECT_NAME} INTERFACE zstd)
endif()

find_package(GSL)
if (GSL_FOUND)
  target_compile_definitions(${PROJECT_NAME} INTERFACE -DENABLE_GSL="TRUE")
  target_link_libraries(${PROJECT_NAME} INTERFACE GSL::gsl)
endif()

add_subdirectory(tools/sz3)
add_subdirectory(tools/sz3c)


option(BUILD_H5Z_FILTER "build the H5Z-SZ3 filter" OFF)
if(BUILD_H5Z_FILTER)
	add_subdirectory(tools/H5Z-SZ3)
endif()

option(BUILD_MDZ "build MDZ (MMD-SZ)" OFF)
if(BUILD_MDZ)
  add_subdirectory(tools/mdz)
endif()


install(TARGETS ${PROJECT_NAME} 
  EXPORT SZ3Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
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
  COMPATIBILITY AnyNewerVersion
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
  )
#export sz3 target for external use
export(TARGETS SZ3 FILE SZ3.cmake)
