include(GNUInstallDirs)

install(TARGETS h5xx EXPORT h5xxTargets
    INCLUDES DESTINATION "include"
)
# copy header files
install(DIRECTORY "${PROJECT_SOURCE_DIR}/h5xx" DESTINATION "include")

include(CMakePackageConfigHelpers)
set(ConfigPackageLocation "${CMAKE_INSTALL_DATADIR}/h5xx/cmake")

install(EXPORT h5xxTargets
    FILE H5XXTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)

configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/H5XXConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/H5XXConfig.cmake"
    INSTALL_DESTINATION ${ConfigPackageLocation}
)

write_basic_package_version_file("H5XXConfigVersion.cmake"
    VERSION ${H5XX_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/H5XXConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/H5XXConfigVersion.cmake"
    DESTINATION ${ConfigPackageLocation}
)
