add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(kernel)
add_subdirectory(tools)
add_subdirectory(ngmx)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/config.h)

