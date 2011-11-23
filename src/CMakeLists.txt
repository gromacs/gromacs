configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/config.h)


add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(kernel)

if(NOT GMX_FAHCORE)
  add_subdirectory(tools)
  add_subdirectory(ngmx)
  add_subdirectory(contrib)
endif(NOT GMX_FAHCORE)
