configure_file(config.h.cmakein config.h)
configure_file(buildinfo.h.cmakein buildinfo.h)

include(../cmake/BuildManPages.cmake)

add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(kernel)

if(NOT GMX_FAHCORE)
  add_subdirectory(tools)
  add_subdirectory(ngmx)
  add_subdirectory(contrib)
endif(NOT GMX_FAHCORE)
