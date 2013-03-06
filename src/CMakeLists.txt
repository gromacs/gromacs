configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/config.h)


add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(kernel)

if(NOT GMX_FAHCORE)
  add_subdirectory(tools)
  add_subdirectory(ngmx)
endif(NOT GMX_FAHCORE)

# Set global variable from the local variable's value
set(TARGETS_NEEDED_BY_MANUAL_BUILD ${TARGETS_NEEDED_BY_MANUAL_BUILD} PARENT_SCOPE)
