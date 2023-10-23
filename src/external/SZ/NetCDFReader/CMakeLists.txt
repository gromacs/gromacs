pkg_search_module(NETCDF  REQUIRED IMPORTED_TARGET netcdf)
add_library(
  netcdfsz
  src/NetCDFReader.c
  )

target_include_directories(netcdfsz
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:netcdfsz>
  )

target_link_libraries(netcdfsz
  PUBLIC SZ PkgConfig::NETCDF
  )

install(TARGETS netcdfsz EXPORT NETCDFSZ
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/netcdfsz)
export(TARGETS netcdfsz FILE HDF5SZ.cmake)
