find_package(HDF5 REQUIRED)
add_library(
  hdf5sz
  src/H5Z_SZ.c
  )
target_link_libraries(
  hdf5sz
  PUBLIC SZ ${HDF5_LIBRARIES}
  )
target_include_directories(
  hdf5sz
  PUBLIC
  ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:hdf5_sz>
  )

if(BUILD_TESTS)
  add_subdirectory(test)
endif()

install(TARGETS hdf5sz EXPORT HDF5SZ
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz)
export(TARGETS hdf5sz FILE HDF5SZ.cmake)
