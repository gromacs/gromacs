#Either use:
#find_package(HDF5 1.12.1 REQUIRED)
#or manually set the hdf5 paths:
#set(HDF5_LIBRARIES "/path/to/hdf5-1.12.1/hdf5/lib/libhdf5.so")
#set(HDF5_INCLUDE_DIRS "/path/to/hdf5-1.12.1/hdf5/include")

find_package(HDF5 REQUIRED)

#message("HDF5LIB PATH:" ${HDF5_LIBRARIES})
#message("HDF5INCLUDE PATH:" ${HDF5_INCLUDE_DIRS})

add_library(
        hdf5sz3
        src/H5Z_SZ3.cpp
)

target_link_libraries(
        hdf5sz3
        PUBLIC SZ3 ${HDF5_LIBRARIES}
)
if(${HDF5_IS_PARALLEL})
    find_package(MPI REQUIRED)
    target_link_libraries(hdf5sz3 PUBLIC MPI::MPI_CXX)
endif()

target_include_directories(
        hdf5sz3
        PUBLIC
        ${HDF5_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:hdf5_sz3>
)

add_subdirectory(test)

install(TARGETS hdf5sz3 EXPORT HDF5SZ3
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_sz3)
export(TARGETS hdf5sz3 FILE HDF5SZ3.cmake)
