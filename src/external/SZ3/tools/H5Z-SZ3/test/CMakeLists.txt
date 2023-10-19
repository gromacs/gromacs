function(build_hdf5_test)
    get_filename_component(test_name ${ARGV0} NAME_WE)
    add_executable(${test_name} ${ARGV})
    target_link_libraries(${test_name} PUBLIC SZ3 hdf5sz3)
endfunction(build_hdf5_test)

build_hdf5_test(print_h5repack_args.cpp)
build_hdf5_test(sz3ToHDF5.cpp)
build_hdf5_test(dsz3FromHDF5.cpp)
build_hdf5_test(convertBinToHDF5.cpp)

