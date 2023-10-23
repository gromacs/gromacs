function(build_hdf5_test)
  get_filename_component(test_name ${ARGV0} NAME_WE)
  add_executable(${test_name} ${ARGV})
  target_link_libraries(${test_name} PUBLIC SZ hdf5sz)
endfunction(build_hdf5_test)

build_hdf5_test(szToHDF5.c) 
build_hdf5_test(dszFromHDF5.c)
build_hdf5_test(convertBinToHDF5.c)
