include(FindPkgConfig)
pkg_search_module(CUNIT REQUIRED IMPORTED_TARGET cunit)

add_library(cunit_extras CUnit_Array.c)
target_link_libraries(cunit_extras PUBLIC PkgConfig::CUNIT)

function(make_sz_cunit_test binary source)
	add_executable(${binary} ${source})
	target_link_libraries(${binary} PUBLIC SZ cunit_extras)
	add_test(NAME ${binary} COMMAND ${binary})
endfunction()

make_sz_cunit_test(test_ByteToolkit test_ByteToolkit.c)
make_sz_cunit_test(test_DynamicByteArray test_DynamicByteArray.c)
make_sz_cunit_test(test_DynamicFloatArray test_DynamicFloatArray.c)
make_sz_cunit_test(test_DynamicIntArray.c test_DynamicIntArray.c)
make_sz_cunit_test(test_dataCompression test_dataCompression.c)
make_sz_cunit_test(test_TypeManager test_TypeManager.c)
#make_sz_cunit_test(test_Consistent test_Consistent.cc)
#make_sz_cunit_test(test_Huffman test_Huffman.c)
#make_sz_cunit_test(test_rw test_rw.c)
#make_sz_cunit_test(test_conf test_conf.c)


