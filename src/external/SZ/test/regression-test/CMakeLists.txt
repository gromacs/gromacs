add_test(NAME regression_tests
	COMMAND test.sh ${CMAKE_SOURCE_DIR}/test/cunit-test/sz.config ${CMAKE_SOURCE_DIR}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
