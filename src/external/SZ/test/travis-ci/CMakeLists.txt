add_test(NAME travis_ci
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test.sh ${CMAKE_BINARY_DIR}/bin
	)

add_custom_target( SZ_data ALL
	${CMAKE_CURRENT_SOURCE_DIR}/getData.sh
	COMMENT "Downloading Data For Testing"
	)

file(COPY
	${CMAKE_CURRENT_SOURCE_DIR}/sz.config
	DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
	)
