link_directories(${TNG_IO_BINARY_DIR}/src/lib)

add_definitions(-DTNG_COMPRESS_FILES_DIR="${CMAKE_BINARY_DIR}/test_tng_compress_files/") # Directory where to write tng_compress test files.

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_tng_compress_files)

set(number 0)
set(numtests 78)

while( number LESS ${numtests})

math( EXPR number "${number} + 1" )

add_executable(test_tng_compress_gen${number} testsuite.c)
target_link_libraries(test_tng_compress_gen${number} tng_io)
if(UNIX)
target_link_libraries(test_tng_compress_gen${number} m)
endif()
list(APPEND gen${number}_build_definitions GEN)
list(APPEND gen${number}_build_definitions TESTPARAM="test${number}.h")
set_target_properties(test_tng_compress_gen${number} PROPERTIES COMPILE_DEFINITIONS "${gen${number}_build_definitions}")
add_dependencies(test_tng_compress_gen${number} test${number}.h)
set_property(TARGET test_tng_compress_gen${number} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/compression_tests)

add_executable(test_tng_compress_read${number} testsuite.c)
target_link_libraries(test_tng_compress_read${number} tng_io)
if(UNIX)
target_link_libraries(test_tng_compress_read${number} m)
endif()
list(APPEND read${number}_build_definitions TESTPARAM="test${number}.h")
set_target_properties(test_tng_compress_read${number} PROPERTIES COMPILE_DEFINITIONS "${read${number}_build_definitions}")
add_dependencies(test_tng_compress_read${number} test${number}.h)
set_property(TARGET test_tng_compress_read${number} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/compression_tests)

endwhile()

if(UNIX)
file(COPY test_tng_compress_write.sh DESTINATION ${CMAKE_BINARY_DIR}/bin/compression_tests)
file(COPY test_tng_compress_read.sh DESTINATION ${CMAKE_BINARY_DIR}/bin/compression_tests)
endif()
if(WIN32)
file(COPY test_tng_compress_write.bat DESTINATION ${CMAKE_BINARY_DIR}/bin/compression_tests)
file(COPY test_tng_compress_read.bat DESTINATION ${CMAKE_BINARY_DIR}/bin/compression_tests)
endif()
