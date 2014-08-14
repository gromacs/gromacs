file(GLOB MMSLAVE_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${MMSLAVE_SOURCES} PARENT_SCOPE)

set(MMSLAVE_PUBLIC_HEADERS
    mmslave.h)
#install(FILES ${MMSLAVE_PUBLIC_HEADERS}
#        DESTINATION ${INCL_INSTALL_DIR}/gromacs
#        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)

add_executable(mmslave_test mmslave_test.c)
target_link_libraries(mmslave_test libgromacs  ${GMX_EXE_LINKER_FLAGS})
set_target_properties(mmslave_test PROPERTIES OUTPUT_NAME "mmslave_test${GMX_BINARY_SUFFIX}")
install(TARGETS mmslave_test DESTINATION ${BIN_INSTALL_DIR} COMPONENT mmslave_test)

