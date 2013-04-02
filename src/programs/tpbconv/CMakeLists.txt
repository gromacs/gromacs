include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

set(exename tpbconv)
configure_file(../main.cpp main.cpp)
file(GLOB TPBCONV_SOURCES tpbconv.c ${CMAKE_CURRENT_BINARY_DIR}/main.cpp)

add_executable(tpbconv ${TPBCONV_SOURCES})
gmx_add_man_page(tpbconv)
target_link_libraries(tpbconv ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
set_target_properties(tpbconv PROPERTIES OUTPUT_NAME "tpbconv${GMX_BINARY_SUFFIX}")
install(TARGETS tpbconv DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
