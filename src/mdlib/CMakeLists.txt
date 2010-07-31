
file(GLOB MDLIB_SOURCES *.c)

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library
file(GLOB_RECURSE NOT_MDLIB_SOURCES *_test.c *\#*)
list(REMOVE_ITEM MDLIB_SOURCES ${NOT_MDLIB_SOURCES})

add_library(md ${MDLIB_SOURCES})
target_link_libraries(md gmx ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES})

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmd.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)
