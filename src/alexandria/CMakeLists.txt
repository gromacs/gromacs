
file(GLOB ALEXANDRIA_SOURCES *.c *.cpp)

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library
file(GLOB_RECURSE NOT_ALEXANDRIA_SOURCES *_test.c *\#*)
list(REMOVE_ITEM ALEXANDRIA_SOURCES ${NOT_ALEXANDRIA_SOURCES})

add_library(alexandria ${ALEXANDRIA_SOURCES})
target_link_libraries(alexandria md gmx ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES})

install(TARGETS alexandria DESTINATION ${LIB_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)
