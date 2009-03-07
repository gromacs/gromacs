
file(GLOB MDLIB_SOURCES *.c)

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library
file(GLOB_RECURSE TEST_DRIVERS *_test.c)
if(${TEST_DRIVERS})
    list(REMOVE_ITEM MDLIB_SOURCES ${TEST_DRIVERS})
endif(${TEST_DRIVERS})


add_library(md ${MDLIB_SOURCES})
target_link_libraries(md gmx)

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR})

