
file(GLOB MDLIB_SOURCES *.c)

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library
file(GLOB_RECURSE NOT_MDLIB_SOURCES *_test.c *\#*)
if(${NOT_MDLIB_SOURCES})
    list(REMOVE_ITEM MDLIB_SOURCES ${NOT_MDLIB_SOURCES})
endif(${NOT_MDLIB_SOURCES})


add_library(md ${MDLIB_SOURCES})
target_link_libraries(md gmx ${GMX_EXTRA_LIBRARIES})

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR})

