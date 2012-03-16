file(GLOB MDLIB_SOURCES *.c)

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library
file(GLOB_RECURSE NOT_MDLIB_SOURCES *_test.c)
list(REMOVE_ITEM MDLIB_SOURCES ${NOT_MDLIB_SOURCES})
set(MDLIB_SOURCES ${MDLIB_SOURCES} PARENT_SCOPE)
