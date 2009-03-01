
file(GLOB MDLIB_SOURCES *.c)

add_library(md ${MDLIB_SOURCES})

target_link_libraries(md gmx)


