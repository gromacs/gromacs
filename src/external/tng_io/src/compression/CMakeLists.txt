add_library(tng_compress bwlzh.c bwt.c coder.c dict.c fixpoint.c huffman.c huffmem.c lz77.c merge_sort.c mtf.c rle.c tng_compress.c vals16.c warnmalloc.c widemuldiv.c xtc2.c xtc3.c)
if(UNIX)
target_link_libraries(tng_compress m)
endif()

set_property(TARGET tng_compress PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set_property(TARGET tng_compress PROPERTY ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

install(TARGETS tng_compress
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
