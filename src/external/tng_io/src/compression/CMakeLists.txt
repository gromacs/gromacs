set(source_files bwlzh.c bwt.c coder.c dict.c fixpoint.c huffman.c huffmem.c lz77.c merge_sort.c mtf.c rle.c tng_compress.c vals16.c warnmalloc.c widemuldiv.c xtc2.c xtc3.c)
add_library(tng_compress ${source_files})
add_library(tng_compress_objlib OBJECT ${source_files})
if(UNIX)
target_link_libraries(tng_compress m)
endif()

install(TARGETS tng_compress
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
