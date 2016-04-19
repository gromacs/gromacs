set(TNG_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
file(RELATIVE_PATH TNG_ROOT_BINARY_DIR ${CMAKE_SOURCE_DIR} ${TNG_ROOT_SOURCE_DIR})
set(TNG_ROOT_BINARY_DIR ${CMAKE_BINARY_DIR}/${TNG_ROOT_BINARY_DIR})

function (TNG_GENERATE_VERSION_H)
    set(TNG_MAJOR_VERSION "1")
    set(TNG_MINOR_VERSION "7")
    set(TNG_VERSION_PATCH_LEVEL "9")
    set(TNG_IO_VERSION "${TNG_MAJOR_VERSION}.${TNG_MINOR_VERSION}.${TNG_VERSION_PATCH_LEVEL}")
    set(TNG_API_VERSION "7")
    configure_file(${TNG_ROOT_SOURCE_DIR}/include/tng/version.h.in
                   ${TNG_ROOT_BINARY_DIR}/include/tng/version.h)

    set(TNG_MAJOR_VERSION ${TNG_MAJOR_VERSION} PARENT_SCOPE)
    set(TNG_IO_VERSION ${TNG_IO_VERSION} PARENT_SCOPE)
endfunction()

tng_generate_version_h()

include(TestBigEndian)
test_big_endian(TNG_INTEGER_BIG_ENDIAN)
include(CheckIncludeFile)
check_include_file(inttypes.h TNG_HAVE_INTTYPES_H)

# Get a list of source files for compiling TNG, and perhaps its zlib dependency
#
# Parameters:
#   BUILD_OWN_ZLIB   (input)  Boolean value whether to build zlib internally
#   TNG_SOURCELIST   (output) Name of variable into which to store a list of source
#                             files for compling TNG (plus perhaps zlib)
#   TNG_COMPILEDEFS  (output) Name of variable into which to store the set of required defines
#
macro(TNG_GET_SOURCE_LIST BUILD_OWN_ZLIB TNG_SOURCELIST TNG_COMPILEDEFS)
    include_directories(BEFORE ${TNG_ROOT_SOURCE_DIR}/include)
    include_directories(BEFORE ${TNG_ROOT_BINARY_DIR}/include)
    set(_tng_compression_sources bwlzh.c bwt.c coder.c dict.c fixpoint.c huffman.c huffmem.c lz77.c merge_sort.c mtf.c rle.c tng_compress.c vals16.c warnmalloc.c widemuldiv.c xtc2.c xtc3.c)
    set(_tng_io_sources tng_io.c md5.c)
    set(${TNG_SOURCELIST})
    set(${TNG_COMPILEDEFS})
    foreach(_file ${_tng_compression_sources})
        list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/src/compression/${_file})
    endforeach()
    foreach(_file ${_tng_io_sources})
        list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/src/lib/${_file})
    endforeach()
    if(TNG_BUILD_FORTRAN)
      list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io_fortran.c)
    endif()
    if(BUILD_OWN_ZLIB)
        # Add minimal necessary number of TNG source files
        foreach(_file adler32.c compress.c crc32.c deflate.c inffast.c inflate.c inftrees.c trees.c uncompr.c zutil.c)
            list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/external/zlib/${_file})
        endforeach()
        include_directories(BEFORE ${TNG_ROOT_SOURCE_DIR}/external/zlib)
    endif()
    if (TNG_HAVE_INTTYPES_H)
        list(APPEND ${TNG_COMPILEDEFS} USE_STD_INTTYPES_H)
    endif()
endmacro()

macro(TNG_SET_SOURCE_PROPERTIES)
    if (TNG_HAVE_INTTYPES_H)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io.c
                     APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
    endif()
    if (TNG_INTEGER_BIG_ENDIAN)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/md5.c
                     APPEND PROPERTY COMPILE_DEFINITIONS TNG_INTEGER_BIG_ENDIAN)
    endif()
endmacro()
