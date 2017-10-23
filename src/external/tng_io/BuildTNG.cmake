set(TNG_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
file(RELATIVE_PATH TNG_ROOT_BINARY_DIR ${CMAKE_SOURCE_DIR} ${TNG_ROOT_SOURCE_DIR})
if ("${TNG_ROOT_BINARY_DIR}" MATCHES "^\.\.")
    set(TNG_ROOT_BINARY_DIR tng)
endif()
set(TNG_ROOT_BINARY_DIR ${CMAKE_BINARY_DIR}/${TNG_ROOT_BINARY_DIR})

set(TNG_MAJOR_VERSION "1")
set(TNG_MINOR_VERSION "8")
set(TNG_VERSION_PATCH_LEVEL "1")
set(TNG_IO_VERSION "${TNG_MAJOR_VERSION}.${TNG_MINOR_VERSION}.${TNG_VERSION_PATCH_LEVEL}")

function (TNG_GENERATE_VERSION_H)
    set(TNG_API_VERSION "8")
    configure_file(${TNG_ROOT_SOURCE_DIR}/include/tng/version.h.in
                   ${TNG_ROOT_BINARY_DIR}/include/tng/version.h)

    set(TNG_MAJOR_VERSION ${TNG_MAJOR_VERSION} PARENT_SCOPE)
    set(TNG_IO_VERSION ${TNG_IO_VERSION} PARENT_SCOPE)
endfunction()

include(TestBigEndian)
test_big_endian(TNG_INTEGER_BIG_ENDIAN)
include(CheckIncludeFile)
check_include_file(inttypes.h TNG_HAVE_INTTYPES_H)
include(CMakeParseArguments)

function(add_tng_io_library NAME)
    tng_generate_version_h()

    set(_tng_compression_sources
        bwlzh.c bwt.c coder.c dict.c fixpoint.c huffman.c huffmem.c
        lz77.c merge_sort.c mtf.c rle.c tng_compress.c vals16.c
        warnmalloc.c widemuldiv.c xtc2.c xtc3.c)
    set(_tng_io_sources tng_io.c md5.c)
    set(_sources)
    foreach(_file ${_tng_compression_sources})
        list(APPEND _sources ${TNG_ROOT_SOURCE_DIR}/src/compression/${_file})
    endforeach()
    foreach(_file ${_tng_io_sources})
        list(APPEND _sources ${TNG_ROOT_SOURCE_DIR}/src/lib/${_file})
    endforeach()
    if(TNG_BUILD_FORTRAN)
        list(APPEND _sources ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io_fortran.c)
    endif()

    set(_options OBJECT OWN_ZLIB)
    cmake_parse_arguments(ARG "${_options}" "" "" ${ARGN})

    set(_build_target ${NAME})
    set(_link_type PRIVATE)
    if (ARG_OBJECT)
        set(_build_target tng_io_obj)
        set(_link_type INTERFACE)
        add_library(${_build_target} OBJECT ${_sources})
        # PIC is only on by default for SHARED libraries, but in case the
        # object library is going to get used in such a library, the objects
        # should be compiled with PIC as well.
        if (BUILD_SHARED_LIBS)
            set_target_properties(${_build_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
        endif()
        add_library(${NAME} INTERFACE)
        target_sources(${NAME} INTERFACE $<TARGET_OBJECTS:tng_io_obj>)
    else()
        add_library(${NAME} ${_sources})
        set_target_properties(${NAME} PROPERTIES
                              VERSION ${TNG_IO_VERSION}
                              SOVERSION ${TNG_MAJOR_VERSION})
        target_include_directories(${NAME} INTERFACE $<INSTALL_INTERFACE:include>)
    endif()
    target_include_directories(${_build_target} PRIVATE
                               $<BUILD_INTERFACE:${TNG_ROOT_SOURCE_DIR}/include>
                               $<BUILD_INTERFACE:${TNG_ROOT_BINARY_DIR}/include>)
    target_include_directories(${NAME} INTERFACE
                               $<BUILD_INTERFACE:${TNG_ROOT_SOURCE_DIR}/include>
                               $<BUILD_INTERFACE:${TNG_ROOT_BINARY_DIR}/include>)

    if (UNIX)
        target_link_libraries(${NAME} ${_link_type} m)
    endif()

    if (ARG_OWN_ZLIB)
        set(_zlib_dir ${TNG_ROOT_SOURCE_DIR}/external/zlib)
        set(_zlib_sources)
        # Add minimal necessary number of TNG source files
        foreach(_file adler32.c compress.c crc32.c deflate.c inffast.c inflate.c inftrees.c trees.c uncompr.c zutil.c)
            list(APPEND _zlib_sources ${_zlib_dir}/${_file})
        endforeach()
        add_library(tng_io_zlib OBJECT ${_zlib_sources})
        if (BUILD_SHARED_LIBS)
            set_target_properties(tng_io_zlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
        endif()
        target_include_directories(tng_io_zlib PUBLIC ${_zlib_dir})
        target_include_directories(${_build_target} PRIVATE ${_zlib_dir})
        target_sources(${NAME} ${_link_type} $<TARGET_OBJECTS:tng_io_zlib>)
    else()
        target_link_libraries(${NAME} ${_link_type} ZLIB::ZLIB)
    endif()

    if (TNG_HAVE_INTTYPES_H)
        target_compile_definitions(${NAME} INTERFACE USE_STD_INTTYPES_H)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io.c
                     APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
    endif()
    if (TNG_INTEGER_BIG_ENDIAN)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/md5.c
                     APPEND PROPERTY COMPILE_DEFINITIONS TNG_INTEGER_BIG_ENDIAN)
    endif()
endfunction()
