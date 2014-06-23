set(source_files tng_io.c md5.c
    ../compression/bwlzh.c ../compression/bwt.c ../compression/coder.c
    ../compression/dict.c ../compression/fixpoint.c ../compression/huffman.c
    ../compression/huffmem.c ../compression/lz77.c ../compression/merge_sort.c
    ../compression/mtf.c ../compression/rle.c ../compression/tng_compress.c
    ../compression/vals16.c ../compression/warnmalloc.c ../compression/widemuldiv.c
    ../compression/xtc2.c ../compression/xtc3.c)

if(TNG_BUILD_FORTRAN)
  list(APPEND source_files tng_io_fortran.c)
endif()

add_library(tng_io ${source_files})

# Now add source-file compilation properties to the source-file
# targets
if(HAVE_INTTYPES_H)
  set_property(SOURCE tng_io.c APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
endif()
if(TNG_BUILD_WITH_ZLIB)
  set_property(SOURCE tng_io.c APPEND PROPERTY COMPILE_DEFINITIONS USE_ZLIB)
endif()

# This test is for md5. The TNG library itself determines the actual byte order -
# not just if it is small or big endian.
include(TestBigEndian)
test_big_endian(TNG_INTEGER_BIG_ENDIAN)
if(TNG_INTEGER_BIG_ENDIAN)
  set_property(SOURCE md5.c APPEND PROPERTY COMPILE_DEFINITIONS TNG_INTEGER_BIG_ENDIAN)
endif()

if(TNG_BUILD_WITH_ZLIB)
  target_link_libraries(tng_io ${ZLIB_LIBRARIES})
endif()

if(UNIX)
  target_link_libraries(tng_io m)
endif()

set_target_properties(tng_io PROPERTIES VERSION ${TNG_IO_VERSION} SOVERSION ${MAJOR_VERSION})

# Create the tng_ioConfig.cmake and tng_ioConfigVersion.cmake files for the install tree
set(CONF_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
configure_file(                tng_io-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-config.cmake" @ONLY)
configure_file(                tng_io-configVersion.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-configVersion.cmake" @ONLY)

# Use GNUInstallDirst to set paths on multiarch systems
include(GNUInstallDirs)

# Install the tng_ioConfig.cmake and tng_ioConfigVersion.cmake
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-configVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/tng_io")

install(TARGETS tng_io
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
