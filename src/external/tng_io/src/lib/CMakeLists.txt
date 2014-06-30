set(source_files tng_io.c md5.c)
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

# Append the required library dependencies
target_link_libraries(tng_io tng_compress)
if(TNG_BUILD_WITH_ZLIB)
  target_link_libraries(tng_io ${ZLIB_LIBRARIES})
endif()

# Use GNUInstallDirst to set paths on multiarch systems
include(GNUInstallDirs)

install(TARGETS tng_io
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
