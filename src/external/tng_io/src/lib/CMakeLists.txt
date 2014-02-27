if(TNG_BUILD_FORTRAN)
  add_library(tng_io tng_io.c md5.c tng_io_fortran.c)
else()
  add_library(tng_io tng_io.c md5.c)
endif()

set_property(TARGET tng_io PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set_property(TARGET tng_io PROPERTY ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

install(TARGETS tng_io
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

if(HAVE_INTTYPES_H)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
endif()

# This test is for md5. The TNG library itself determines the actual byte order -
# not just if it is small or big endian.
include(TestBigEndian)
test_big_endian(TNG_INTEGER_BIG_ENDIAN)
if(TNG_INTEGER_BIG_ENDIAN)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS TNG_INTEGER_BIG_ENDIAN)
endif()

if(TNG_BUILD_WITH_ZLIB)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS USE_ZLIB)
  target_link_libraries(tng_io tng_compress ${ZLIB_LIBRARIES})
else()
  target_link_libraries(tng_io tng_compress)
endif()
