if (BUILD_SHARED_LIBS)
  set (BUILD_EXT_LIBS_TYPE "SHARED")
else ()
  set (BUILD_EXT_LIBS_TYPE "STATIC")
endif ()

add_library(ZLIB ${BUILD_EXT_LIBS_TYPE} 
  ./gzclose.c
  ./uncompr.c
  ./trees.c
  ./gzread.c
  ./adler32.c
  ./gzwrite.c
  ./compress.c
  ./inftrees.c
  ./crc32.c
  ./inffast.c
  ./zutil.c
  ./gzlib.c
  ./infback.c
  ./inflate.c
  ./deflate.c
  )

target_include_directories(ZLIB
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )
if (${BUILD_EXT_LIBS_TYPE} MATCHES "SHARED")
  set (LIB_RELEASE_NAME "ZLIB")
  set (LIB_DEBUG_NAME "ZLIB${CMAKE_DEBUG_POSTFIX}")
else ()
  if (WIN32 AND NOT MINGW)
    set (LIB_RELEASE_NAME "libZLIB")
    set (LIB_DEBUG_NAME "libZLIB${CMAKE_DEBUG_POSTFIX}")
  else ()
    set (LIB_RELEASE_NAME "ZLIB")
    set (LIB_DEBUG_NAME "ZLIB${CMAKE_DEBUG_POSTFIX}")
  endif ()
endif ()
set_target_properties (ZLIB PROPERTIES
    OUTPUT_NAME                ${LIB_RELEASE_NAME}
#    OUTPUT_NAME_DEBUG          ${LIB_DEBUG_NAME}
    OUTPUT_NAME_RELEASE        ${LIB_RELEASE_NAME}
    OUTPUT_NAME_MINSIZEREL     ${LIB_RELEASE_NAME}
    OUTPUT_NAME_RELWITHDEBINFO ${LIB_RELEASE_NAME}
  )

install(TARGETS ZLIB EXPORT ZLIBConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES zlib.h zconf.h DESTINATION include)
export(TARGETS ZLIB FILE ZLIBConfig.cmake)
