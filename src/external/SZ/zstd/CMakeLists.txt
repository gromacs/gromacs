if (BUILD_SHARED_LIBS)
  set (BUILD_EXT_LIBS_TYPE "SHARED")
else ()
  set (BUILD_EXT_LIBS_TYPE "STATIC")
endif ()

add_library(zstd ${BUILD_EXT_LIBS_TYPE} 
  ./common/entropy_common.c
  ./common/pool.c
  ./common/threading.c
  ./common/debug.c
  ./common/xxhash.c
  ./common/fse_decompress.c
  ./common/zstd_common.c
  ./common/error_private.c
  ./compress/zstd_ldm.c
  ./compress/zstd_lazy.c
  ./compress/huf_compress.c
  ./compress/zstd_opt.c
  ./compress/zstd_double_fast.c
  ./compress/zstd_compress.c
  ./compress/zstd_fast.c
  ./compress/fse_compress.c
  ./compress/zstdmt_compress.c
  ./compress/hist.c
  ./decompress/zstd_decompress.c
  ./decompress/huf_decompress.c
  ./deprecated/zbuff_common.c
  ./deprecated/zbuff_compress.c
  ./deprecated/zbuff_decompress.c
  ./legacy/zstd_v05.c
  ./legacy/zstd_v04.c
  ./legacy/zstd_v06.c
  ./legacy/zstd_v07.c
  ./legacy/zstd_v03.c
  ./legacy/zstd_v02.c
  ./legacy/zstd_v01.c
  ./dictBuilder/cover.c
  ./dictBuilder/divsufsort.c
  ./dictBuilder/zdict.c
  )

target_include_directories(zstd
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/compress
    ${CMAKE_CURRENT_SOURCE_DIR}/decompress
    ${CMAKE_CURRENT_SOURCE_DIR}/deprecated
    ${CMAKE_CURRENT_SOURCE_DIR}/dictBuilder
    ${CMAKE_CURRENT_SOURCE_DIR}/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/legacy
  )
if (${BUILD_EXT_LIBS_TYPE} MATCHES "SHARED")
  set (LIB_RELEASE_NAME "zstd")
  set (LIB_DEBUG_NAME "zstd${CMAKE_DEBUG_POSTFIX}")
else ()
  if (WIN32 AND NOT MINGW)
    set (LIB_RELEASE_NAME "libzstd")
    set (LIB_DEBUG_NAME "libzstd${CMAKE_DEBUG_POSTFIX}")
  else ()
    set (LIB_RELEASE_NAME "zstd")
    set (LIB_DEBUG_NAME "zstd${CMAKE_DEBUG_POSTFIX}")
  endif ()
endif ()
set_target_properties (zstd PROPERTIES
    OUTPUT_NAME                ${LIB_RELEASE_NAME}
#    OUTPUT_NAME_DEBUG          ${LIB_DEBUG_NAME}
    OUTPUT_NAME_RELEASE        ${LIB_RELEASE_NAME}
    OUTPUT_NAME_MINSIZEREL     ${LIB_RELEASE_NAME}
    OUTPUT_NAME_RELWITHDEBINFO ${LIB_RELEASE_NAME}
  )


install(TARGETS zstd EXPORT ZSTDConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} 
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES 
  deprecated/zbuff.h 
  dictBuilder/zdict.h 
  zstd.h 
  common/zstd_errors.h 
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/
  )
export(TARGETS zstd FILE ZSTDConfig.cmake)


