if (BUILD_SHARED_LIBS)
  set (BUILD_EXT_LIBS_TYPE "SHARED")
else ()
  set (BUILD_EXT_LIBS_TYPE "STATIC")
endif ()

add_library (SZ ${BUILD_EXT_LIBS_TYPE}
  src/ArithmeticCoding.c
  src/ByteToolkit.c
  src/CacheTable.c
  src/callZlib.c
  src/CompressElement.c
  src/conf.c
  src/dataCompression.c
  src/dictionary.c
  src/DynamicByteArray.c
  src/DynamicDoubleArray.c
  src/DynamicFloatArray.c
  src/DynamicIntArray.c
  src/Huffman.c
  src/iniparser.c
  src/MultiLevelCacheTable.c
  src/MultiLevelCacheTableWideInterval.c
  src/pastri.c
  src/exafelSZ.c
  src/rw.c
  src/rwf.c
  src/sz.c
  src/szd_double.c
  src/szd_double_pwr.c
  src/szd_double_ts.c
  src/szd_float.c
  src/szd_float_pwr.c
  src/szd_float_ts.c
  src/szd_int16.c
  src/szd_int32.c
  src/szd_int64.c
  src/szd_int8.c
  src/sz_double.c
  src/sz_double_pwr.c
  src/sz_double_ts.c
  src/szd_uint16.c
  src/szd_uint32.c
  src/szd_uint64.c
  src/szd_uint8.c
  src/szf.c
  src/sz_float.c
  src/sz_float_pwr.c
  src/sz_float_ts.c
  src/sz_int16.c
  src/sz_int32.c
  src/sz_int64.c
  src/sz_int8.c
  src/sz_omp.c
  src/sz_uint16.c
  src/sz_uint32.c
  src/sz_uint64.c
  src/sz_uint8.c
  src/TightDataPointStorageD.c
  src/TightDataPointStorageF.c
  src/TightDataPointStorageI.c
  src/TypeManager.c
  src/utility.c
  src/VarSet.c
  src/sz_stats.c
)

target_include_directories(SZ 
  PRIVATE 
  ${CMAKE_BINARY_DIR}
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/sz>
  )
if (${BUILD_EXT_LIBS_TYPE} MATCHES "SHARED")
  set (LIB_RELEASE_NAME "SZ")
  set (LIB_DEBUG_NAME "SZ${CMAKE_DEBUG_POSTFIX}")
else ()
  if (WIN32 AND NOT MINGW)
    set (LIB_RELEASE_NAME "libSZ")
    set (LIB_DEBUG_NAME "libSZ${CMAKE_DEBUG_POSTFIX}")
  else ()
    set (LIB_RELEASE_NAME "SZ")
    set (LIB_DEBUG_NAME "SZ${CMAKE_DEBUG_POSTFIX}")
  endif ()
endif ()
set_target_properties (SZ PROPERTIES
    OUTPUT_NAME                ${LIB_RELEASE_NAME}
#    OUTPUT_NAME_DEBUG          ${LIB_DEBUG_NAME}
    OUTPUT_NAME_RELEASE        ${LIB_RELEASE_NAME}
    OUTPUT_NAME_MINSIZEREL     ${LIB_RELEASE_NAME}
    OUTPUT_NAME_RELWITHDEBINFO ${LIB_RELEASE_NAME}
  )

target_link_libraries (SZ PUBLIC ${ZLIB_dep} ${ZSTD_dep} m)

target_compile_options(SZ
	PRIVATE $<$<CONFIG:Debug>:-Wall -Wextra -Wpedantic -Wno-unused-parameter>
	)

if(BUILD_OPENMP)
  target_link_libraries(SZ PRIVATE OpenMP::OpenMP_C)
endif()
if(BUILD_PASTRI)
  target_compile_definitions(SZ PUBLIC HAVE_PASTRI)
endif()
if(BUILD_TIMECMPR)
  target_compile_definitions(SZ PUBLIC HAVE_TIMECMPR)
endif()
if(BUILD_RANDOMACCESS)
  target_compile_definitions(SZ PUBLIC HAVE_RANDOMACCESS)
endif()
if(BUILD_FORTRAN)
  enable_language(Fortran)
  target_sources(SZ PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rw_interface.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sz_interface.F90
  )
endif()
if(BUILD_STATS)
  target_compile_definitions(SZ PUBLIC HAVE_WRITESTATS)
endif()

install (TARGETS SZ ${thirdparty_export} EXPORT SZConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

install (EXPORT SZConfig DESTINATION share/SZ/cmake)
install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sz)
export(TARGETS SZ FILE SZ.cmake)
