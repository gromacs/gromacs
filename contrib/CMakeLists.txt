#TODO
# - hardcoded lib suffix in ${FFTW}_LIBRARIES
# - rename libfftw to libgmx_fftw to avoid collision with system lib

if(GMX_FFTW)
  # fftw limitation -- can't be built in
  # a directory with whitespace in its name.
  if(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[ \t].*")
    message(FATAL_ERROR
      "Can't build fftw in a directory with whitespace in its name")
  endif(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[ \t].*")
  set(FFTW_LIBNAME libfftw3)
  # build fftw as an external project
  if(NOT GMX_DOUBLE)
    set(FFTW_LIBNAME "${FFTW_LIBNAME}f")
    set(GMX_FFTW_PREC --enable-float)
  endif(NOT GMX_DOUBLE)
  if(BUILD_SHARED_LIBS)
    set(GMX_FFTW_SHARED_FLAG --enable-shared --disable-static)
  else(BUILD_SHARED_LIBS)
    set(GMX_FFTW_SHARED_FLAG --disable-shared --enable-static)
  endif(BUILD_SHARED_LIBS)
  if(${GMX_ACCELERATION} STREQUAL "SSE2")
    set(GMX_FFTW_OPTIMIZATION_CONFIGURATION --enable-sse2)
  elseif(${GMX_ACCELERATION} STREQUAL "AVX_128_FMA" OR ${GMX_ACCELERATION} STREQUAL "AVX_256")
    set(GMX_FFTW_OPTIMIZATION_CONFIGURATION --enable-avx)
  endif()
  include(ExternalProject)
  ExternalProject_add(fftw
    PREFIX fftw
    URL "http://www.fftw.org/fftw-3.3.2.tar.gz"
    URL_MD5 6977ee770ed68c85698c7168ffa6e178
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/fftw/configure
    ${GMX_FFTW_SHARED_FLAG}
    ${GMX_FFTW_OPTIMIZATION_CONFIGURATION}
    ${GMX_FFTW_PREC}
    --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw
    )
  set(${FFTW}_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/fftw/lib/${FFTW_LIBNAME}.so PARENT_SCOPE)
  set(${FFTW}_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/fftw/lib/${FFTW_LIBNAME}.so)
  set(${FFTW}_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/fftw/include PARENT_SCOPE)
  set(${FFTW}_FOUND TRUE PARENT_SCOPE)
  set(${FFTW}_HAVE_SIMD TRUE PARENT_SCOPE)
  set(GMX_FFTW_OPTIMIZATION_CONFIGURATION "" CACHE INTERNAL "optimazation flags for fftw compilation")
  file(GLOB GMX_FFTW_LIBRARIES ${${FFTW}_LIBRARIES}*) #get all .so files
  install(FILES ${GMX_FFTW_LIBRARIES} DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
endif(GMX_FFTW)
