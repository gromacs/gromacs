if(GMX_FFTW)
  # fftw limitation -- can't be built in
  # a directory with whitespace in its name.
  if(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[[:space:]].*")
    message(FATAL_ERROR
      "Can't build fftw in a directory with whitespace in its name")
  endif(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[[:space:]].*")
  set(FFTW_LIBNAME libfftw3)
  # build fftw as an external project
  if(NOT GMX_DOUBLE)
    set(FFTW_LIBNAME "${FFTW_LIBNAME}f")
    set(GMX_FFTW_PREC --enable-float)
  endif(NOT GMX_DOUBLE)
  #always build a static lib, so it gets added to gmx libs and don't need to be installed
  set(GMX_FFTW_SHARED_FLAG --disable-shared --enable-static)
  if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS)
      set(GMX_FFTW_SHARED_FLAG ${GMX_FFTW_SHARED_FLAG} --with-pic)
  endif()
  if(${GMX_ACCELERATION} STREQUAL "SSE2")
    set(GMX_FFTW_OPTIMIZATION_CONFIGURATION --enable-sse2 CACHE INTERNAL "optimization flags for fftw compilation")
  elseif(${GMX_ACCELERATION} STREQUAL "AVX_128_FMA" OR ${GMX_ACCELERATION} STREQUAL "AVX_256")
    set(GMX_FFTW_OPTIMIZATION_CONFIGURATION --enable-avx CACHE INTERNAL "optimization flags for fftw compilation")
  endif()
  include(ExternalProject)
  ExternalProject_add(fftwBuild
    PREFIX fftw
    URL "http://www.fftw.org/fftw-3.3.2.tar.gz"
    URL_MD5 6977ee770ed68c85698c7168ffa6e178
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/fftwBuild/configure
    ${GMX_FFTW_SHARED_FLAG}
    ${GMX_FFTW_OPTIMIZATION_CONFIGURATION}
    ${GMX_FFTW_PREC}
    --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw
    )
  add_library(gmxfftw STATIC IMPORTED GLOBAL)
  add_dependencies(gmxfftw fftwBuild)
  set_target_properties(gmxfftw PROPERTIES IMPORTED_LOCATION  ${CMAKE_CURRENT_BINARY_DIR}/fftw/lib/${FFTW_LIBNAME}.a)
  set(${FFTW}_LIBRARIES gmxfftw PARENT_SCOPE)
  set(${FFTW}_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/fftw/include PARENT_SCOPE)
  set(${FFTW}_FOUND TRUE PARENT_SCOPE)
  set(${FFTW}_HAVE_SIMD TRUE PARENT_SCOPE)
  set(GMX_FFTW_OPTIMIZATION_CONFIGURATION "" CACHE INTERNAL "optimization flags for fftw compilation")
endif(GMX_FFTW)
