if(GMX_CONTRIBUTED_FFTW)
  string(TOUPPER "${FFTW}" FFTW)
  set(GMX_CONTRIBUTED_FFTW_OPTIMIZATION_CONFIGURATION "" CACHE INTERNAL "optimization flags for fftw compilation")
  if(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[[:space:]].*") # fftw can't be built in a dir with whitespace
    message(FATAL_ERROR "Can't build fftw in a directory with whitespace in its name")
  endif(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[[:space:]].*")
  set(GMX_CONTRIBUTED_FFTW_LIBNAME libfftw3)
  if(NOT GMX_DOUBLE)
    set(GMX_CONTRIBUTED_FFTW_LIBNAME "${GMX_CONTRIBUTED_FFTW_LIBNAME}f")
    set(GMX_CONTRIBUTED_FFTW_PREC --enable-float)
  endif(NOT GMX_DOUBLE)
  #always build a static lib, so it gets added to gmx libs and don't need to be installed
  set(GMX_CONTRIBUTED_FFTW_SHARED_FLAG --disable-shared --enable-static)
  if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS) #fftw don't use pic by default
      set(GMX_CONTRIBUTED_FFTW_SHARED_FLAG ${GMX_CONTRIBUTED_FFTW_SHARED_FLAG} --with-pic)
  endif (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS)
  if(${GMX_CPU_ACCELERATION} MATCHES "^(SSE|AVX)")
    set(GMX_CONTRIBUTED_FFTW_OPTIMIZATION_CONFIGURATION --enable-sse2 CACHE INTERNAL "optimization flags for fftw compilation")
  endif(${GMX_CPU_ACCELERATION} MATCHES "^(SSE|AVX)")
  include(ExternalProject)
  ExternalProject_add(fftwBuild
    URL "http://www.fftw.org/fftw-3.3.2.tar.gz" URL_MD5 6977ee770ed68c85698c7168ffa6e178
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fftwBuild-prefix/src/fftwBuild/configure
    ${GMX_CONTRIBUTED_FFTW_SHARED_FLAG} ${GMX_CONTRIBUTED_FFTW_OPTIMIZATION_CONFIGURATION}
    ${GMX_CONTRIBUTED_FFTW_PREC} --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftwBuild-prefix)
  add_library(gmxfftw STATIC IMPORTED GLOBAL)
  add_dependencies(gmxfftw fftwBuild)
  set_target_properties(gmxfftw PROPERTIES IMPORTED_LOCATION  ${CMAKE_CURRENT_BINARY_DIR}/fftwBuild-prefix/lib/${GMX_CONTRIBUTED_FFTW_LIBNAME}.a)
  set(${FFTW}_LIBRARIES gmxfftw PARENT_SCOPE)
  set(${FFTW}_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/fftwBuild-prefix/include PARENT_SCOPE)
  set(${FFTW}_FOUND TRUE PARENT_SCOPE)
  set(${FFTW}_HAVE_SIMD TRUE PARENT_SCOPE)
  message("Contributed fftw (must not be distributed binarily!) will be build with the following optimizations: ${GMX_CONTRIBUTED_FFTW_OPTIMIZATION_CONFIGURATION}")
endif(GMX_CONTRIBUTED_FFTW)
