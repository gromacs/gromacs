# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (CMINPACK)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

include(${PROJECT_SOURCE_DIR}/cmake/cminpack_utils.cmake)
# Set version and OS-specific settings
set(CMINPACK_VERSION 1.3.4 CACHE STRING "CMinpack version")
set(CMINPACK_SOVERSION 1 CACHE STRING "CMinpack API version")
DISSECT_VERSION()
GET_OS_INFO()

if (OS_LINUX)
  option (USE_FPIC "Use the -fPIC compiler flag." ON)
else (OS_LINUX)
  option (USE_FPIC "Use the -fPIC compiler flag." OFF)
endif (OS_LINUX)

option (SHARED_LIBS "Build shared libraries instead of static." OFF)
if (SHARED_LIBS)
  message (STATUS "Building shared libraries.")
  set (LIB_TYPE SHARED)
else (SHARED_LIBS)
  message (STATUS "Building static libraries.")
  set (LIB_TYPE STATIC)
  if(WIN32)
    add_definitions(-DCMINPACK_NO_DLL)
  endif(WIN32)
endif (SHARED_LIBS)

add_subdirectory (cmake)

set (cminpack_srcs
  cminpack.h cminpackP.h
  chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c
  dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c
  dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c covar1.c
  minpack.h
  chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c
  dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c
  dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c
  )
set (cminpack_hdrs
    cminpack.h minpack.h)

add_library (cminpack ${LIB_TYPE} ${cminpack_srcs})

install (TARGETS cminpack 
   LIBRARY DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
   ARCHIVE DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
   RUNTIME DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library)
install (FILES ${cminpack_hdrs} DESTINATION ${CMINPACK_INCLUDE_INSTALL_DIR}
    COMPONENT cminpack_hdrs)

if (USE_FPIC AND NOT SHARED_LIBS)
  set_target_properties (cminpack PROPERTIES COMPILE_FLAGS -fPIC)
endif (USE_FPIC AND NOT SHARED_LIBS)

set_target_properties(cminpack PROPERTIES VERSION ${CMINPACK_VERSION} SOVERSION ${CMINPACK_SOVERSION})

