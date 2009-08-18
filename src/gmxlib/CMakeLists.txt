
# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.

file(GLOB GMXLIB_SOURCES *.c 
     statistics/*.c nonbonded/*.c nonbonded/nb_kernel_c/*.c)

if(GMX_SSE2)
  if(GMX_DOUBLE)
    file(GLOB GMX_SSE2_SOURCES nonbonded/nb_kernel_sse2_double/*.c)
  else(GMX_DOUBLE)
    file(GLOB GMX_SSE2_SOURCES nonbonded/nb_kernel_sse2_single/*.c)
  endif(GMX_DOUBLE)
endif(GMX_SSE2)

if(NOT GMX_EXTERNAL_BLAS)
  file(GLOB BLAS_SOURCES gmx_blas/*.c)
endif(NOT GMX_EXTERNAL_BLAS)

if(NOT GMX_EXTERNAL_LAPACK)
  file(GLOB LAPACK_SOURCES gmx_lapack/*.c)
endif(NOT GMX_EXTERNAL_LAPACK)

# This would be the standard way to include thread_mpi, but we want libgmx
# to link the functions directly
#if(GMX_THREADS)
#    add_subdirectory(thread_mpi)
#endif(GMX_THREADS)
#target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_MPI_LIB})

# Files called xxx_test.c are test drivers with a main() function for module xxx.c,
# so they should not be included in the library
file(GLOB_RECURSE NOT_GMXLIB_SOURCES *_test.c *\#*)
list(REMOVE_ITEM GMXLIB_SOURCES ${NOT_GMXLIB_SOURCES})  

add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${GMX_SSE2_SOURCES} ${THREAD_MPI_SRC})
target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES} )

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR})

