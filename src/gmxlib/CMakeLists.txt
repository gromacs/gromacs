
# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.

file(GLOB GMXLIB_SOURCES *.c 
     selection/*.c trajana/*.c
     statistics/*.c nonbonded/*.c nonbonded/nb_kernel_c/*.c)

if(GMX_IA32_SSE)
  if(GMX_ASM_USEASM-ATT)
    if(GMX_ASM_USECCOMPILER)
      SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
    endif(GMX_ASM_USECCOMPILER)
    enable_language(ASM-ATT)
    file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_ia32_sse/*.c nonbonded/nb_kernel_ia32_sse/*.s)
    file(GLOB SKIP_INTELSYNTAX nonbonded/nb_kernel_ia32_sse/*intel_syntax*)
    list(REMOVE_ITEM GMX_MORESSE_SOURCES ${SKIP_INTELSYNTAX})
    if(GMX_ASM_USECCOMPILER)
      set_source_files_properties(${GMX_MORESSE_SOURCES} PROPERTIES COMPILE_FLAGS "-c -m32")
    endif(GMX_ASM_USECCOMPILER)
  else(GMX_ASM_USEASM-ATT)
    enable_language(ASM)
    file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_ia32_sse/*.c nonbonded/nb_kernel_ia32_sse/*.intel_syntax.s)
  endif(GMX_ASM_USEASM-ATT)
endif(GMX_IA32_SSE)

if(GMX_X86_64_SSE)
  if(GMX_ASM_USEASM-ATT)
    if(GMX_ASM_USECCOMPILER)
      SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
    endif(GMX_ASM_USECCOMPILER)
    enable_language(ASM-ATT)
    file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_x86_64_sse/*.c nonbonded/nb_kernel_x86_64_sse/*.s)
    file(GLOB SKIP_INTELSYNTAX nonbonded/nb_kernel_x86_64_sse/*intel_syntax*)
    list(REMOVE_ITEM GMX_MORESSE_SOURCES ${SKIP_INTELSYNTAX})
    if(GMX_ASM_USECCOMPILER)
      set_source_files_properties(${GMX_MORESSE_SOURCES} PROPERTIES COMPILE_FLAGS "-c")
    endif(GMX_ASM_USECCOMPILER)
  else(GMX_ASM_USEASM-ATT)
    enable_language(ASM)
    file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_x86_64_sse/*.c nonbonded/nb_kernel_x86_64_sse/*.intel_syntax.s)
  endif(GMX_ASM_USEASM-ATT)
endif(GMX_X86_64_SSE)

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

add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${GMX_MORESSE_SOURCES} ${GMX_SSE2_SOURCES} ${THREAD_MPI_SRC})
target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES}  ${THREAD_LIB})

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR})

