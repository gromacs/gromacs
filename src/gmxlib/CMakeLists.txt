include_directories(${CMAKE_CURRENT_SOURCE_DIR})


# add target that generates version.c every time a make is run
# only do this if we generate the version
if(USE_VERSION_H)
    add_custom_target(gmx_version ALL
            COMMAND ${CMAKE_COMMAND} 
                -D Git_EXECUTABLE="${Git_EXECUTABLE}"
                -D Git_VERSION="${Git_VERSION}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_C_CMAKEIN="${CMAKE_SOURCE_DIR}/src/gmxlib/version.c.cmakein"
                -D VERSION_C_OUT="${CMAKE_CURRENT_BINARY_DIR}/version.c"
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake 
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/gmxlib 
            DEPENDS ${CMAKE_SOURCE_DIR}/src/gmxlib/version.c.cmakein
            COMMENT "Generating version information")
endif()

# The nonbonded directory contains subdirectories that are only
# conditionally built, so we cannot use a GLOB_RECURSE here.
file(GLOB GMXLIB_SOURCES *.c 
     selection/*.c trajana/*.c
     statistics/*.c nonbonded/*.c nonbonded/nb_kernel_c/*.c)

# This source file is generated
file(GLOB VERSION_SOURCE version.c)

if(VERSION_SOURCE)
	list(REMOVE_ITEM GMXLIB_SOURCES ${VERSION_SOURCE})
endif(VERSION_SOURCE)

# add version.c to the list of sources and tell cmake that it is generated
if(USE_VERSION_H)
LIST(APPEND GMXLIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/version.c) # auto-generated
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/version.c 
                PROPERTIES GENERATED true)
endif()

if(GMX_IA32_SSE OR GMX_IA32_SSE2)
  if(GMX_ASM_USEASM-NASM)
    enable_language(ASM-NASM)
    if(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_ia32_sse2/*.c nonbonded/nb_kernel_ia32_sse2/*intel_syntax.s)
    else(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_ia32_sse/*.c nonbonded/nb_kernel_ia32_sse/*intel_syntax.s)
    endif(GMX_DOUBLE)
  else(GMX_ASM_USEASM-NASM)
    if(GMX_ASM_USECCOMPILER)
      SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
    endif(GMX_ASM_USECCOMPILER)
    enable_language(ASM-ATT)
    if(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_ia32_sse2/*.c nonbonded/nb_kernel_ia32_sse2/*sse2.s nonbonded/nb_kernel_ia32_sse2/*asm.s)
    else(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_ia32_sse/*.c nonbonded/nb_kernel_ia32_sse/*sse.s nonbonded/nb_kernel_ia32_sse/*asm.s)
    endif(GMX_DOUBLE)
    if(GMX_ASM_USECCOMPILER)
      set_source_files_properties(${GMX_MORESSE_SOURCES} PROPERTIES COMPILE_FLAGS "-c -m32")
    endif(GMX_ASM_USECCOMPILER)
  endif(GMX_ASM_USEASM-NASM)
endif(GMX_IA32_SSE OR GMX_IA32_SSE2)

if(GMX_X86_64_SSE OR GMX_X86_64_SSE2)
  if(GMX_ASM_USEASM-NASM)
    enable_language(ASM-NASM)
    if(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_x86_64_sse2/*.c nonbonded/nb_kernel_x86_64_sse2/*intel_syntax.s)
    else(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_x86_64_sse/*.c nonbonded/nb_kernel_x86_64_sse/*intel_syntax.s)
    endif(GMX_DOUBLE)
  else(GMX_ASM_USEASM-NASM)
    if(GMX_ASM_USECCOMPILER)
      SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
    endif(GMX_ASM_USECCOMPILER)
    enable_language(ASM-ATT)
    if(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_x86_64_sse2/*.c nonbonded/nb_kernel_x86_64_sse2/*sse2.s nonbonded/nb_kernel_x86_64_sse2/*asm.s)
    else(GMX_DOUBLE)
      file(GLOB GMX_MORESSE_SOURCES nonbonded/nb_kernel_x86_64_sse/*.c nonbonded/nb_kernel_x86_64_sse/*sse.s nonbonded/nb_kernel_x86_64_sse/*asm.s)
    endif(GMX_DOUBLE)
    if(GMX_ASM_USECCOMPILER)
      set_source_files_properties(${GMX_MORESSE_SOURCES} PROPERTIES COMPILE_FLAGS "-c")
    endif(GMX_ASM_USECCOMPILER)
  endif(GMX_ASM_USEASM-NASM)
endif(GMX_X86_64_SSE OR GMX_X86_64_SSE2)

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
# Selection has test_ instead of _test.  Removing here for special case, perhaps make general?
if(GMX_FAHCORE)
    file(GLOB SELECTION_TEST selection/test*)
    list(REMOVE_ITEM GMXLIB_SOURCES ${SELECTION_TEST})
endif(GMX_FAHCORE)

add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${GMX_MORESSE_SOURCES} ${GMX_SSE2_SOURCES} ${THREAD_MPI_SRC})
target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES}  ${THREAD_LIB})
add_dependencies(gmx gmx_version) 
set_target_properties(gmx PROPERTIES OUTPUT_NAME "gmx${GMX_BINARY_SUFFIX}")

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmx.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)
