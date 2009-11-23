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

if(GMX_DOUBLE)
  set(SSETYPE sse2)
else()
  set(SSETYPE sse)
endif()

if(GMX_IA32_ASM)
  file(GLOB GMX_SSEKERNEL_C_SRC   nonbonded/nb_kernel_ia32_${SSETYPE}/*.c)
  if(GMX_ASM_USEASM-NASM)
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_ia32_${SSETYPE}/*intel_syntax*.s)    
  else()
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_ia32_${SSETYPE}/*${SSETYPE}.s nonbonded/nb_kernel_ia32_${SSETYPE}/*asm.s)
  endif()
endif(GMX_IA32_ASM)

if(GMX_X86_64_ASM)
  file(GLOB GMX_SSEKERNEL_C_SRC   nonbonded/nb_kernel_x86_64_${SSETYPE}/*.c)
  if(GMX_ASM_USEASM-NASM)
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_x86_64_${SSETYPE}/*intel_syntax*.s)
  else()
    file(GLOB GMX_SSEKERNEL_ASM_SRC nonbonded/nb_kernel_x86_64_${SSETYPE}/*${SSETYPE}.s nonbonded/nb_kernel_x86_64_${SSETYPE}/*asm.s)
  endif()
endif(GMX_X86_64_ASM)

if(GMX_FORTRAN)
  if (GMX_DOUBLE)
    file(GLOB FORTRAN_SOURCES nonbonded/nb_kernel_f77_double/*.[cf])
  else(GMX_DOUBLE)
    file(GLOB FORTRAN_SOURCES nonbonded/nb_kernel_f77_single/*.[cf])
  endif(GMX_DOUBLE)
endif(GMX_FORTRAN)

if(GMX_POWER6)
  file(GLOB FORTRAN_SOURCES nonbonded/nb_kernel_power6/*.[cF])
endif(GMX_POWER6)

if(GMX_BLUEGENE)
# This ordering speeds up parallel compilation
# a lot, by doing the longer ones first (more or less).
  set(GMX_BLUEGENE_C_SRC
	nonbonded/nb_kernel_bluegene/nb_kernel334_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel332_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel333_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel331_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel330_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel324_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel322_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel323_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel321_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel320_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel314_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel312_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel313_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel311_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel310_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel304_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel302_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel303_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel301_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel300_bluegene.c

	nonbonded/nb_kernel_bluegene/nb_kernel234_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel232_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel233_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel231_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel230_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel224_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel222_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel223_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel221_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel220_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel214_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel212_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel213_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel211_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel210_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel204_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel202_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel203_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel201_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel200_bluegene.c

	nonbonded/nb_kernel_bluegene/nb_kernel134_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel132_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel133_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel131_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel130_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel124_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel122_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel123_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel121_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel120_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel114_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel112_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel113_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel111_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel110_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel104_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel102_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel103_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel101_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel100_bluegene.c

	nonbonded/nb_kernel_bluegene/nb_kernel010_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel020_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel030_bluegene.c

	nonbonded/nb_kernel_bluegene/nb_kernel400_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel410_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel420_bluegene.c
	nonbonded/nb_kernel_bluegene/nb_kernel430_bluegene.c

	nonbonded/nb_kernel_bluegene/nb_kernel_bluegene.c
  )
endif(GMX_BLUEGENE)

if(GMX_PPC_ALTIVEC)
  file(GLOB GMX_PPC_ALTIVEC_SRC nonbonded/nb_kernel_ppc_altivec/*.c)
endif(GMX_PPC_ALTIVEC)

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

# only fiddle with assembly kernels if we're not doing OpenMM build
if(NOT GMX_OPENMM) 
if(GMX_ASM_USEASM-NASM)
  enable_language(ASM-NASM)

  # if NASM is used, we need a special build command for windows...
  FOREACH(SRC ${GMX_SSEKERNEL_ASM_SRC})
    GET_FILENAME_COMPONENT(FILE_BASE ${SRC} NAME_WE)
    SET(OBJ ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}${CMAKE_C_OUTPUT_EXTENSION})

    ADD_CUSTOM_COMMAND(OUTPUT ${OBJ}
                       MAIN_DEPENDENCY ${SRC}
                       COMMAND ${CMAKE_ASM-NASM_COMPILER} -f ${CMAKE_ASM-NASM_OBJECT_FORMAT} -o ${OBJ} ${SRC})

    SET(ALL_ASM_OBJS ${ALL_ASM_OBJS} ${OBJ})
  ENDFOREACH(SRC ${GMX_SSEKERNEL_ASM_SRC})
  set(GMX_SSEKERNEL_ASM_SRC ${ALL_ASM_OBJS})

else(GMX_ASM_USEASM-NASM)

  enable_language(ASM-ATT)
  SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
  if(GMX_IA32_ASM)
    set_source_files_properties(${GMX_SSEKERNEL_ASM_SRC} PROPERTIES COMPILE_FLAGS "-c -m32")
  else()
    set_source_files_properties(${GMX_SSEKERNEL_ASM_SRC} PROPERTIES COMPILE_FLAGS "-c -m64")
  endif()

endif(GMX_ASM_USEASM-NASM)
endif(NOT GMX_OPENMM)

add_library(gmx ${GMXLIB_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES} ${GMX_SSEKERNEL_C_SRC} ${GMX_SSEKERNEL_ASM_SRC} ${FORTRAN_SOURCES} ${GMX_BLUEGENE_C_SRC} ${GMX_PPC_ALTIVEC_SRC} ${THREAD_MPI_SRC})
target_link_libraries(gmx ${GMX_EXTRA_LIBRARIES}  ${THREAD_LIB})
if(USE_VERSION_H)
	add_dependencies(gmx gmx_version) 
endif()
set_target_properties(gmx PROPERTIES OUTPUT_NAME "gmx${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

install(TARGETS gmx DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmx.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmx.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmx${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)

