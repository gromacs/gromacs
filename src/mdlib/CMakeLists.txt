
file(GLOB MDLIB_SOURCES *.c nbnxn_kernels/*.c)

include_directories(nbnxn_kernels)

if(GMX_GPU)
    add_subdirectory(nbnxn_cuda)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/gmxlib/cuda_tools nbnxn_cuda)
    set(GMX_GPU_LIBRARIES nbnxn_cuda)
endif()

# due to a compiler bug gcc 4.4.x crashses on force.c and constr.c with -O3,
# but strangely not with -O2 + all extra -O3 flags
if(GMX_OPENMP AND CMAKE_COMPILER_IS_GNUCC AND
   C_COMPILER_VERSION VERSION_GREATER "4.3.999" AND C_COMPILER_VERSION VERSION_LESS "4.5.999")
    if(GMX_NO_GCC44_BUG_WORKAROUND)
        message(WARNING "gcc ${C_COMPILER_VERSION} detected, but workaround for optimization bug is disabled")
    else()
        message(STATUS "gcc ${C_COMPILER_VERSION} detected, applying workaround for optimization bug (disable with GMX_NO_GCC44_BUG_WORKAROUND)")
        # O3 only in release mode
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            set_source_files_properties("force.c" PROPERTIES COMPILE_FLAGS
                "-O2 -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -ftree-vectorize -fipa-cp-clone")
            set_source_files_properties("constr.c" PROPERTIES COMPILE_FLAGS
                "-O2 -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -ftree-vectorize -fipa-cp-clone")
        endif()
    endif()
endif()

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library

add_library(md ${MDLIB_SOURCES})

target_link_libraries(md ${GMX_GPU_LIBRARIES} gmx ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES})

set_target_properties(md PROPERTIES OUTPUT_NAME "md${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmd.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libmd${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
