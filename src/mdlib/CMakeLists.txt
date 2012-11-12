# The nbnxn kernel files take the longest time to compile, so we sneak
# them in first to take best advantage of make-time parallelisation.
file(GLOB MDLIB_SOURCES nbnxn_kernels/*.c *.c)

if(GMX_GPU)
    add_subdirectory(nbnxn_cuda)
    set(GMX_GPU_LIBRARIES nbnxn_cuda)
endif()

# apply gcc 4.4.x bug workaround
if(GMX_USE_GCC44_BUG_WORKAROUND)
   include(gmxGCC44O3BugWorkaround)
   gmx_apply_gcc44_bug_workaround("force.c")
   gmx_apply_gcc44_bug_workaround("constr.c")
endif()

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library

if(NOT GMX_FFT_FFTPACK)
list(REMOVE_ITEM MDLIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/fftpack.c)
endif()

add_library(md ${MDLIB_SOURCES})

target_link_libraries(md ${GMX_GPU_LIBRARIES} gmx ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${OpenMP_SHARED_LINKER_FLAGS})
if(GMX_BUILD_OWN_FFTW)
    # This dependency has to be made here rather than the CMakeLists.txt that
    # does the FFTW build, because of the order in which
    # add_subdirectory() calls are made in the top-level CMakeLists.txt; the
    # md library target does not necessarily exist yet. Also enabling and
    # disabling GMX_BUILD_OWN_FFTW changes dependencies correctly.
    add_dependencies(md gmxfftw)
endif()
set_target_properties(md PROPERTIES OUTPUT_NAME "md${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}"
    COMPILE_FLAGS "${OpenMP_C_FLAGS}")

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmd.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libmd${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
