
file(GLOB MDLIB_SOURCES *.c *.cpp nbnxn_kernels/*.c)
if(GMX_FFT_FFTPACK)
list(APPEND MDLIB_SOURCES ${CMAKE_SOURCE_DIR}/src/external/fftpack/fftpack.c)
endif()
set(MDLIB_SOURCES ${MDLIB_SOURCES} PARENT_SCOPE)
if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)

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
