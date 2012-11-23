
file(GLOB MDLIB_SOURCES nbnxn_kernels/*.c *.c *.cpp)
if(GMX_FFT_FFTPACK)
list(APPEND MDLIB_SOURCES ${CMAKE_SOURCE_DIR}/src/external/fftpack/fftpack.c)
endif()
set(MDLIB_SOURCES ${MDLIB_SOURCES} PARENT_SCOPE)
if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)

if(GMX_GPU)
    add_subdirectory(nbnxn_cuda)
    set(GMX_GPU_LIBRARIES ${GMX_GPU_LIBRARIES} nbnxn_cuda PARENT_SCOPE)
endif()
