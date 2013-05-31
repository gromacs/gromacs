
file(GLOB MDLIB_SOURCES nbnxn_kernels/*.c *.c *.cpp)
set(MDLIB_SOURCES ${MDLIB_SOURCES} PARENT_SCOPE)

if(GMX_GPU)
    add_subdirectory(nbnxn_cuda)
    set(GMX_GPU_LIBRARIES ${GMX_GPU_LIBRARIES} nbnxn_cuda PARENT_SCOPE)
endif()
