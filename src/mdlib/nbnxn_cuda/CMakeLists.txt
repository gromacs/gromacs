if(GMX_GPU)
    remove_definitions("-DMPI")

    file(GLOB CUDA_NB_SOURCES *.cu)
    CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/gmxlib/cuda_tools)
    CUDA_ADD_LIBRARY(nbnxn_cuda STATIC ${CUDA_NB_SOURCES}
            OPTIONS
            RELWITHDEBINFO -g
            DEBUG -g -D_DEBUG_=1)
endif()
