if(GMX_GPU)
    file(GLOB CUDA_TOOLS_SOURCES *.cu)
    CUDA_ADD_LIBRARY(cuda_tools STATIC ${CUDA_TOOLS_SOURCES}
            OPTIONS
            RELWITHDEBINFO -g
            DEBUG -g -D_DEBUG_=1)
endif()
