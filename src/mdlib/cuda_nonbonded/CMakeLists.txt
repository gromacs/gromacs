if(GMX_GPU)
    remove_definitions("-DMPI")

    file(GLOB CUDA_NB_SOURCES *.cu)
    CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} 
                            "${CMAKE_CURRENT_SOURCE_DIR}/../cuda_tools"
                            "${CMAKE_SOURCE_DIR}/include/types")
    CUDA_ADD_LIBRARY(cuda_nonbonded STATIC ${CUDA_NB_SOURCES}
            OPTIONS
            RELWITHDEBINFO -g
            DEBUG -g -D_DEBUG_=1)
endif()
