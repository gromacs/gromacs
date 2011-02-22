if(GMX_GPU)
    file(GLOB CUDA_NB_SOURCES *.cu)
    CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../cuda_tools")
    set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)        
    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
        CUDA_ADD_LIBRARY(cuda_nonbonded STATIC
            ${CUDA_NB_SOURCES}
            DEBUG -g -D_DEBUG_=1)
    else()
        CUDA_ADD_LIBRARY(cuda_nonbonded STATIC
            ${CUDA_NB_SOURCES})
    endif()
endif()
