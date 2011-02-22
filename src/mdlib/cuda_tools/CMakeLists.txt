if(GMX_GPU)
    file(GLOB CUDA_TOOLS_SOURCES *.cu)
    CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)        
    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
        CUDA_ADD_LIBRARY(cuda_tools STATIC
            ${CUDA_TOOLS_SOURCES}
            DEBUG -g -D_DEBUG_=1)
    else()
        CUDA_ADD_LIBRARY(cuda_tools STATIC
            ${CUDA_TOOLS_SOURCES})
    endif()
endif()
