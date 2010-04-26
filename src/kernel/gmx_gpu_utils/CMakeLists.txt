

CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)        
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    CUDA_ADD_LIBRARY(gmx_gpu_utils  
        gmx_gpu_utils.cu 
        memtestG80_core.cu 
        OPTIONS "-DLINUX"
        DEBUG -g -D_DEBUG_=1 
        )
else()
    CUDA_ADD_LIBRARY(gmx_gpu_utils 
        gmx_gpu_utils.cu 
        memtestG80_core.cu 
        OPTIONS "-DLINUX"
        )
endif()
CUDA_BUILD_CLEAN_TARGET()

install(TARGETS gmx_gpu_utils DESTINATION ${LIB_INSTALL_DIR})   

