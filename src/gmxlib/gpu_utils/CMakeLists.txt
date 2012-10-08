# (slightly sloppy) OS definitions required by memtestG80
set(_os_def)
if(UNIX)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(_os_def "-DOSX")
    else() # everything that's UNIX & UNIX-like except OS X
        set(_os_def "-DLINUX")
    endif()
else()
    if(WIN32)
        set(_os_def "-DWINDOWS")
    else()
        message(FATAL_ERROR " Could not detect OS required for memtestG80.")
    endif()
endif()

CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)        
file(GLOB GPU_UTILS_SOURCES *.cu dummy.cpp)
CUDA_ADD_LIBRARY(gpu_utils STATIC ${GPU_UTILS_SOURCES}
                 OPTIONS ${_os_def}
                 RELWITHDEBINFO -g
                 DEBUG -g -D_DEBUG_=1 )

CUDA_BUILD_CLEAN_TARGET()
