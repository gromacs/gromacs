
file(GLOB MDLIB_SOURCES *.c)

if(GMX_GPU)
  file(GLOB MDLIB_GPU_SOURCES *.cu)
  CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)        
    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
        CUDA_ADD_LIBRARY(md_gpu 
            ${MDLIB_GPU_SOURCES}
            DEBUG -g -D_DEBUG_=1)
    else()
        CUDA_ADD_LIBRARY(md_gpu 
            ${MDLIB_GPU_SOURCES})
    endif()

    # add_subdirectory("../kernel/gmx_gpu_utils")   
    #include_directories("${CMAKE_SOURCE_DIR}/src/kernel/gmx_gpu_utils")
    set(GMX_GPU_LIBRARIES gmx_gpu_utils)
endif()

# Files	called xxx_test.c are test drivers with a main() function for 
# module xxx.c, so they should not be included in the library
file(GLOB_RECURSE NOT_MDLIB_SOURCES *_test.c *\#*)
list(REMOVE_ITEM MDLIB_SOURCES ${NOT_MDLIB_SOURCES})

add_library(md ${MDLIB_SOURCES})

if(GMX_GPU)
    target_link_libraries(md md_gpu gmx ${GMX_GPU_LIBRARIES} ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES})
else()
    target_link_libraries(md gmx ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES})
endif()

set_target_properties(md PROPERTIES OUTPUT_NAME "md${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION})

install(TARGETS md DESTINATION ${LIB_INSTALL_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmd.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmd.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig RENAME "libmd${GMX_LIBS_SUFFIX}.pc")
