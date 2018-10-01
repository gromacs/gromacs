project( clFFT )

# get sources
file(GLOB CLFFT_SOURCES statTimer/*.cpp library/*.cpp library/*.c)
file(GLOB CLFFT_SOURCES_UNWANTED statTimer/dllmain.cpp library/dllmain.cpp)
list(REMOVE_ITEM CLFFT_SOURCES ${CLFFT_SOURCES_UNWANTED})

add_library(clFFT OBJECT ${CLFFT_SOURCES})

# disable all warnings
target_compile_options(clFFT PRIVATE "-w")
if (BUILD_SHARED_LIBS)
    set_target_properties(clFFT PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        )
endif()

# Windows doesn't need anything special to load the dynamic libraries
# that the AMD clFFT implementation uses, but *nix and BSD do.
if (NOT WIN32)
    include(CheckCXXSymbolExists)
    check_cxx_symbol_exists(dlopen dlfcn.h HAVE_DLOPEN)
    if (NOT HAVE_DLOPEN)
        message(FATAL_ERROR "Cannot use clFFT for OpenCL support unless dlopen is available")
    endif()
endif()

target_include_directories(clFFT
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/library
        ${CMAKE_CURRENT_SOURCE_DIR}/statTimer
        )
