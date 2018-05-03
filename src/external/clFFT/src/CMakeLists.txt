project( clFFT )

# get sources
file(GLOB CLFFT_SOURCES statTimer/*.cpp library/*.cpp library/*.c)
file(GLOB CLFFT_SOURCES_UNWANTED statTimer/dllmain.cpp library/dllmain.cpp)
list(REMOVE_ITEM CLFFT_SOURCES ${CLFFT_SOURCES_UNWANTED})
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/include)

# disable all warnings
set_source_files_properties(${CLFFT_SOURCES} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -w")

# link
add_library(clFFT STATIC ${CLFFT_SOURCES})
set_property(TARGET clFFT PROPERTY POSITION_INDEPENDENT_CODE ON)
