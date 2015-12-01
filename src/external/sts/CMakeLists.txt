cmake_minimum_required(VERSION 2.8.12)
project(sts)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package( Threads )

set(SOURCE_FILES ex2.cpp sts/sts.cpp sts/range.cpp sts/thread.cpp)
add_executable(sts ${SOURCE_FILES})
target_link_libraries( sts ${CMAKE_THREAD_LIBS_INIT} )