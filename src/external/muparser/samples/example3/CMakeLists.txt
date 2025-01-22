if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.28)
endif()
project(muparser-example3)

# find muparser target already installed
find_package(muparser 2.0 REQUIRED)

add_executable(example3 example3.cpp)
target_link_libraries(example3 muparser::muparser)

include(CTest)
add_test(example3 example3)

