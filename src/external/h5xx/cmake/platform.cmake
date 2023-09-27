set(CMAKE_BUILD_TYPE_INIT "Release")

if(DEFINED CMAKE_CXX_COMPILER_ID)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++11 -pedantic -Wno-long-long")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++98 -pedantic -Wno-long-long -Wno-c99-extensions")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++98 -pedantic -Wno-long-long")

  else()

    message(WARNING "Unsupported C++ compiler: ${CMAKE_CXX_COMPILER_ID}")

  endif()
endif()
