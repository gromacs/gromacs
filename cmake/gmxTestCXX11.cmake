include(CheckCXXSourceCompiles)
MACRO(GMX_TEST_CXX11 VARIABLE FLAG)
    IF(NOT DEFINED HAVE_${VARIABLE})
        MESSAGE(STATUS "Checking for C++11 support")
        if(NOT WIN32)
            set(CXX11_FLAG "-std=c++0x")
        else()
            set(CXX11_FLAG "/Qstd=c++0x")
        endif()
        CHECK_CXX_COMPILER_FLAG("${CXX11_FLAG}" CXXFLAG_STD_CXX0X)
        if(NOT CXXFLAG_STD_CXX0X)
            set(CXX11_FLAG "")
        endif()
        set(CMAKE_REQUIRED_DEFINITIONS "${CXX11_FLAG}")
        check_cxx_source_compiles(
"#include <vector>
#include <memory>
#include <utility>
struct A {
  std::unique_ptr<int> p;
};
int main() {
  typedef std::unique_ptr<int> intPointer;
  intPointer p(new int(10));
  std::vector<intPointer> v;
  v.push_back(std::move(p));
  std::vector<A> v2;
  v2.push_back(A());  //requires default move constructor
}" HAVE_${VARIABLE})
        set(CMAKE_REQUIRED_DEFINITIONS "")
        if(HAVE_${VARIABLE})
            set(${VARIABLE} 1 CACHE INTERNAL "Result of C++11 support test" FORCE)
            set(${FLAG} ${CXX11_FLAG} CACHE INTERNAL "Compiler flag for C++11 support" FORCE)
            MESSAGE(STATUS "Checking for C++11 support - yes")
        else()
            set(${VARIABLE} 0 CACHE INTERNAL "Result of C++11 support test" FORCE)
            set(${FLAG} "" CACHE INTERNAL "Compiler flag for C++11 support" FORCE)
            MESSAGE(STATUS "Checking for C++11 support - no")
        endif()
    ENDIF(NOT DEFINED HAVE_${VARIABLE})
ENDMACRO()
