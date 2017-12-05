#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

include(CheckCXXSourceCompiles)
include(FindThreads)

# Check whether both a suitable C++11-compatible compiler and standard
# library is available, and give a fatal error if not.
#
# Any required compiler flag for C++11 support is returned in
# ${FLAG}. The other parameters are only inputs, naming variables that
# contain flags that may have been detected, or set by the user.
function(GMX_TEST_CXX11 CXX11_CXX_FLAG_NAME STDLIB_CXX_FLAG_NAME STDLIB_LIBRARIES_NAME)

    # First check that the compiler is OK, and find the appropriate flag.

    if(WIN32 AND NOT MINGW)
        set(CXX11_CXX_FLAG "/Qstd=c++11")
    elseif(CYGWIN)
        set(CXX11_CXX_FLAG "-std=gnu++11") #required for strdup
    else()
        set(CXX11_CXX_FLAG "-std=c++11")
    endif()
    CHECK_CXX_COMPILER_FLAG("${CXX11_CXX_FLAG}" CXXFLAG_STD_CXX0X)
    if(NOT CXXFLAG_STD_CXX0X)
        set(CXX11_CXX_FLAG "")
    endif()
    set(CMAKE_REQUIRED_FLAGS "${CXX11_CXX_FLAG}")
    check_cxx_source_compiles(
"// Permit testing typeid keyword
#include <typeinfo>
// Test that a subclass has a proper copy constructor
struct a {
  a() {};
  a(const a&) {};
  a(a&&) = delete;
};
class b: public a
{
};
b bTest() {
  return b();
}
// ICC requires that a suitable GCC is available. It is using its standard library and emulates
// GCC behaviour based on its version. Relevant here it emulates the implementation of the move
// constructor. This compiler check should only fail based on the compiler not GCC. The GCC version
// is checked by the following STL check. It is known that all ICC>=15 have the proper move
// constructor. Thus this check is disabled for ICC.
#if !((defined __INTEL_COMPILER && __INTEL_COMPILER >= 1500) || (defined __ICL && __ICL >= 1500))
// Test that a subclass has a proper move constructor
struct c {
  c() {};
  c(const c&) = delete;
  c(c&&) {};
};
struct d : public c {
};
d dTest() {
  return d();
}
#endif
struct e {
  // Test that operator bool() works
  explicit operator bool() {return true;}
  // Test that an in-class initializer works
  int x = 1;
  // Test that a default constructor is generated
  e() = default;
};
// Test that constexpr works
constexpr int factorial(int n)
{
    return n <= 1? 1 : (n * factorial(n - 1));
}
// Test that r-value references work
void checkRvalueReference(int &&);
// Test that extern templates work
template <typename T> void someFunction();
extern template void someFunction<int>();
// Test using statement
using myInt = int;
// Test template using statement
template<class T> using myPointer = T*;
myPointer<int> x;
int main() {
  // Test nullptr
  double *x = nullptr;
  (void)x; // Suppressing unused variable warning
  // Test range-based for loops
  int array[5] = { 1, 2, 3, 4, 5 };
  for (int& x : array)
    x *= 2;
  // Test alignas
  alignas(4*sizeof(int)) int y;
  // Test typeid
  const std::type_info &intType = typeid(int);
  // Test static assertions do compile
  static_assert(true, \"if you see this, true somehow isn't\");
  // Test a lambda
  [=]{};
}" CXX11_SUPPORTED)
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.1")
            message(FATAL_ERROR "GROMACS requires version 4.8.1 or later of the GNU C++ compiler for complete C++11 support")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
            message(FATAL_ERROR "GROMACS requires version 3.3 or later of the Clang C++ compiler for complete C++11 support")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0")
            message(FATAL_ERROR "GROMACS requires version 15.0 or later of the Intel C++ compiler for complete C++11 support")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.0.23026")
            message(FATAL_ERROR "GROMACS requires version 2015 (19.0.23026) or later of the MSVC C++ compiler for complete C++11 support")
        endif()
    endif()
    if(CXX11_SUPPORTED)
        set(${CXX11_CXX_FLAG_NAME} ${CXX11_CXX_FLAG} PARENT_SCOPE)
    else()
        message(FATAL_ERROR "This version of GROMACS requires a C++11 compiler. Please use a newer compiler or use the GROMACS 5.1.x release. See the installation guide for details.")
    endif()

    # Now check the standard library is OK

    set(CMAKE_REQUIRED_FLAGS "${CXX11_CXX_FLAG} ${${STDLIB_CXX_FLAG_NAME}}")
    set(CMAKE_REQUIRED_LIBRARIES "${${STDLIB_LIBRARIES_NAME}} ${CMAKE_THREAD_LIBS_INIT}")
    check_cxx_source_compiles(
"#include <algorithm>
#include <array>
#include <chrono>
#include <iterator>
#include <map>
#include <memory>
#include <thread>
#include <type_traits>
#include <tuple>
#include <utility>
#include <string>
#include <vector>
int main() {
  // Test for std::vector
  std::vector<double> doubles(100);
  // Test for std::array
  std::array<int, 3> someInts;
  // Test std::for_each and a lambda
  std::for_each(std::begin(doubles), std::end(doubles), [&](double &d) { d = 2.3; });
  // Test std::unique_ptr
  typedef std::unique_ptr<int> intPointer;
  // Test using std::unique_ptr
  intPointer p(new int(10));
  // Test std::map
  std::map<int, std::unique_ptr<int>> m;
  // Test std::make_pair
  m.insert(std::make_pair(5, std::move(p)));
  // Test std::chrono (was missing before gcc 4.8.1)
  auto start = std::chrono::steady_clock::now();
  if (std::chrono::steady_clock::now() - start < std::chrono::seconds(2))
  {
      // Test std::thread
      std::thread t;
  }
  // Test std::is_pod
  static_assert(std::is_pod<int>::value, \"int isn't pod\");
  // Test std::tuple
  auto theTuple = std::make_tuple<int, double>(3, 4.2);
  // Test std::tie
  int tupleInt;
  double tupleDouble;
  std::tie(tupleInt, tupleDouble) = theTuple;
  // Test std::string
  std::string message(\"hello\");
}" CXX11_STDLIB_PRESENT)
    if(NOT CXX11_STDLIB_PRESENT)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0.3")
            message(FATAL_ERROR "GROMACS requires that the Intel C++ compiler use a compatible C++ standard library, for complete C++11 support, however not all compiler versions support all possible standard library implementations. In particular, before icc version 16.0.3, the gcc version 5 standard library was not supported. If you are affected by such a case you should probably update your compiler version. Consult the GROMACS installation guide to check exactly what is supported and how to direct the use of a standard library from an older gcc version (but at least version 4.8.1 is needed).")
        else()
            message(FATAL_ERROR "This version of GROMACS requires C++11-compatible standard library. Please use a newer compiler, and/or a newer standard library, or use the GROMACS 5.1.x release. Consult the installation guide for details before upgrading components.")
        endif()
    endif()
endfunction()
