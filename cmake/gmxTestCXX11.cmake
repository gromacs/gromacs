#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

# Check whether both a suitable C++11-compatible compiler and standard
# library is available, and give a fatal error if not.
#
# Any required compiler flag for C++11 support is returned in
# ${FLAG}. The other parameters are only inputs, naming variables that
# contain flags that may have been detected, or set by the user.
function(GMX_TEST_CXX11 CXX11_CXX_FLAG_NAME STDLIB_CXX_FLAG_NAME STDLIB_LIBRARIES_NAME)

    # First check that the compiler is OK, and find the appropriate flag.

    if(WIN32 AND NOT MINGW)
        set(CXX11_CXX_FLAG "/Qstd=c++0x")
    elseif(CYGWIN)
        set(CXX11_CXX_FLAG "-std=gnu++0x") #required for strdup
    else()
        set(CXX11_CXX_FLAG "-std=c++0x")
    endif()
    CHECK_CXX_COMPILER_FLAG("${CXX11_CXX_FLAG}" CXXFLAG_STD_CXX0X)
    if(NOT CXXFLAG_STD_CXX0X)
        set(CXX11_CXX_FLAG "")
    endif()
    set(CMAKE_REQUIRED_FLAGS "${CXX11_CXX_FLAG}")
    check_cxx_source_compiles(
"// Test that a subclass has a proper copy constructor
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
// Early patch versions of icc 16 (and perhaps earlier versions)
// have an issue with this test, but the GROMACS tests pass,
// so we disable this test in that sub-case.
#if (defined __INTEL_COMPILER && __INTEL_COMPILER >= 1700) || (defined __ICL && __ICL >= 1700) || (defined __INTEL_COMPILER_UDPATE && __INTEL_COMPILER_UPDATE >= 3)
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
// Test that operator bool() works
struct e {
  explicit operator bool() {return true;}
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
int main() {
  // Test nullptr
  double *x = nullptr;
  // Test range-based for loops
  int array[5] = { 1, 2, 3, 4, 5 };
  for (int& x : array)
    x *= 2;
}" CXX11_SUPPORTED)
    if(CXX11_SUPPORTED)
        set(${CXX11_CXX_FLAG_NAME} ${CXX11_CXX_FLAG} PARENT_SCOPE)
    else()
        message(FATAL_ERROR "This version of GROMACS requires a C++11 compiler. Please use a newer compiler or use the GROMACS 5.1.x release. See the installation guide for details.")
    endif()

    # Now check the standard library is OK

    set(CMAKE_REQUIRED_FLAGS "${CXX11_CXX_FLAG} ${${STDLIB_CXX_FLAG_NAME}}")
    set(CMAKE_REQUIRED_LIBRARIES "${${STDLIB_LIBRARIES_NAME}}")
    check_cxx_source_compiles(
"#include <chrono>
#include <map>
#include <memory>
#include <thread>
#include <utility>
int main() {
  typedef std::unique_ptr<int> intPointer;
  intPointer p(new int(10));
  std::map<int, std::unique_ptr<int>> m;
  m.insert(std::make_pair(5, std::move(p)));
  auto start = std::chrono::system_clock::now();
  if (std::chrono::system_clock::now() - start < std::chrono::seconds(2))
  {
      std::thread t;
  }
}" CXX11_STDLIB_PRESENT)
    if(NOT CXX11_STDLIB_PRESENT)
        message(FATAL_ERROR "This version of GROMACS requires C++11-compatible standard library. Please use a newer compiler, or a newer standard library, or use the GROMACS 5.1.x release. See the installation guide for details.")
    endif()
endfunction()
