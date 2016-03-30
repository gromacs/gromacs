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
MACRO(GMX_TEST_CXX11 VARIABLE FLAG)
    if(WIN32 AND NOT MINGW)
        set(CXX11_FLAG "/Qstd=c++0x")
    elseif(CYGWIN)
        set(CXX11_FLAG "-std=gnu++0x") #required for strdup
    else()
        set(CXX11_FLAG "-std=c++0x")
    endif()
    CHECK_CXX_COMPILER_FLAG("${CXX11_FLAG}" CXXFLAG_STD_CXX0X)
    if(NOT CXXFLAG_STD_CXX0X)
        set(CXX11_FLAG "")
    endif()
    if (GMX_STDLIB_CXX_FLAGS)
        CHECK_CXX_COMPILER_FLAG("${GMX_STDLIB_CXX_FLAGS}" CXXFLAG_STDLIB)
    else()
        set(CXXFLAG_STDLIB "")
    endif()
    set(CMAKE_REQUIRED_FLAGS "${CXX11_FLAG} ${CXXFLAG_STDLIB}")
    set(CMAKE_REQUIRED_LIBRARIES "${GMX_STDLIB_LIBRARIES}")
    check_cxx_source_compiles(
"#include <map>
#include <memory>
#include <utility>
class a { explicit operator bool() {return true;} };
int main() {
  typedef std::unique_ptr<int> intPointer;
  intPointer p(new int(10));
  std::map<int, std::unique_ptr<int>> m;
  m.insert(std::make_pair(5, std::move(p)));
}" ${VARIABLE})
    set(CMAKE_REQUIRED_FLAGS "")
    if(${VARIABLE})
        set(${FLAG} ${CXX11_FLAG})
    endif()
ENDMACRO()
