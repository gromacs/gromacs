#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

# - Define macro to check if compiling and linking to an external
# Boost actually works, because the find_package macro is content if
# the headers are found. This can fail (for example) with a C++11
# compiler (like clang 3.4) when CMake picks up Boost intended to be
# used with gcc 4.8. This fails messily and with highly cryptic error
# messages, probably for reasons that are related to the language
# changes.
#
#  GMX_FIND_BOOST(VARIABLE)
#
#  VARIABLE will be set to true if Boost support is present

include(CheckCSourceCompiles)
include(gmxOptionUtilities)
function(GMX_FIND_BOOST VARIABLE USE_CXX11 CXX11_FLAGS)
    set(minimum_version 1.44.0)
    if (NOT DEFINED GMX_EXTERNAL_BOOST OR GMX_EXTERNAL_BOOST)
        find_package(Boost ${minimum_version})
        if(Boost_FOUND AND Boost_VERSION VERSION_LESS "104400")
            set(Boost_FOUND FALSE)
        endif()
    endif()

    if(GMX_EXTERNAL_BOOST AND NOT Boost_FOUND)
        message(FATAL_ERROR "Boost >= ${minimum_version} not found. You can set GMX_EXTERNAL_BOOST=OFF to compile against the minimal version of Boost included with GROMACS.")
    endif()

    if(Boost_FOUND)
        gmx_check_if_changed(_do_boost_recompile BOOST_INCLUDE_DIRS)
        if(_do_boost_recompile)
            unset(BOOST_COMPILES_OK CACHE)
        endif()
        if(USE_CXX11)
            # This is important because Boost implementations are not
            # necessarily portable across C++ language versions. We
            # need to do the compiler test the same way we will use
            # it.
            set(CMAKE_REQUIRED_FLAGS "${CXX11_FLAGS}")
        endif()
        set(CMAKE_REQUIRED_INCLUDES "${BOOST_INCLUDE_DIRS}")
        check_cxx_source_compiles(
"/* Include all Boost headers that are used by GROMACS, e.g. with
 * git grep -h 'include.*boost' -- src/gromacs src/programs | sort | uniq */
#include <boost/current_function.hpp>
#include <boost/exception/detail/attribute_noreturn.hpp>
#include <boost/exception/errinfo_api_function.hpp>
#include <boost/exception/errinfo_errno.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/get_error_info.hpp>
#include <boost/exception/info.hpp>
#include <boost/exception_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/throw_exception.hpp>
int main(void) { return 0; }"
            BOOST_COMPILES_OK)

        if(GMX_EXTERNAL_BOOST AND NOT BOOST_COMPILES_OK)
            message(FATAL_ERROR "Boost >= ${minimum_version} found, but a test program could not be compiled and linked with it. You can set GMX_EXTERNAL_BOOST=OFF to compile against the minimal version of Boost included with GROMACS.")
        endif()

        set(${VARIABLE} ${BOOST_COMPILES_OK} PARENT_SCOPE)
    else()
        set(${VARIABLE} OFF PARENT_SCOPE)
    endif()

endfunction()



