#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# For compilers which might require libstdc++ (Clang and
# vendor-specific forks thereof, such as Intel and Fujitsu), and do
# not already work, find it and set CMAKE_CXX_FLAGS.
#
# Does nothing if compiler includes std-library (e.g. GCC), or already
# works, or compiler uses different std-library
# (either because of different defaults (e.g. on MacOS) or user flags (e.g. -stdlib=libc++)).
# The heuristic by the compiler of how to find libstdc++ is honored. Any user-provided flags in
# e.g. CXXFLAGS for the location of libstdc++ are honored. The user can choose the libstdc++ by setting
# GMX_GPLUSPLUS_PATH, PATH or CMAKE_PREFIX_PATH to make sure the correct the g++ is found.
# Gives error if no g++ is found or warns if the g++ found isn't new enough.
# The location of g++ is cached as GMX_GPLUSPLUS_PATH making sure that the same libstdc++ is used
# for builds at different times using the same cache file (so that e.g. loading a module is
# not required for a reproducible build). Note that GMX_GPLUSPLUS_PATH is ignored if it is
# not needed because the compiler already found a std library via some other mechanism.

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel") # Compilers supported
    return()
endif()

include(CheckCXXSourceCompiles)

# Test that required 2017 standard library features work.
# Note that this check also requires linking to succeed, and
# GROMACS does not support using std::experimental::filesystem.
set (SAMPLE_CODE_TO_TEST_CXX17_WITHOUT_STD_FILESYSTEM "
#include <string>
#include <string_view>
#include <optional>
int main(int argc, char **argv) {
  std::optional<std::string> input(argv[0]);
  std::string_view view(input.value());
  return int(view[0]);
}")
set (SAMPLE_CODE_TO_TEST_CXX17 "
#include <filesystem>
${SAMPLE_CODE_TO_TEST_CXX17_WITHOUT_STD_FILESYSTEM}
size_t foo()
{
  auto path = std::filesystem::current_path();
  return path.string().size();
}")
check_cxx_source_compiles("${SAMPLE_CODE_TO_TEST_CXX17}" CXX17_COMPILES_SIMPLY)

if (CXX17_COMPILES_SIMPLY)
    # The compiler has been set up properly to find a standard
    # library, and if so GROMACS should leave it alone.
    return()
endif()

# The compiler couldn't use the standard libary for an unknown reason.
# See if the compiler is using libstdc++ (via libstc++ heuristics). If so,
# then we may be able to help the compiler find the standard library.
check_cxx_source_compiles("#include <new>
int main() { return __GLIBCXX__; }" USING_LIBSTDCXX)

if (USING_LIBSTDCXX)
    # If the existing compiler options already set a gcc toolchain then we should not
    # attempt to add that, regardless of whether the user set GMX_GPLUSPLUS.
    if (DEFINED GMX_GPLUSGPLUS_PATH)
        set(EXTRA_MESSAGE ", ignoring the value of GMX_GPLUSPLUS_PATH")
    endif()
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if ("${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_cmake_build_type}}" MATCHES "--gcc-toolchain")
        message(STATUS "The --gcc-toolchain option is already present in the CMAKE_CXX_FLAGS "
            "(or perhaps those specific to the CMAKE_BUILD_TYPE), and the GROMACS build "
            "will use that one${EXTRA_MESSAGE}.")
    else()
        set(TRY_TO_FIND_GPLUSPLUS TRUE)
    endif()
endif()

if (TRY_TO_FIND_GPLUSPLUS)
    # Find a gcc (perhaps already specified by the user in
    # GMX_GPLUSPLUS_PATH) and prepare to reproducibly use its libstdc++.
    find_program(GMX_GPLUSPLUS_PATH g++)
    if (NOT EXISTS "${GMX_GPLUSPLUS_PATH}")
        message(FATAL_ERROR "Couldn't find g++. Please set GMX_GPLUSPLUS_PATH, PATH or CMAKE_PREFIX_PATH "
            "accordingly for cmake to find it. Note that clang supports using libc++ with "
            "-DCMAKE_CXX_FLAGS=--stdlib=libc++, and if so there will be no need to find g++.")
    endif()

    # Ensure that a suitable version of g++ was found, caching the
    # result for future configuration stages.
    # Note that this approach does not work when cross compiling.
    if (NOT GMX_GPLUSPLUS_VERSION)
        execute_process(COMMAND ${GMX_GPLUSPLUS_PATH} -dumpfullversion -dumpversion OUTPUT_VARIABLE GMX_GPLUSPLUS_VERSION
            ERROR_VARIABLE GMX_GPLUSPLUS_VERSION_ERROR
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (NOT "${GMX_GPLUSPLUS_VERSION}" MATCHES "^[0-9]+\\.[0-9]+\\.?[0-9]?$") #Should never happen
            message(FATAL_ERROR "Couldn't detect g++ version for ${GMX_GPLUSPLUS_PATH}. Version output: ${GMX_GPLUSPLUS_VERSION} "
                ", error: ${GMX_GPLUSPLUS_VERSION_ERROR}. Please report to developers.")
        endif()
        # Cache this, so future configurations won't have to run g++ again.
        set(GMX_GPLUSPLUS_VERSION ${GMX_GPLUSPLUS_VERSION} CACHE STRING "Version of g++ from which libstdc++ is obtained")
    endif()
    if (${GMX_GPLUSPLUS_VERSION} VERSION_LESS 9)
        message(WARNING "Found g++ at ${GMX_GPLUSPLUS_PATH}. Its version is ${GMX_GPLUSPLUS_VERSION}. "
            "GROMACS encourages at least version 9. "
            "If you see problems, please specify a different g++ using GMX_GPLUSPLUS_PATH, PATH or CMAKE_PREFIX_PATH.")
    endif()

    # Now make some sanity checks on the compiler using libstdc++.
    # We assume the compiler is based on clang and thus it will take a --gcc-toolchain option, which is checked later.
    get_filename_component(_gplusplus_path "${GMX_GPLUSPLUS_PATH}" REALPATH)
    get_filename_component(_gplusplus_path "${_gplusplus_path}" DIRECTORY) #strip g++
    get_filename_component(_gplusplus_path "${_gplusplus_path}" DIRECTORY) #strip bin
    set(GMX_GPLUSPLUS_INCLUDE_PATH_FOUND FALSE)
    # On some Cray systems, we have a symlink /foo/bin/g++, while the real installation is in /foo/snos/
    foreach (INFIX "/" "/snos/")
        foreach (INCLUDE_DIR in "include/c++" "include/g++")
            if (EXISTS "${_gplusplus_path}${INFIX}${INCLUDE_DIR}")
                set(GMX_GPLUSPLUS_INCLUDE_PATH_FOUND TRUE)
                break()
            endif()
        endforeach()
        if(GMX_GPLUSPLUS_INCLUDE_PATH_FOUND)
            set(_gplusplus_path "${_gplusplus_path}${INFIX}")
            break()
        endif()
    endforeach()
    if(NOT GMX_GPLUSPLUS_INCLUDE_PATH_FOUND)
        message(FATAL_ERROR "Directory ${_gplusplus_path}/include/c++ found via ${GMX_GPLUSPLUS_PATH} "
            "doesn't exist even though it should. Please report this to the GROMACS developers.")
    endif()

    set(_gcc_toolchain_flag_to_use "--gcc-toolchain=${_gplusplus_path}")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(${_gcc_toolchain_flag_to_use} CXX_COMPILER_ACCEPTS_GCC_TOOLCHAIN_FLAG)

    if (CXX_COMPILER_ACCEPTS_GCC_TOOLCHAIN_FLAG)
        # Set up to use the libstdc++ from that g++. Note that we checked
        # the existing contents of CMAKE_CXX_FLAGS* variables earlier, so
        # we will not override any user settings here.
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_gcc_toolchain_flag_to_use}")
    else()
        message(STATUS "Your compiler is using the libstdc++ C++ standard library but does not "
            "recognize the --gcc-toolchain flag. The GROMACS build system does not know how to "
            "help find a standard library for this compiler. If you experience further problems "
            "related to this, please report this to the GROMACS developers.")
    endif()
endif()

if (NOT USING_LIBSTDCXX)
    # Just linking an extra library for std::filesystem didn't help,
    # so let's try to narrow down what fails.
    check_cxx_source_compiles("${SAMPLE_CODE_TO_TEST_CXX17_WITHOUT_STD_FILESYSTEM}" CXX17_COMPILES_WITHOUT_STD_FILESYSTEM)
    if (NOT CXX17_COMPILES_WITHOUT_STD_FILESYSTEM)
        message(FATAL_ERROR "The C++ compiler cannot find a working standard library. "
            "The compiler was not trying to use libstdc++. "
            "The GROMACS build system cannot handle this case. "
            "Please use a working C++17 compiler and standard library.")
    else()
        message(FATAL_ERROR "The C++ compiler cannot find a working standard library "
            "that supports std::filesystem. "
            "The compiler was not trying to use libstdc++. "
            "The GROMACS build system cannot handle this case. "
            "Please use a working C++17 compiler and standard library.")
    endif()
endif()

# Now run a sanity check on the compiler using libstdc++, regardless
# of how it was specified or found.
check_cxx_source_compiles("${SAMPLE_CODE_TO_TEST_CXX17}" CXX17_COMPILES_WITH_HELP)

if (NOT CXX17_COMPILES_WITH_HELP)
    if (TRY_TO_FIND_GPLUSPLUS)
        set (OTHER_EXTRA_MESSAGE " The g++ found at ${_gplusplus_path} found via ${GMX_GPLUSPLUS_PATH} had a suitable version, so something else must be the problem.")
    else()
        set (OTHER_EXTRA_MESSAGE " Check your toolchain documentation or environment flags so that they will find a suitable C++17 standard library.")
    endif()
    message(FATAL_ERROR "GROMACS requires C++17, but a test of such functionality in the C++ standard "
        "library failed to compile.${OTHER_EXTRA_MESSAGE}")
endif()
