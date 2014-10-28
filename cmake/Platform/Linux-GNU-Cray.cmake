
#=============================================================================
# Copyright 2010 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# This version is taken from cmake 3.0.2 Modules/Platform/Linux-GNU.cmake

# This module is shared by multiple languages; use include blocker.
if(__LINUX_COMPILER_GNU)
  return()
endif()
set(__LINUX_COMPILER_GNU 1)

macro(__linux_compiler_gnu lang)
  # Gromacs modification here, so that compiler tests are not forced
  # to use -rdynamic, which is not the intended use of the Cray
  # toolchain, even though the underlying compiler might support it.
  set(CMAKE_SHARED_LIBRARY_LINK_${lang}_FLAGS "")
endmacro()
