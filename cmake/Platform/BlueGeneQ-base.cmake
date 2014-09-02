#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#
#=============================================================================
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2011 Kitware, Inc., Insight Software Consortium
# Copyright 2010 Todd Gamblin <tgamblin@llnl.gov>
# Copyright 2012 Julien Bigot <julien.bigot@cea.fr>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

#
# BlueGeneQ base platform file.
#
# NOTE: Do not set your platform to "BlueGeneQ-base".  This file is included
# by the real platform files.  Use one of these two platforms instead:
#
#     BlueGeneQ-dynamic  For dynamically linked builds
#     BlueGeneQ-static   For statically linked builds
#
# This platform file tries its best to adhere to the behavior of the MPI
# compiler wrappers included with the latest BG/Q drivers.
#


#
# For BG/Q builds, we're cross compiling, but we don't want to re-root things
# (e.g. with CMAKE_FIND_ROOT_PATH) because users may have libraries anywhere on
# the shared filesystems, and this may lie outside the root.  Instead, we set the
# system directories so that the various system BG/Q CNK library locations are
# searched first.  This is not the clearest thing in the world, given IBM's driver
# layout, but this should cover all the standard ones.
#
set(CMAKE_SYSTEM_LIBRARY_PATH
  /bgsys/drivers/ppcfloor/comm/xl/lib                       # default comm layer (used by mpi compiler wrappers)
  /bgsys/drivers/ppcfloor/spi/lib/                          # other low-level stuff
  /bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/lib # CNK Linux image -- standard runtime libs, pthread, etc.
)

#
# This adds directories that find commands should specifically ignore for cross compiles.
# Most of these directories are the includeand lib directories for the frontend on BG/Q systems.
# Not ignoring these can cause things like FindX11 to find a frontend PPC version mistakenly.
# We use this on BG instead of re-rooting because backend libraries are typically strewn about
# the filesystem, and we can't re-root ALL backend libraries to a single place.
#
set(CMAKE_SYSTEM_IGNORE_PATH
  /lib             /lib64             /include
  /usr/lib         /usr/lib64         /usr/include
  /usr/local/lib   /usr/local/lib64   /usr/local/include
  /usr/X11/lib     /usr/X11/lib64     /usr/X11/include
  /usr/lib/X11     /usr/lib64/X11     /usr/include/X11
  /usr/X11R6/lib   /usr/X11R6/lib64   /usr/X11R6/include
  /usr/X11R7/lib   /usr/X11R7/lib64   /usr/X11R7/include
)

#
# Indicate that this is a unix-like system
#
set(UNIX 1)

#
# Library prefixes, suffixes, extra libs.
#
set(CMAKE_LINK_LIBRARY_SUFFIX "")
set(CMAKE_STATIC_LIBRARY_PREFIX "lib")     # lib
set(CMAKE_STATIC_LIBRARY_SUFFIX ".a")      # .a

set(CMAKE_SHARED_LIBRARY_PREFIX "lib")     # lib
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")     # .so
set(CMAKE_EXECUTABLE_SUFFIX "")            # .exe
set(CMAKE_DL_LIBS "dl")

#
# This macro needs to be called for dynamic library support.  Unfortunately on BG/Q,
# We can't support both static and dynamic links in the same platform file.  The
# dynamic link platform file needs to call this explicitly to set up dynamic linking.
#
macro(__BlueGeneQ_set_dynamic_flags compiler_id lang)
  if (${compiler_id} STREQUAL XL)
    # Flags for XL compilers if we explicitly detected XL
    set(CMAKE_SHARED_LIBRARY_${lang}_FLAGS           "-qpic")
    set(CMAKE_SHARED_LIBRARY_CREATE_${lang}_FLAGS    "-qmkshrobj -qnostaticlink")
    set(BG/Q_${lang}_DYNAMIC_EXE_FLAGS                "-qnostaticlink -qnostaticlink=libgcc")
  else()
    # Assume flags for GNU compilers (if the ID is GNU *or* anything else).
    set(CMAKE_SHARED_LIBRARY_${lang}_FLAGS           "-fPIC")
    set(CMAKE_SHARED_LIBRARY_CREATE_${lang}_FLAGS    "-shared")
    set(BG/Q_${lang}_DYNAMIC_EXE_FLAGS                "-dynamic")
  endif()

  # Both toolchains use the GNU linker on BG/Q, so these options are shared.
  set(CMAKE_SHARED_LIBRARY_RUNTIME_${lang}_FLAG      "-Wl,-rpath,")
  set(CMAKE_SHARED_LIBRARY_RPATH_LINK_${lang}_FLAG   "-Wl,-rpath-link,")
  set(CMAKE_SHARED_LIBRARY_SONAME_${lang}_FLAG       "-Wl,-soname,")
  set(CMAKE_EXE_EXPORTS_${lang}_FLAG                 "-Wl,--export-dynamic")
  set(CMAKE_SHARED_LIBRARY_LINK_${lang}_FLAGS        "")  # +s, flag for exe link to use shared lib
  set(CMAKE_SHARED_LIBRARY_RUNTIME_${lang}_FLAG_SEP  ":") # : or empty

  set(BG/Q_${lang}_DEFAULT_EXE_FLAGS
    "<FLAGS> <CMAKE_${lang}_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> <LINK_LIBRARIES>")
  set(CMAKE_${lang}_LINK_EXECUTABLE
    "<CMAKE_${lang}_COMPILER> ${BG/Q_${lang}_DYNAMIC_EXE_FLAGS} ${BG/Q_${lang}_DEFAULT_EXE_FLAGS}")
endmacro()

#
# This macro needs to be called for static builds.  Right now it just adds -Wl,-relax
# to the link line.
#
macro(__BlueGeneQ_set_static_flags compiler_id lang)
  set(BG/Q_${lang}_DEFAULT_EXE_FLAGS
    "<FLAGS> <CMAKE_${lang}_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> <LINK_LIBRARIES>")
  set(CMAKE_${lang}_LINK_EXECUTABLE
    "<CMAKE_${lang}_COMPILER> ${BG/Q_${lang}_DEFAULT_EXE_FLAGS}")

  if(CMAKE_BUILD_TYPE MATCHES "Deb" AND ${compiler_id} STREQUAL "XL")
      # Work around an unknown compiler bug triggered in
      # compute_globals(). Using -O0 disables -qhot and this seems
      # to break the normal OpenMP flag -qsmp unless qualified with
      # noauto.
      set(OpenMP_C_FLAGS "-qsmp=noauto" CACHE STRING "Compiler flag for OpenMP parallelization")
      set(OpenMP_CXX_FLAGS "-qsmp=noauto" CACHE STRING "Compiler flag for OpenMP parallelization")
  endif()

endmacro()
