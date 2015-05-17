#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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
#   GMX_FIX_OSX_RPATH(TARGET)
#
# - Define macro to fix incorrect rpath references for OS X. For this OS,
#   dynamic libaries (dylib) can have a compiled-in identifier specifiying
#   that their location is specified relative to an "rpath" location.
#   Unfortunately this is broken for some libraries, in particular Intel TBB.
#   This does not matter when libraries are installed in the default locations
#   since those are included in DYLD_LIBRARY_PATH, but for libraries in
#   non-standard directories we will get errors at execution time unless their
#   directory is specified explicitly in DYLD_LIBRARY_PATH.
#   This can be fixed by altering the id field in the library, but since that
#   requires root permission we instead alter the library specification in
#   the libraries and binaries we generate, to add the rpath prefix.
#
#   This macro takes a target as argument, it uses otool -L to list all
#   linked-in libraries, and if a library is missing a prefix (bad) we add
#   an rpath prefix. Note that you must also add the base directory to the
#   rpath when building, but without this macro even that won't help.'
#
#   The command will not have any effect for non-OSX-systems, or if the
#   paths are already correct.
#
macro(gmx_fix_osx_rpath target)
    if(APPLE)
        # For each library listed in the otool -L output, if the library has no
        # prefix we change the reference to include the @rpath/ prefix.
        add_custom_command(TARGET ${target} POST_BUILD COMMAND
                           for i in `otool -L $<TARGET_FILE:${target}> | awk '$$1~/^lib.*dylib/{print $$1}'`\;
                           do install_name_tool -change $$i @rpath/$$i $<TARGET_FILE:${target}>\; done)
    endif()
endmacro()
