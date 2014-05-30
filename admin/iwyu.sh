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

#/bin/bash

# This script needs to be run from a build folder as a subfolder of the source.
# The first argument is the file to analysis. include-what-you-use needs to
# be in the path. Add --apply if you want changes to be applied.  Any extra 
# arguments are added as is to the command (can be used for extra defines).

filename=$1
shift

cmd="include-what-you-use -DHAVE_CONFIG_H  -I../src -I../src/external/thread_mpi/include -Isrc -Isrc/gromacs/utility \
    -Xiwyu --mapping_file=../admin/iwyu.imp -mavx"

# We cannot detect wether it is a C++ or C header. Should be find to always use C++
if [ "${filename##*.}" == "h" ]; then
    cmd="$cmd -x c++"
fi

cmd="$cmd $filename"

# Always use C++11. This is always the standard for clang.
if [ "${filename##*.}" == "cpp" -o "${filename##*.}" == "h" ]; then
    cmd="$cmd -std=c++11"
fi

# Read all special arguments and add others to the command
apply=0
for arg in "$@"; do
    if [ $arg == "--apply" ]; then
	apply=1
    else
	cmd="$cmd $arg"
    fi
done

if [ $apply -eq 1 ] ; then
    cmd="$cmd 2>&1 | fix_includes.py --nosafe_headers"
fi

`$cmd`
echo $cmd
