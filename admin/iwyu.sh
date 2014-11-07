#!/bin/bash
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

# The build and source folder can be specified with -B and -S respectively.
# By default it assume "-B. -S..".
# include-what-you-use needs to be in the path. Add --apply if you want
# changes to be applied.  Any extra arguments are added as is to the
# command (can be used for extra defines or include paths).

filename=
build_path=.
src_path=..
cmd="include-what-you-use -DHAVE_CONFIG_H -mavx"

# Read all special arguments and add others to the command
apply=0
for arg in "$@"; do
    if [ $arg == "--apply" ]; then
	apply=1
    elif [[ $arg == -[SB] ]]; then
	echo -S and -B require an argument
	exit 1
    elif [[ $arg == -B* ]]; then
	build_path=${arg:2}
    elif [[ $arg == -S* ]]; then
	src_path=${arg:2}
    elif [[ $arg != -* ]]; then
	if [ "$filename" == "" ]; then
	    filename=$arg
	else
	    echo "This script can only be run on one file at a time"
	    exit 1
	fi
    else
	cmd="$cmd $arg"
    fi
done

if [ "$filename" == "" ]; then
    echo "No file specified"
    exit 1
fi

# We cannot detect wether it is a C++ or C header. Should be fine to always use C++
if [ "${filename##*.}" == "h" ]; then
    cmd="$cmd -x c++"
fi

cmd="$cmd $filename"

# Always use C++11.
if [ "${filename##*.}" == "cpp" -o "${filename##*.}" == "h" ]; then
    cmd="$cmd -std=c++11"
fi

# keep gmxpre.h for source files
if [ "${filename##*.}" == "cpp" -o "${filename##*.}" == "c" ]; then
    cmd="$cmd  -Xiwyu --pch_in_code -Xiwyu --prefix_header_includes=keep"
fi

if [ $src_path == "." ]; then
    src_folder="src" # ./src confuses IWYU
else
    src_folder="$src_path/src"
fi

cmd="$cmd -I${src_folder} -I${src_folder}/external/thread_mpi/include
     -I$build_path/src -I${src_folder}/external/boost
     -Xiwyu --mapping_file=${src_path}/admin/iwyu.imp"

if [ $apply -eq 1 ] ; then
    cmd="$cmd 2>&1 | fix_includes.py --nosafe_headers ||
         ${src_path}/docs/doxygen/includesorter.py $filename -B$build_path -S$src_path"
fi

eval $cmd
