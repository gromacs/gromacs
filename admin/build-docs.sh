#!/bin/bash
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

# This script is used by Gerrit to build the documentation in the
# Documentation_* jobs.
#
# The script takes the type of build as the first argument, and any additional
# arguments are passed to CMake to configure the build system.  These
# additional parameters can be used to, e.g., specify the locations of the
# needed executables like Doxygen.
#
# The main purpose of this script is to isolate Jenkins from the details of
# building the documentation, as the build system may still evolve.  All code
# here is tightly tied to the build system.  Communication back to Jenkins is
# through
#  * stdout (any instance of "FAILED" marks the build unstable)
#  * HTML documentation produced at build/docs/html/,
#  * some log files, and
#  * return value of this script (non-zero marks the build failed).
# This interface should be kept as static as possible to avoid the need to
# rebase changes because of job configuration changes.

function usage() {
    echo "usage: build-docs.sh gerrit|nightly [-D...=...]*"
    echo "All additional options (-D etc.) are passed to CMake"
}

GERRIT_MODE=
NIGHTLY_MODE=
case "$1" in
    gerrit)
        GERRIT_MODE=1
        ;;
    nightly) 
        NIGHTLY_MODE=1
        ;;
    *)
        usage
        exit 2
        ;;
esac
shift

set -x

srcdir=`git rev-parse --show-toplevel`
cd $srcdir
if [ ! -f "admin/build-docs.sh" ] ; then
    echo "Failed to find root of the source tree"
    exit 1
fi

# Must build out of source
rm -rf build
mkdir build
cd build

# Make a configuration that is fast, just for building docs.
cmake ..  "$@" \
    -DGMX_BUILD_HELP=yes -DGMX_BUILD_MANUAL=yes \
    -DCMAKE_BUILD_TYPE=Debug -DGMX_OPENMP=off -DGMX_SIMD=None -DGMX_GPU=no

# Need to make gmx for the html output from tools
make gmx -j 2 || echo "FAILED to make gmx"

# webpage target makes all the documentation components, but not sure this will work stably in parallel yet
make webpage || echo "FAILED to make webpage"

grep "LaTeX Warning: Reference .* on page .* undefined" docs/manual/gromacs.log && echo "FAILED - undefined references in manual"

# run check-source
if [[ $GERRIT_MODE ]] ; then
    make check-source || echo "FAILED: check-source found errors"
fi

# Ideally, we would also make the other components individually, to check
# that the targets still work, but most of them duplicate the build
cd ..

if [ -f build/docs/html/index.html ] ; then
    linkchecker build/docs/html/index.html -f docs/linkcheckerrc \
        --ignore-url html-full --ignore-url html-user --ignore-url html-lib \
        --ignore-url .tar.gz --ignore-url _sources
      # add this to previous line once content stabilizes
      # || echo "FAILED linkchecker"
else
    echo "No build/docs/html/index.html was made; can't check links!"
fi

# TODO: Consider copying the logs to a static location/hierarchy that allow
# changing the names/locations of the logs and adding new ones without
# reconfiguring the job.

# For now, don't give errors when linkchecker fails
exit 0
