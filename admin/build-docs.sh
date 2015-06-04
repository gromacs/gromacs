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
# Documentation_* jobs.  The main purpose of this script is to isolate Jenkins
# from the details of building the documentation, as the build system may still
# evolve.  All code here is tightly tied to the build system, and the interface
# to Jenkins tries to be as static as possible to avoid the need to rebase
# changes because of job configuration changes.
#
# Interface from Jenkins to this script is:
#  * the location of this script file,
#  * the first argument is the build type ("gerrit" or "nightly"),
#  * any additional arguments should be of form -Dxxx=yyy and are passed
#    to CMake to specify details related to Jenkins configuration, such as the
#    locations of the needed executables (e.g., DOXYGEN_EXECUTABLE).
#
# Communication back to Jenkins is through
#  * return value of this script (non-zero marks the build failed).
#  * stdout (any instance of "FAILED" marks the build unstable),
#  * HTML documentation produced at build/docs/html/, and
#  * log files under build/docs/logs/:
#     - fail-reason.log contains a short description of what failed/was unstable
#     - all .log files from this directory are published as Jenkins artifacts
#     - all .log files under doxygen/ subdir are scanned for Doxygen warnings
#     - all .log files under sphinx/ subdir are scanned for Sphinx warnings

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

FAILED=
function report_failure() {
    echo "$1"
    echo "$1" >> docs/logs/fail-reason.log
    FAILED=1
}
function report_unstable() {
    echo "FAILED: $1"
    echo "$1" >> docs/logs/fail-reason.log
}
function exit_if_failed() {
    if [[ $FAILED ]] ; then
        echo "Documentation build FAILED:"
        cat docs/logs/fail-reason.log
        exit 1
    fi
}

set -x

srcdir=`git rev-parse --show-toplevel`
cd $srcdir
if [ ! -f "admin/build-docs.sh" ] ; then
    echo "Failed to find root of the source tree"
    exit 1
fi

# Some of the documentation targets can only be built out of source.
rm -rf build
mkdir build
cd build

# Copy the output logs to a static hierarchy to allow changing the exact file
# names and adding new logs without changing the Jenkins job configuration.
mkdir docs
mkdir docs/logs

cmake_args="-DGMX_BUILD_HELP=ON -DGMX_BUILD_MANUAL=ON"
if [[ $GERRIT_MODE ]] ; then
    cmake_args+=" -DGMX_COMPACT_DOXYGEN=ON"
fi
# Make a configuration that is fast, just for building docs.
cmake_args+=" -DCMAKE_BUILD_TYPE=Debug -DGMX_OPENMP=OFF -DGMX_SIMD=None -DGMX_GPU=OFF"
cmake ..  ${cmake_args} "$@" || report_failure "CMake configuration failed"
exit_if_failed

# Need to make gmx to export rst help.
make gmx -j 2 || report_failure "Building gmx failed"

# Do various parts individually to report the errors.
make manual || report_failure "PDF manual failed to build"
cp docs/manual/gromacs.log docs/logs/
grep "LaTeX Warning: Reference .* on page .* undefined" docs/logs/gromacs.log && report_unstable "undefined references in manual"

make doxygen-all || report_failure "Doxygen documentation failed to build"

# run check-source
if [[ $GERRIT_MODE ]] ; then
    make check-source || report_failure "check-source failed to run"
fi
mkdir docs/logs/doxygen
cp docs/doxygen/*.log docs/logs/doxygen/
# TODO: Check also the other logs.
[[ -s docs/logs/check-source.log ]] && report_unstable "check-source found issues"
exit_if_failed

# Generate Sphinx input.
make sphinx-input || report_failure "Generating Sphinx input failed"
make sphinx-programs || report_failure "Running gmx help -export rst failed"
exit_if_failed

# Run various Sphinx commands
make webpage-sphinx || report_failure "HTML generation with Sphinx failed"
make man || report_failure "Man page generation with Sphinx failed"
make install-guide || report_failure "INSTALL generation with Sphinx failed"
mkdir docs/logs/sphinx
cp docs/sphinx-*.log docs/logs/sphinx/
# TODO: Report unstability reason if any of these logs are non-empty.
exit_if_failed

# webpage target makes some final work for the documentation.
make webpage || report_failure "webpage target failed to build"

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

exit_if_failed
# Exit with zero code to avoid failing builds because of whatever was the last
# command executed.
exit 0
