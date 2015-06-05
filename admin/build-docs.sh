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
# Command line options to use an existing build directory are provided for more
# easily testing the script outside Jenkins.
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

# Parse arguments
function usage() {
    echo "usage: build-docs.sh gerrit|nightly"
    echo "       [--use-build-dir=<existing build dir>] [--skip-cmake]"
    echo "       [--build-cmd=<build command>] [-D...=...]*"
    echo "All other additional options (-D etc.) are passed to CMake"
}

declare -a cmake_args

build_dir=build
build_cmd=make
do_cmake=1
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
for arg in "$@" ; do
    if [[ "$arg" == --use-build-dir=* ]] ; then
        build_dir=${arg#--use-build-dir=}
    elif [[ "$arg" == --build-cmd=* ]] ; then
        build_cmd=${arg#--build-cmd=}
    elif [[ "$arg" == --skip-cmake ]] ; then
        do_cmake=
    else
        cmake_args+=("$arg")
    fi
done

log_output_dir=docs/logs
unsuccessful_log=$log_output_dir/unsuccessful-reason.log

# Utilities for failure reporting
FAILED=
function report_failure() {
    echo "$1"
    echo "$1" >> $unsuccessful_log
    FAILED=1
}
function report_unstable() {
    echo "FAILED: $1"
    echo "$1" >> $unsuccessful_log
}
function exit_if_failed() {
    if [[ $FAILED ]] ; then
        echo "Documentation build FAILED:"
        cat $unsuccessful_log
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
mkdir -p $build_dir
cd $build_dir

# Copy the output logs to a static hierarchy to allow changing the exact file
# names and adding new logs without changing the Jenkins job configuration.
mkdir -p $log_output_dir
rm -f $unsuccessful_log

cmake_args+=("-DGMX_BUILD_HELP=ON" "-DGMX_BUILD_MANUAL=ON")
if [[ $GERRIT_MODE ]] ; then
    cmake_args+=("-DGMX_COMPACT_DOXYGEN=ON")
fi
# Make a configuration that is fast, just for building docs.
cmake_args+=("-DCMAKE_BUILD_TYPE=Debug" "-DGMX_OPENMP=OFF" "-DGMX_SIMD=None" "-DGMX_GPU=OFF")
if [[ $do_cmake ]] ; then
    cmake $srcdir  "${cmake_args[@]}" || report_failure "CMake configuration failed"
    exit_if_failed
else
    echo "Skipped cmake; args ${cmake_args[@]}"
fi

# Need to make gmx to export rst help.
$build_cmd gmx -j 2 || report_failure "Building gmx failed"

# Do various parts individually to report the errors.
$build_cmd manual || report_failure "PDF manual failed to build"
cp docs/manual/gromacs.log $log_output_dir/
grep "LaTeX Warning: Reference .* on page .* undefined" docs/manual/gromacs.log && report_unstable "undefined references in manual"

$build_cmd doxygen-all || report_failure "Doxygen documentation failed to build"

# run check-source
if [[ $GERRIT_MODE ]] ; then
    $build_cmd check-source || report_failure "check-source failed to run"
fi
mkdir $log_output_dir/doxygen
cp docs/doxygen/*.log $log_output_dir/doxygen/
for target in {check-source,doxygen-xml,doxygen-full,doxygen-lib,doxygen-user} ; do
    if [ -s $log_output_dir/doxygen/$target.log ] ; then
        report_unstable "$target produced warnings"
    fi
done
exit_if_failed

# Generate Sphinx input.
$build_cmd sphinx-input || report_failure "Generating Sphinx input failed"
$build_cmd sphinx-programs || report_failure "Running gmx help -export rst failed"
exit_if_failed

# Run various Sphinx commands
$build_cmd webpage-sphinx || report_failure "Sphinx: HTML generation failed"
$build_cmd man || report_failure "Sphinx: man page generation failed"
$build_cmd install-guide || report_failure "Sphinx: INSTALL generation failed"
mkdir ${log_output_dir}/sphinx
cp docs/sphinx-*.log ${log_output_dir}/sphinx/
for log in {html,man,install} ; do
    if [ -s $log_output_dir/sphinx/sphinx-$log ] ; then
        case $log in
            html)
                format="HTML"
                ;;
            man)
                format="man page"
                ;;
            install)
                format="INSTALL"
                ;;
        esac
        report_unstable "Sphinx: $format generation produced warnings"
    fi
done
exit_if_failed

# webpage target makes some final work for the documentation.
$build_cmd webpage || report_failure "webpage target failed to build"

cd $srcdir

if [ -f $build_dir/docs/html/index.html ] ; then
    linkchecker $build_dir/docs/html/index.html -f docs/linkcheckerrc \
        --ignore-url html-full --ignore-url html-user --ignore-url html-lib \
        --ignore-url .tar.gz --ignore-url _sources
      # add this to previous line once content stabilizes
      # || echo "FAILED linkchecker"
else
    echo "No docs/html/index.html in $build_dir was made; can't check links!"
fi

cd $build_dir

[[ -s $unsuccessful_log ]] && cat $unsuccessful_log
[[ $FAILED ]] && exit 1
exit 0
