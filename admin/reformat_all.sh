#!/bin/bash
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
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

# This script runs clang-format, copyright header checks, or include sorter on
# all applicable files in the source tree.
#
# See `reformat_all.sh -h` for a brief usage, and docs/dev-manual/code-formatting.rst
# for more details (docs/dev-manual/gmxtree.rst for include sorter).

function usage() {
    echo "usage: reformat_all.sh [-f|--force]"
    echo "           [--filter=(complete_formatting|clangformat|copyright)]"
    echo "           [--pattern=<pattern>] [<action>] [-B=<build dir>]"
    echo "<action>: (list-files|clang-format*|copyright) (*=default)"
}

filter=default
force=
patterns=()
action=clang-format
for arg in "$@" ; do
    if [[ "$arg" == "list-files" ||
          "$arg" == "clang-format" || "$arg" == "copyright" ]] ; then
        action=$arg
    elif [[ "$arg" == --filter=* ]] ; then
        filter=${arg#--filter=}
    elif [[ "$arg" == --pattern=* ]] ; then
        patterns[${#patterns[@]}]=${arg#--pattern=}
    elif [[ "$arg" == -B=* ]] ; then
        builddir=${arg#-B=}
    elif [[ "$arg" == "-f" || "$arg" == "--force" ]] ; then
        force=1
    elif [[ "$arg" == "-h" || "$arg" == "--help" ]] ; then
        usage
        exit 0
    else
        echo "Unknown option: $arg"
        echo
        usage
        exit 2
    fi
done

if [[ ! "$force" && "$action" != "list-files" ]] ; then
    if ! git diff-files --quiet ; then
        echo "Modified files found in work tree. Use -f to override."
        exit 1
    fi
fi

case "$action" in
    list-files)
        command=cat
        ;;
    clang-format)
        if [ -z "$CLANG_FORMAT" ] ; then
            CLANG_FORMAT=clang-format-7
        fi
        if ! which "$CLANG_FORMAT" 1>/dev/null ; then
            echo "clang-format not found. Specify one with CLANG_FORMAT"
            exit 2
        fi
        command="xargs $CLANG_FORMAT -i"
        ;;
    copyright)
        command="xargs admin/copyright.py --check"
        ;;
    *)
        echo "Unknown action: $action"
        exit 2
esac

if [[ "$filter" == "default" ]] ; then
    if [[ "$action" == "clang-format" ]] ; then
        filter=clangformat
    else
        filter=$action
    fi
fi

case "$filter" in
    copyright)
        filter_re="(complete_formatting|copyright)"
        ;;
    clangformat)
        filter_re="(complete_formatting|clangformat)"
        ;;
    complete_formatting)
        filter_re="(complete_formatting|clangformat)"
        ;;
    *)
        echo "Unknown filter mode: $filter"
        echo
        usage
        exit 2
esac

cd `git rev-parse --show-toplevel`

if ! git ls-files "${patterns[@]}" | git check-attr --stdin filter | \
    sed -nEe "/${filter_re}$/ {s/:.*//;p;}" | $command ; then
    echo "The reformatting command failed! Please check the output."
    exit 1
fi
