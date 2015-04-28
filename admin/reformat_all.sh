#!/bin/bash
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

# This script runs uncrustify, copyright header checks, or include sorter on
# all applicable files in the source tree.
#
# See `reformat_all.sh -h` for a brief usage, and docs/dev-manual/uncrustify.rst
# for more details (docs/dev-manual/gmxtree.rst for include sorter).

function usage() {
    echo "usage: reformat_all.sh [-f|--force]"
    echo "           [--filter=(uncrustify|copyright|includesort)]"
    echo "           [--pattern=<pattern>] [<action>] [-B=<build dir>]"
    echo "<action>: (list-files|uncrustify*|copyright|includesort) (*=default)"
}

filter=default
force=
patterns=()
action=uncrustify
for arg in "$@" ; do
    if [[ "$arg" == "list-files" || "$arg" == "uncrustify" ||
          "$arg" == "copyright" || "$arg" == "includesort" ]] ; then
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
    uncrustify)
        # Check that uncrustify is present
        if [ -z "$UNCRUSTIFY" ] ; then
            echo "Please set the path to uncrustify using UNCRUSTIFY."
            echo "Note that you need a custom version of uncrustify."
            echo "See comments in the script file for how to get one."
            exit 2
        fi
        if ! which "$UNCRUSTIFY" 1>/dev/null ; then
            echo "Uncrustify not found: $UNCRUSTIFY"
            exit 2
        fi
        command="xargs $UNCRUSTIFY -c admin/uncrustify.cfg --no-backup"
        ;;
    copyright)
        command="xargs admin/copyright.py --check"
        ;;
    includesort)
        if [ -z "${builddir}" ] ; then
            echo "Build directory must be set with -B for includesort."
            exit 2
        fi
        command="docs/doxygen/includesorter.py -S . -B ${builddir} -F -"
        ;;
    *)
        echo "Unknown action: $action"
        exit 2
esac

if [[ "$filter" == "default" ]] ; then
    filter=$action
fi

case "$filter" in
    includesort)
        filter_re="(uncrustify|includesort)"
        ;;
    uncrustify)
        filter_re="(uncrustify|uncrustify_only)"
        ;;
    copyright)
        filter_re="(uncrustify|copyright|includesort)"
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
