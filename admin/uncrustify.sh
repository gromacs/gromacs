#!/bin/bash
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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

# This script runs uncrustify on modified files and reports/applies the
# results.  By default, the current HEAD commit is compared to the work tree,
# and files that
#   1. are different between these two trees and
#   2. change under uncrustify
# are reported.  This behavior can be changed by
#  1. Specifying an --rev=REV argument, which uses REV instead of HEAD as
#     the base of the comparison.
#  2. Specifying an action:
#       check-* reports the files that uncrustify changes
#       diff-* prints the actual diff of what would change
#       update-* applies the changes to the repository
#       *-workdir operates on the working directory (files on disk)
#       *-index operates on the index of the repository
# By default, update-* refuses to update "dirty" files (i.e., that differ
# between the disk and the index) to make it easy to revert the changes.
# This can be overridden by adding a -f/--force option.
#
# The location of the uncrustify executable must be specified using an
# environment variable UNCRUSTIFY.  Note that to produce the indentation used
# in the source tree, you need a custom version of uncrustify.
#
# To identify which files to run through uncrustify, the script uses git
# filters, specified in .gitattributes files.  Only files that have the filter
# set as "uncrustify" (or something ending in "uncrustify") are processed: if
# other files have been changed, they are ignored by the script.
#
# If you want to run uncrustify automatically for changes you make, there are
# two options:
#  1. Copy the git-pre-commit script in this directory to .git/hooks/pre-commit
#     and set
#       git config hooks.uncrustifymode check
#       git config hooks.uncrustifypath /path/to/uncrustify
#     See comments in the hook script for more information.
#  2. Configure a git filter (doesn't require this script, only the
#     .gitattributes files):
#       git config filter.uncrustify.clean \
#           "/path/to/uncrustify -c admin/uncrustify.cfg -q -l cpp"
# The pre-commit hook + manually running the script gives better/more intuitive
# control (with the filter, it is possible to have a work tree that is
# different from HEAD and still have an empty 'git diff') and provides better
# performance for changes that modify many files.
# The filter allows one to transparently merge branches that have not been run
# through uncrustify, and is applied more consistently (the pre-commit hook is
# not run for every commit, e.g., during a rebase).

# Parse command-line arguments
function usage() {
    echo "usage: uncrustify.sh [-f|--force] [--rev=REV] [action]"
    echo "actions: (check|diff|update)-(index|workdir) (default:check-workdir)"
}

action="check-workdir"
baserev="HEAD"
force=
for arg in "$@" ; do
    if [[ "$arg" == "check-index" || "$arg" == "check-workdir" || \
          "$arg" == "diff-index" || "$arg" == "diff-workdir" || \
          "$arg" == "update-index" || "$arg" == "update-workdir" ]]
    then
        action=$arg
    elif [[ "$arg" == "-f" || "$arg" == "--force" ]] ; then
        force=1
    elif [[ "$arg" == --rev=* ]] ; then
        baserev=${arg#--rev=}
    elif [[ "$arg" == "-h" ]] ; then
        usage
        exit 0
    else
        echo "Unknown option: $arg"
        echo
        usage
        exit 2
    fi
done

# Check that uncrustify is present
if [ -z "$UNCRUSTIFY" ]
then
    echo "Please set the path to uncrustify using UNCRUSTIFY."
    echo "Note that you need a custom version of uncrustify."
    exit 2
fi
if ! which -s "$UNCRUSTIFY"
then
    echo "Uncrustify not found: $UNCRUSTIFY"
    exit 2
fi

# Switch to the root of the source tree and check the config file
admin_dir=`dirname "$0"`
cd $admin_dir/..
srcdir=`pwd`
admin_dir=$srcdir/admin
cfg_file=$admin_dir/uncrustify.cfg
if [ ! -f "$cfg_file" ]
then
    echo "Uncrustify configuration file not found: $cfg_file"
    exit 2
fi

# Actual processing starts: create a temporary directory
tmpdir=`mktemp -d -t gmxuncrust`

# Produce a list of changed files
# Only include files that have uncrustify set as filter in .gitattributes
diff_args=
if [[ $action == *-index ]]
then
    diff_args="--cached"
fi
git diff-index $diff_args --diff-filter=ACMR $baserev >$tmpdir/difflist
cut -f2 <$tmpdir/difflist | \
    git check-attr --stdin filter | \
    sed -e 's/.*: filter: //' | \
    paste $tmpdir/difflist - | \
    grep 'uncrustify$' >$tmpdir/filtered
cut -f2 <$tmpdir/filtered >$tmpdir/filelist
git diff-files --name-only | grep -Ff $tmpdir/filelist >$tmpdir/localmods

# Extract changed files to a temporary directory
mkdir $tmpdir/org
mkdir $tmpdir/uncrustify
if [[ $action == *-index ]] ; then
    git checkout-index --prefix=$tmpdir/org/ --stdin <$tmpdir/filelist
else
    rsync --files-from=$tmpdir/filelist $srcdir $tmpdir/org
fi

# Run uncrustify on the temporary directory
cd $tmpdir/org

if ! $UNCRUSTIFY -c $cfg_file -F $tmpdir/filelist --prefix=../uncrustify/ -q ; then
    echo "Reformatting failed!"
    rm -rf $tmpdir
    exit 2
fi
# TODO: Consider checking whether rerunning uncrustify causes additional changes

cd $tmpdir

# If a diff was requested, show it and we are done
if [[ $action == diff-* ]] ; then
    # TODO: It would be nicer to have prefixes of equal length, or to get rid
    # of them altogether.
    # TODO: Make it possible to pass options to alter the diff output.
    git diff --no-index --no-prefix org/ uncrustify/
    rm -rf $tmpdir
    exit 0
fi

# Find the changed files
touch $tmpdir/messages
changes=
set -o pipefail
if ! git diff --no-index --name-only --exit-code org/ uncrustify/ | \
         sed -e 's#uncrustify/##' > $tmpdir/changed
then
    changes=1
    awk '{print $0 ": needs uncrustify"}' $tmpdir/changed \
        >> $tmpdir/messages
fi
# Check if changed files have changed outside the index
if grep -Ff $tmpdir/localmods $tmpdir/changed > $tmpdir/conflicts
then
    awk '{print $0 ": has changes in work tree"}' $tmpdir/conflicts \
        >> $tmpdir/messages
    if [[ ! $force && $action == update-* ]] ; then
        sort $tmpdir/messages
        echo "Modified files found in work tree, skipping update."
        rm -rf $tmpdir
        exit 2
    fi
fi

# Update the index/work tree if requested
if [[ $action == update-index ]] ; then
    grep -Ff $tmpdir/changed $tmpdir/filtered > $tmpdir/tohash
    cd $tmpdir/uncrustify
    IFS='
'
    for change in `cut -f2 $tmpdir/tohash | \
                   git --git-dir=$srcdir/.git hash-object -w --stdin-paths --no-filters | \
                   paste - $tmpdir/tohash`
    do
        # NOTE: the patterns below contain literal tabs
        sha1=${change%%	*}
        rest=${change#*	}
        mode=${rest:8:6}
        path=${rest#*	}
        path=${path%%	*}
        # Contains a literal tab
        echo "$mode $sha1	$path" >> $tmpdir/toindex
    done
    unset IFS
    git --git-dir=$srcdir/.git update-index --index-info < $tmpdir/toindex
elif [[ $action == update-workdir ]] ; then
    rsync --files-from=$tmpdir/changed $tmpdir/uncrustify/ $srcdir/
fi

# Report what was done
if [[ $action == update-* ]] ; then
    sort $tmpdir/messages | sed -e 's/needs uncrustify/uncrustified/'
else
    sort $tmpdir/messages
fi

rm -rf $tmpdir
exit $changes
