#!/bin/bash
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

# This script runs uncrustify on modified files and reports/applies the
# results.  It also checks for outdated copyright years or headers.
# By default, the current HEAD commit is compared to the work tree,
# and files that
#  1. are different between these two trees and
#  2. change under uncrustify and/or have outdated copyright header
# are reported.  This behavior can be changed by
#  1. Specifying an --rev=REV argument, which uses REV instead of HEAD as
#     the base of the comparison.
#  2. Specifying an action:
#       check-*:   reports the files that uncrustify changes
#       diff-*:    prints the actual diff of what would change
#       update-*:  applies the changes to the repository
#       *-workdir: operates on the working directory (files on disk)
#       *-index:   operates on the index of the repository
#     For convenience, if you omit the workdir/index suffix, workdir is assumed
#     (i.e., diff equals diff-workdir).
#  3. Specifying --uncrustify=off, which does not run uncrustify.
#  4. Specifying --copyright=<mode>, which alters the level of copyright
#     checking is done:
#       off:     does not check copyright headers at all
#       year:    only update copyright year in new-format copyright headers
#       add:     in addition to 'year', add copyright headers to files that
#                don't have any
#       update:  in addition to 'year' and 'add', also update new-format
#                copyright headers if they are broken or outdated.
#       replace: replace any copyright header with a new-format copyright
#                header
#       full:    do all of the above
# By default, update-* refuses to update "dirty" files (i.e., that differ
# between the disk and the index) to make it easy to revert the changes.
# This can be overridden by adding a -f/--force option.
#
# The location of the uncrustify executable must be specified using an
# environment variable UNCRUSTIFY.  Note that to produce the indentation used
# in the source tree, you need a custom version of uncrustify.
# To get and configure the currently used version, you can do the following:
#  1. Run
#       git clone -b gromacs git://github.com/rolandschulz/uncrustify.git
#       cd uncrustify
#       ./configure
#       make
#  2. Copy the binary src/uncrustify into a directory of your choice.
#  3. Set the UNCRUSTIFY environment variable to point to the copied binary.
#
# To identify which files to run through uncrustify, the script uses git
# filters, specified in .gitattributes files.  Only files that have the filter
# set as "uncrustify" (or something ending in "uncrustify") are processed: if
# other files have been changed, they are ignored by the script.  Files passed
# to uncrustify, as well as files with filter "copyright", get their copyright
# header checked.  To only run uncrustify for a file, set the filter to
# "uncrustify_only".
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
# performance for changes that modify many files.  It is the only way that
# currently also checks the copyright headers.
# The filter allows one to transparently merge branches that have not been run
# through uncrustify, and is applied more consistently (the pre-commit hook is
# not run for every commit, e.g., during a rebase).

# Parse command-line arguments
function usage() {
    echo "usage: uncrustify.sh [-f|--force] [--rev=REV]"
    echo "           [--uncrustify=(off|check)] [--copyright=<cmode>] [<action>]"
    echo "<action>: (check*|diff|update)[-(index|workdir*)] (*=default)"
    echo "<cmode>:  off|add|update*|replace|full"
}

action="check-workdir"
declare -a diffargs
baserev="HEAD"
force=
uncrustify_mode=check
copyright_mode=update
for arg in "$@" ; do
    if [[ "$arg" == "check-index" || "$arg" == "check-workdir" || \
          "$arg" == "diff-index" || "$arg" == "diff-workdir" || \
          "$arg" == "update-index" || "$arg" == "update-workdir" ]]
    then
        action=$arg
    elif [[ "$arg" == "check" || "$arg" == "diff" || "$arg" == "update" ]] ; then
        action=$arg-workdir
    elif [[ "$action" == diff-* ]] ; then
        diffargs+=("$arg")
    elif [[ "$arg" == --uncrustify=* ]] ; then
        uncrustify_mode=${arg#--uncrustify=}
        if [[ "$uncrustify_mode" != "off" && "$uncrustify_mode" != "check" ]] ; then
            echo "Unknown option: $arg"
            echo
            usage
            exit 2
        fi
    elif [[ "$arg" == --copyright=* ]] ; then
        copyright_mode=${arg#--copyright=}
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
if [[ "$uncrustify_mode" != "off" ]]
then
    if [ -z "$UNCRUSTIFY" ]
    then
        echo "Please set the path to uncrustify using UNCRUSTIFY."
        echo "Note that you need a custom version of uncrustify."
        echo "See comments in the script file for how to get one."
        exit 2
    fi
    if ! which "$UNCRUSTIFY" 1>/dev/null
    then
        echo "Uncrustify not found: $UNCRUSTIFY"
        exit 2
    fi
fi

# Switch to the root of the source tree and check the config file
srcdir=`git rev-parse --show-toplevel`
cd $srcdir
admin_dir=$srcdir/admin
cfg_file=$admin_dir/uncrustify.cfg
if [ ! -f "$cfg_file" ]
then
    echo "Uncrustify configuration file not found: $cfg_file"
    exit 2
fi

# Actual processing starts: create a temporary directory
tmpdir=`mktemp -d -t gmxuncrust.XXXXXX`

# Produce a list of changed files
# Only include files that have proper filter set in .gitattributes
internal_diff_args=
if [[ $action == *-index ]]
then
    internal_diff_args="--cached"
fi
git diff-index $internal_diff_args --diff-filter=ACMR $baserev >$tmpdir/difflist
cut -f2 <$tmpdir/difflist | \
    git check-attr --stdin filter | \
    sed -e 's/.*: filter: //' | \
    paste $tmpdir/difflist - | \
    grep -E '(uncrustify|uncrustify_only|copyright)$' >$tmpdir/filtered
cut -f2 <$tmpdir/filtered >$tmpdir/filelist_all
grep -E '(uncrustify|uncrustify_only)$' <$tmpdir/filtered | \
    cut -f2 >$tmpdir/filelist_uncrustify
grep -E '(uncrustify|copyright)$' <$tmpdir/filtered | \
    cut -f2 >$tmpdir/filelist_copyright
git diff-files --name-only | grep -Ff $tmpdir/filelist_all >$tmpdir/localmods

# Extract changed files to a temporary directory
mkdir $tmpdir/org
if [[ $action == *-index ]] ; then
    git checkout-index --prefix=$tmpdir/org/ --stdin <$tmpdir/filelist_all
else
    rsync --files-from=$tmpdir/filelist_all $srcdir $tmpdir/org
fi
# Duplicate the original files to a separate directory, where all changes will
# be made.
cp -r $tmpdir/org $tmpdir/new

# Create output file for what was done (in case no messages get written)
touch $tmpdir/messages

# Run uncrustify on the temporary directory
cd $tmpdir/new
if [[ $uncrustify_mode != "off" ]] ; then
    if ! $UNCRUSTIFY -c $cfg_file -F $tmpdir/filelist_uncrustify --no-backup >$tmpdir/uncrustify.out 2>&1 ; then
        echo "Reformatting failed. Check uncrustify output below for errors:"
        cat $tmpdir/uncrustify.out
        rm -rf $tmpdir
        exit 2
    fi
    # Find the changed files if necessary
    if [[ $action != diff-* ]] ; then
        msg="needs uncrustify"
        if [[ $action == update-* ]] ; then
            msg="uncrustified"
        fi
        git diff --no-index --name-only ../org/ . | \
            awk -v msg="$msg" '{sub(/.\//,""); print $0 ": " msg}' >> $tmpdir/messages
    fi
    # TODO: Consider checking whether rerunning uncrustify causes additional changes
fi

# Update the copyright headers using the requested mode
if [[ $copyright_mode != "off" ]] ; then
    cpscript_args="--update-year"
    case "$copyright_mode" in
        year)
            ;;
        add)
            cpscript_args+=" --add-missing"
            ;;
        update)
            cpscript_args+=" --add-missing --update-header"
            ;;
        replace)
            cpscript_args+=" --replace-header"
            ;;
        full)
            cpscript_args+=" --add-missing --update-header --replace-header"
            ;;
        *)
            echo "Unknown copyright mode: $copyright_mode"
            exit 2
    esac
    if [[ $action == check-* ]] ; then
        cpscript_args+=" --check"
    fi
    # TODO: Probably better to invoke python explicitly through a customizable
    # variable.
    if ! $admin_dir/copyright.py -F $tmpdir/filelist_copyright $cpscript_args >>$tmpdir/messages
    then
        echo "Copyright checking failed!"
        rm -rf $tmpdir
        exit 2
    fi
fi

cd $tmpdir

# If a diff was requested, show it and we are done
if [[ $action == diff-* ]] ; then
    git diff --no-index --no-prefix "${diffargs[@]}" org/ new/
    rm -rf $tmpdir
    exit 0
fi

# Find the changed files
changes=
set -o pipefail
if ! git diff --no-index --name-only --exit-code org/ new/ | \
         sed -e 's#new/##' > $tmpdir/changed
then
    changes=1
fi

# Check if changed files have changed outside the index
if grep -Ff $tmpdir/localmods $tmpdir/changed > $tmpdir/conflicts
then
    awk '{print $0 ": has changes in work tree"}' $tmpdir/conflicts \
        >> $tmpdir/messages
    if [[ ! $force && $action == update-* ]] ; then
        echo "Modified files found in work tree, skipping update. Use -f to override."
        echo "The following would have been done:"
        sort $tmpdir/messages
        rm -rf $tmpdir
        exit 2
    fi
fi

# Update the index/work tree if requested
if [[ $action == update-index ]] ; then
    grep -Ff $tmpdir/changed $tmpdir/filtered > $tmpdir/tohash
    cd $tmpdir/new
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
    rsync --files-from=$tmpdir/changed $tmpdir/new/ $srcdir/
fi

# Report what was done
sort $tmpdir/messages

rm -rf $tmpdir
exit $changes
