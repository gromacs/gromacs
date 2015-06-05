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

# This script is used by Gerrit to run uncrustify in the uncrustify_* jobs.
# The main purpose of this script is to isolate Jenkins from the details of
# the uncrustify.sh script, and to provide a facade that can keep the interface
# to Jenkins the same even if uncrustify.sh changes.
# It also provides some logic to produce better messages for Jenkins.
#
# Interface from Jenkins to this script is:
#  * location of the script
#  * UNCRUSTIFY environment variable to specify the location of uncrustify
#
# Communication back to Jenkins is through
#  * return value of this script (non-zero marks the build failed).
#  * stdout (any instance of "FAILED" marks the build unstable),
#  * unsuccessful-reason.log contains a short description of what failed
#    or was unstable
#  * uncrustify.log contains messages about what files needed changes

srcdir=`git rev-parse --show-toplevel`
cd $srcdir
if [ ! -f "admin/run-uncrustify.sh" ] ; then
    echo "Failed to find root of the source tree"
    exit 1
fi

warning_log=uncrustify.log
unsuccessful_log=unsuccessful-reason.log

admin/uncrustify.sh check --rev=HEAD^ --warnings=$warning_log
stat=$?
if [ $stat -eq 1 ] ; then
    echo "FAILED: uncrustify.sh found issues"
    warning_count=`wc -l <$warning_log`
    if [ $warning_count -lt 5 ] ; then
        cp -f $warning_log $unsuccessful_log
    else
        rm -f $unsuccessful_log
        uncrustify_count=`grep "needs uncrustify" $warning_log | wc -l`
        cpyear_count=`grep "copyright year" $warning_log | wc -l`
        cpheader_count=`grep "copyright header" $warning_log | wc -l`
        [ $uncrustify_count -gt 0 ] && echo "formatting issues in" $uncrustify_count "file(s)" >> $unsuccessful_log
        [ $cpyear_count -gt 0 ] && echo "copyright year missing in" $cpyear_count "file(s)" >> $unsuccessful_log
        [ $cpheader_count -gt 0 ] && echo "copyright header issues in" $cpheader_count "file(s)" >> $unsuccessful_log
        cat $unsuccessful_log
    fi
    exit 0
elif [ $stat -ne 0 ] ; then
    echo "FAILED: uncrustify.sh failed to run"
    echo "uncrustify.sh failed to run" > $unsuccessful_log
    exit 1
fi
exit 0
