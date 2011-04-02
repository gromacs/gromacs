#!/bin/sh

# usage: genversion.sh PKGVERSION top_srcdir
# Requires git 1.5.1 or later

PKGVERSION=$1
SRCDIR=$2
GITDIR=$SRCDIR/.git

if which git >/dev/null && test -d $GITDIR ; then
    # Get either -dirty or an empty string, depending on whether there are local changes.
    dirtystr=`(cd $SRCDIR && git diff-index --quiet HEAD) || echo "-dirty"`
    # Get date of the head commit as YYYYMMDD (commit date).
    # Git before 1.5.3 does not support any sensible date format,
    # so we need to massage the output.
    if git --git-dir=$GITDIR rev-list -n1 --pretty=format:%ci HEAD | grep '[0-9]\{4\}' >/dev/null 2>&1; then
        date=`git --git-dir=$GITDIR rev-list -n1 --pretty=format:%ci HEAD | sed -ne '/commit/!{s/ .*$//;s/-//g;p;}'`
    else
        date=`git --git-dir=$GITDIR rev-list -n1 --pretty=format:%cD HEAD | \
              sed -ne '/commit/!{s/^.*, *\([ 0-9][0-9]\) \([a-zA-Z]*\) \([0-9]*\) .*$/\3\2\1/;y/ /0/;\
                   s/Jan/01/;s/Feb/02/;s/Mar/03/;s/Apr/04/;s/May/05/;s/Jun/06/;
                   s/Jul/07/;s/Aug/08/;s/Sep/09/;s/Oct/10/;s/Nov/11/;s/Dec/12/;
                   p;}'`
    fi
    # Get a 7-character hash for the HEAD commit.
    shorthash=`git --git-dir=$GITDIR rev-parse --short=7 HEAD 2>/dev/null`
    # Get the full hash for the HEAD commit.
    fullhash=`git --git-dir=$GITDIR rev-parse HEAD 2>/dev/null`
    # Generate a version string like 4.0.99-dev-YYYYMMDD-1234abc-dirty.
    # If PKGVERSION ends in a date, replace it with the head commit.
    version=`echo $PKGVERSION | sed -e 's/-dev-[0-9]*$/-dev/'`-$date-$shorthash$dirtystr
    #version=$PKGVERSION-$date-$shorthash$dirtystr
    # Find the name of the remote which has the git.gromacs.org:gromacs in its url.
    # If not found, just set baserev to "unknown".
    gmxremote=`git --git-dir=$GITDIR config --get-regexp 'remote\..*\.url' 'git\.gromacs\.org[:|/]gromacs' | sed -e 's/remote\.\(.*\)\.url.*/\1/'`
    if test "x$gmxremote" = "x" ; then
        baserev="unknown"
    else
        # Find the most recent ancestral commit that appears in $gmxremote.
        # Gets the hash and the number of local commits that are newer than
        # the found commit.
        baserev=`git --git-dir=$GITDIR rev-list HEAD | git --git-dir=$GITDIR name-rev --stdin --refs=refs/remotes/$gmxremote/* | awk 'NF > 1 {print $1 " (" NR-1 " newer local commits)"; exit 0}'`
        # Extract the base hash
        basehash=`expr "$baserev" : '\([0123456789abcdef]*\) '`
        # Do not report the base revision if it is the same as the most recent commit
        if test "$basehash" = "$fullhash" ; then
            baserev=
        fi
    fi
else
    version=$PKGVERSION
    fullhash="unknown"
    dirtystr=
    baserev=
fi

# Write out to a temporary file, to compare with current version.c.
cat > version.c.tmp << END
#include "version.h"
const char _gmx_ver_string[] = "VERSION $version";
const char _gmx_full_git_hash[] = "$fullhash$dirtystr";
const char _gmx_central_base_hash[] = "$baserev";
END

# Replace version.c with version.c.tmp if they differ.
[ -f version.c ] || touch version.c
cmp -s version.c.tmp version.c || mv -f version.c.tmp version.c
rm -f version.c.tmp

