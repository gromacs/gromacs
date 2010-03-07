#!/bin/sh

# usage: genversion.sh PKGVERSION top_srcdir

PKGVERSION=$1
SRCDIR=$2
GITDIR=$SRCDIR/.git

if which git >/dev/null && test -d $GITDIR ; then
    # Get either -dirty or an empty string, depending on whether there are local changes.
    dirtystr=`(cd $SRCDIR && git diff-index --quiet HEAD) || echo "-dirty"`
    # Get date of the head commit as YYYYMMDD (commit date).
    date=`git --git-dir=$GITDIR rev-list -n1 --pretty=format:%ci HEAD | sed -ne '/commit/!{s/-\| .*$//g;p}'`
    # Get a 7-character hash for the HEAD commit.
    shorthash=`git --git-dir=$GITDIR rev-parse --short=7 HEAD 2>/dev/null`
    # Get the full hash for the HEAD commit.
    fullhash=`git --git-dir=$GITDIR rev-parse HEAD 2>/dev/null`$dirtystr
    # Generate a version string like 4.0.99-dev-YYYYMMDD-1234abc-dirty.
    # If PKGVERSION ends in a date, replace it with the head commit.
    version=`echo $PKGVERSION | sed -e 's/-dev-[0-9]*$/-dev/'`-$date-$shorthash$dirtystr
    #version=$PKGVERSION-$date-$shorthash$dirtystr
    # Find the name of the remote which has the git.gromacs.org:gromacs in its url.
    # If not found, just set baserev to "unknown".
    gmxremote=`git --git-dir=$GITDIR config --get-regexp 'remote\..*\.url' 'git\.gromacs\.org:gromacs' | sed -e 's/remote\.\(.*\)\.url.*/\1/'`
    if test "x$gmxremote" = "x" ; then
        baserev="unknown"
    else
        # Find the most recent ancestral commit that appears in $gmxremote.
        # Gets the hash and the number of local commits that are newer than
        # the found commit.
        baserev=`git --git-dir=$GITDIR rev-list HEAD | git --git-dir=$GITDIR name-rev --stdin --refs=refs/remotes/$gmxremote/* | awk 'NF > 1 {print $1 " (" NR-1 " newer local commits)"; nextfile}'`
    fi
else
    version=$PKGVERSION
    fullhash="unknown"
    baserev="unknown"
fi

# write out to a temporary file, to compare with current version.c
echo "#include \"version.h\"" > version.c.tmp
echo "const char _gmx_ver_string[] = \"VERSION $version\";" >> version.c.tmp
echo "const char _gmx_full_git_hash[] = \"$fullhash\";" >> version.c.tmp
echo "const char _gmx_central_base_hash[] = \"$baserev\";" >> version.c.tmp

# Write contents into version.c.tmp if they differ from the current version.c.
[ -f version.c ] || touch version.c
cmp -s version.c.tmp version.c || mv -f version.c.tmp version.c
rm -f version.c.tmp

