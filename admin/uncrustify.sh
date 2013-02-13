#!/bin/bash

# Parse command-line arguments
action="check-workdir"
baserev="HEAD"
force=
for arg in "$@" ; do
    if [[ "$arg" == "check-index" || "$arg" == "check-workdir" || \
          "$arg" == "update-index" || "$arg" == "update-workdir" ]]
    then
        action=$arg
    elif [[ "$arg" == "-f" || "$arg" == "--force" ]] ; then
        force=1
    elif [[ "$arg" == --rev=* ]] ; then
        baserev=${arg#--rev=}
    else
        echo "Unknown option: $arg"
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
cfg_file=`pwd`/admin/uncrustify.cfg
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
if [[ $action == "check-index" || $action == "update-index" ]]
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

# Extract changed files to a temporary directory
mkdir $tmpdir/changed
IFS='
'
if [[ $action == *-index ]] ; then
    for change in `cat $tmpdir/filtered` ; do
        # NOTE: the patterns below contain literal tabs
        path=${change#*	}
        path=${path%%	*}
        sha1=${change:56:40}
        mkdir -p $tmpdir/changed/`dirname $path`
        git cat-file -p $sha1 > $tmpdir/changed/$path
    done
else
    for path in `cat $tmpdir/filelist` ; do
        mkdir -p $tmpdir/changed/`dirname $path`
        cp $path $tmpdir/changed/$path
    done
fi

# Run uncrustify on the temporary directory
pushd $tmpdir/changed >/dev/null 2>&1
if ! $UNCRUSTIFY -c $cfg_file -F $tmpdir/filelist -q ; then
    echo "Reformatting failed!"
    rm -rf $tmpdir
    exit 2
fi
popd >/dev/null 2>&1

# Check the differences or apply them
error=
changes=0
case $action in
    check-index|check-workdir)
        for path in `cat $tmpdir/filelist` ; do
            if ! diff $tmpdir/changed/$path $tmpdir/changed/$path.uncrustify >/dev/null ; then
                echo "$path: requires reformatting"
                changes=1
            fi
        done
        ;;
    update-index)
        for change in `cat $tmpdir/filtered` ; do
            # NOTE: the patterns below contain literal tabs
            path=${change#*	}
            path=${path%%	*}
            if ! diff $tmpdir/changed/$path.uncrustify $path >/dev/null ; then
                if [ ! $force ] && ! git diff-files --quiet $path ; then
                    echo "$path: has modifications in work tree"
                    error=1
                    continue
                fi
                mode=${change:8:6}
                sha1=${change:56:40}
                new_sha1=`git hash-object -w --stdin --path="$path" <$tmpdir/changed/$path.uncrustify`
                git update-index --cacheinfo $mode $new_sha1 $path
                echo "$path: reformatted"
                changes=1
            fi
        done
        ;;
    update-workdir)
        for path in `cat $tmpdir/filelist` ; do
            if ! diff $tmpdir/changed/$path.uncrustify $path >/dev/null ; then
                if [ ! $force ] && ! git diff-files --quiet $path ; then
                    echo "$path: has modifications in work tree"
                    error=1
                    continue
                fi
                cp $tmpdir/changed/$path.uncrustify $path
                echo "$path: reformatted"
                changes=1
            fi
        done
        ;;
esac

unset IFS
rm -rf $tmpdir
if [ $error ] ; then exit 2; fi
exit $changes
