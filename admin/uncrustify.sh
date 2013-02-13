#!/bin/bash

action="check-workdir"
baserev="HEAD"
force=
for arg in "$@"
do
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

admin_dir=`dirname "$0"`
cfg_file="$admin_dir/uncrustify.cfg"

if [ ! -f "$cfg_file" ]
then
    echo "Uncrustify configuration file not found: $cfg_file"
    exit 2
fi

diff_args=
if [[ $action == "check-index" || $action == "update-index" ]]
then
    diff_args="--cached"
fi
IFS='
'
error=0
for change in `git diff-index $diff_args --diff-filter=ACMR $baserev`
do
    mode=`echo $change | cut -d' ' -f2`
    sha1=`echo $change | cut -d' ' -f4`
    path=`echo $change | cut -f2`
    filename=`basename $path`
    ext=${filename##*.}
    if [[ "$ext" != "c" && "$ext" != "cpp" && "$ext" != "h" ]]
    then
        continue
    fi
    if [[ "$path" == src/external/* || \
          "$path" == src/contrib/* || \
          "$path" == */gmx_blas/* || \
          "$path" == */gmx_lapack/* || \
          "$path" == */nb_kernel_Elec* || \
          "$path" == src/gromacs/selection/parser.cpp || \
          "$path" == src/gromacs/selection/parser.h || \
          "$path" == src/gromacs/selection/scanner.cpp || \
          "$path" == src/gromacs/selection/scanner_flex.h ]]
    then
        continue
    fi
    local_mods=
    if ! git diff-files --quiet $path
    then
        local_mods=1
    fi
    case $action in
        check-index)
            new_sha1=`git cat-file -p $sha1 | \
                $UNCRUSTIFY -c "$cfg_file" -l cpp 2>/dev/null | \
                git hash-object --stdin --path="$path"`
            if [ $? -ne 0 ] ; then
                echo "$path: error in checking"
                error=1
            fi
            if [ $sha1 != $new_sha1 ] ; then
                echo "$path requires reformatting"
                error=1
            fi
            ;;
        check-workdir)
            sha1=`git hash-object $path`
            if [ $? -ne 0 ] ; then
                echo "$path: error in checking"
                error=1
            fi
            new_sha1=`$UNCRUSTIFY -c "$cfg_file" -f $path 2>/dev/null | \
                git hash-object --stdin --path="$path"`
            if [ $? -ne 0 ] ; then
                echo "$path: error in checking"
                error=1
            fi
            if [ $sha1 != $new_sha1 ] ; then
                echo "$path requires reformatting"
                error=1
            fi
            ;;
        update-index)
            new_sha1=`git cat-file -p $sha1 | \
                $UNCRUSTIFY -c "$cfg_file" -l cpp 2>/dev/null | \
                git hash-object -w --stdin --path="$path"`
            if [ $? -ne 0 ] ; then
                echo "$path: error in reformatting, aborting"
                exit 2
            fi
            if [ $sha1 != $new_sha1 ] ; then
                if [[ $local_mods && ! $force ]] ; then
                    echo "$path: requires reformatting, but has modifications in work tree"
                    exit 2
                fi
                # TODO: Should check that further uncrustify runs do not cause additional changes
                git update-index --cacheinfo $mode $new_sha1 $path
                echo "$path reformatted"
                error=1
            fi
            ;;
        update-workdir)
            if [[ $local_mods && ! $force ]] ; then
                echo "$path: has modifications in work tree"
                exit 2
            fi
            sha1=`git hash-object $path`
            if [ $? -ne 0 ] ; then
                echo "$path: error in reformatting, aborting"
                exit 2
            fi
            $UNCRUSTIFY -c "$cfg_file" --no-backup $path 2>/dev/null
            # TODO: Should check that further uncrustify runs do not cause additional changes
            if [ $? -ne 0 ] ; then
                echo "$path: error in reformatting, aborting"
                exit 2
            fi
            new_sha1=`git hash-object $path`
            if [ $? -ne 0 ] ; then
                echo "$path: error in reformatting, aborting"
                exit 2
            fi
            if [ $sha1 != $new_sha1 ] ; then
                echo "$path reformatted"
                error=1
            fi
            ;;
        *)
            echo "Invalid operation mode '$action'"
            exit 2
            ;;
    esac
done
unset IFS
exit $error
