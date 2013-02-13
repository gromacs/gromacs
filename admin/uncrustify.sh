#!/bin/bash

action=$1
if [[ -z "$action" ]] ; then action="check-index"; fi

if [ -z "$UNCRUSTIFY" ] ; then UNCRUSTIFY="uncrustify"; fi

if ! which -s "$UNCRUSTIFY"
then
    echo "Source file formatting check failed."
    echo "Uncrustify not found: $UNCRUSTIFY"
    exit 1
fi

admin_dir=`dirname $0`

IFS='
'
error=0
for change in `git diff-index --cached --diff-filter=ACMR HEAD`
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
    # TODO: Needs more blacklist entries and/or a nicer approach
    if [[ "$path" == src/external/* ]]
    then
        continue
    fi
    case $action in
        check-index)
            new_sha1=`git cat-file -p $sha1 | \
                $UNCRUSTIFY -c $admin_dir/uncrustify.cfg -l cpp 2>/dev/null | \
                git hash-object --stdin --path="$path"`
            if [ $sha1 != $new_sha1 ]
            then
                echo "$path requires reformatting"
                error=1
            fi
            ;;
        update-index)
            new_sha1=`git cat-file -p $sha1 | \
                $UNCRUSTIFY -c $admin_dir/uncrustify.cfg -l cpp 2>/dev/null | \
                git hash-object -w --stdin --path="$path"`
            if [ $sha1 != $new_sha1 ]
            then
                echo "$path reformatted"
                git update-index --cacheinfo $mode $new_sha1 $path
                error=1
            fi
            ;;
        *)
            echo "Invalid operation mode '$action'"
            exit 1
            ;;
    esac
done
unset IFS
exit $error
