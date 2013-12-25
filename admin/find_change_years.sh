#!/bin/bash

for file in "$@" ; do
    years=`git log --follow --pretty=format:%ci $file | sort -nu | cut -c1-4 | tr '\n' ','`
    years_no_follow=`git log --pretty=format:%ci $file | sort -nu | cut -c1-4 | tr '\n' ','`
    if [ "$years_no_follow" != "$years" ] ; then
        echo -e $file '\t' $years_no_follow '\t' $years
    else
        echo -e $file '\t' $years '\t' -
    fi
done
