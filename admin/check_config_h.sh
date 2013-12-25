#!/bin/bash

awk 'BEGIN {FS=" +|\\("} /#.*define / {print $2}' src/config.h.cmakein | \
    egrep -v '(inline|gmx_restrict|VERSION|FULLINDIRECT|USE_FAH_XDR)' >defines.txt
xargs grep -wFf defines.txt \
    < ../build/doxygen/depgraphs/installed-headers.txt
for def in `cat defines.txt` ; do
    count=`git grep -wF $def | wc -l`
    if [[ $count -eq 1 ]] ; then
        echo "Unused: $def"
    fi
done
