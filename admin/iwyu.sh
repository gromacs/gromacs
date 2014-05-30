#/bin/bash

# This script needs to be run from a build folder as a subfolder of the source.
# The first argument is the file to analysis. include-what-you-use needs to
# be in the path. Add --apply if you want changes to be applied.  Any extra 
# arguments are added as is to the command (can be used for extra defines).

filename=$1
shift

cmd="include-what-you-use -DHAVE_CONFIG_H  -I../src -I../src/external/thread_mpi/include -Isrc -Isrc/gromacs/utility \
    -Xiwyu --mapping_file=../admin/iwyu.imp -mavx"

# We cannot detect wether it is a C++ or C header. Should be find to always use C++
if [ "${filename##*.}" == "h" ]; then
    cmd="$cmd -x c++"
fi

cmd="$cmd $filename"

# Always use C++11. This is always the standard for clang.
if [ "${filename##*.}" == "cpp" -o "${filename##*.}" == "h" ]; then
    cmd="$cmd -std=c++11"
fi

# Read all special arguments and add others to the command
apply=0
for arg in "$@"; do
    if [ $arg == "--apply" ]; then
	apply=1
    else
	cmd="$cmd $arg"
    fi
done

if [ $apply -eq 1 ] ; then
    cmd="$cmd 2>&1 | fix_includes.py --nosafe_headers"
fi

`$cmd`
echo $cmd
