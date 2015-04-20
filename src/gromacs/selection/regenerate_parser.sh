#!/bin/bash
#
# This script runs Bison and/or Flex to regenerate the files as follows:
#   parser.y  -> parser.cpp, parser.h
#   scanner.l -> scanner.cpp, scanner_flex.h
# The commands are run only if the generated files are older than the
# Bison/Flex input files, or if a '-f' flag is provided.

# Note: You can check parser.cpp/scanner.cpp for the exact versions of
# bison/flex that were used in the generation. Some OSs have older versions
# of these tools installed.

FORCE=
if [ "x$1" == "x-f" ] ; then
    FORCE=1
fi

if [ "x$BISON" == "x" ] ; then
    BISON=bison
fi
if [ "x$FLEX" == "x" ] ; then
    FLEX=flex
fi

# For convenience, change to the directory where the files are located
# if the script is run from the root of the source tree.
dirname=src/gromacs/selection
if [[ -f $dirname/parser.y && -f $dirname/scanner.l ]] ; then
    cd $dirname
fi

# We apply some trivial patches to the output to avoid warnings for PGI
# (and maybe other) compilers
[[ $FORCE || parser.y  -nt parser.cpp ]]  && \
    echo Generating parser.cpp and parser.h... && \
    $BISON -t -o parser.cpp --defines=parser.h parser.y && \
    patch -p0 < parser.patch && \
    rm -f parser.cpp.orig
[[ $FORCE || scanner.l -nt scanner.cpp ]] && \
    echo Generating scanner.cpp and scanner_flex.h... && \
    $FLEX -o scanner.cpp scanner.l && \
    patch -p0 < scanner.patch && \
    rm -f scanner.cpp.orig
