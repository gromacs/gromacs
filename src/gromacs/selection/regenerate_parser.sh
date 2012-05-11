#!/bin/bash
#
# This script runs Bison and/or Flex to regenerate the files as follows:
#   parser.y  -> parser.cpp, parser.h
#   scanner.l -> scanner.cpp, scanner_flex.h
# The commands are run only if the generated files are older than the
# Bison/Flex input files, or if a '-f' flag is provided.

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

[[ $FORCE || parser.y  -nt parser.cpp ]]  && $BISON -t -o parser.cpp --defines=parser.h parser.y
[[ $FORCE || scanner.l -nt scanner.cpp ]] && $FLEX -o scanner.cpp scanner.l
