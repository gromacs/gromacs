#!/bin/sh
#
# This script runs Bison and/or Flex to regenerate the files
# parser.c and parser.h from parser.y and scanner.c from scanner.l.
# The commands are run only if the generated files are older than the
# Bison/Flex input files.

[[ parser.y  -nt parser.c ]]  && bison -d -t -o parser.c parser.y
[[ scanner.l -nt scanner.c ]] && flex -t scanner.l >scanner.c
