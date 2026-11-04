/*
 * Force-included by MSVC into CifParserBase.C and CifScannerBase.C.
 *
 * Those files share variables with the C parser and scanner in generated/, but
 * declare or define them without extern "C". Other compilers give C++ variables
 * the same symbol names as C ones, so this goes unnoticed there; MSVC decorates
 * them, and the link fails. Declaring them extern "C" first gives the later
 * declarations and definitions in the vendored files C linkage.
 */
#ifndef GMX_MMCIF_MSVC_CIF_PARSER_LINKAGE_H
#define GMX_MMCIF_MSVC_CIF_PARSER_LINKAGE_H

#include <stdio.h>

extern "C"
{
#include "CifParser.h"

    extern FILE* cifparser_in;
    extern char* cifparser_text;
    extern int   cifparser_leng;

    extern char* Glob_tBufKeyword;
    extern char* Glob_pBufValue;
    extern char* Glob_dataBlockName;
}

#endif
