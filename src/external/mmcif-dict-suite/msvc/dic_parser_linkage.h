/*
 * Force-included by MSVC into DICParserBase.C and DICScannerBase.C, for the
 * same reason as cif_parser_linkage.h.
 */
#ifndef GMX_MMCIF_MSVC_DIC_PARSER_LINKAGE_H
#define GMX_MMCIF_MSVC_DIC_PARSER_LINKAGE_H

#include <stdio.h>

extern "C"
{
#include "DICParser.h"

    extern FILE* dicparser_in;
    extern char* dicparser_text;
    extern int   dicparser_leng;

    extern char* Glob_tBufKeywordSaveDIC;
    extern char* Glob_pBufValueDIC;
    extern char* Glob_dataBlockNameDIC;
}

#endif
