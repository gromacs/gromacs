%{ 
/* 
  PURPOSE:    DDL 2.1 compliant CIF file lexer ...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CifScannerInt.h"
#include "CifParser.h"

/*
#undef  YY_DECL
#define YY_DECL int CifParser::yylex()
#undef   FLEX_LEXER_INCLUDED
#define  FLEX_LEXER_INCLUDED
#include "CifParser.h"
#include "CifScanner.h"

*/

%}

/* %option noyywrap c++ yyclass="CifScanner" prefix="Cif" */
%option noinput
%option nounput

CR      \r
LF      \n
EOL     (({CR}{LF})|{LF}|{CR}|({LF}{CR}))

%%

^;[^\n\r]*{EOL} { 
    int ret;

    ret = ProcessNoneFromScanner();
    if (ret == LSITEMVALUE_CIF)
        return ret;
       }

([ \t]|{EOL})+      {                                    /* white space */ 
    ProcessWhiteSpaceFromScanner();

    }

[Dd][Aa][Tt][Aa][_][^\n\r\t ]*   {                 /* data_ */ 
    int ret;

    ret = ProcessDataFromScanner();
    if (ret == DATABLOCK_CIF)
        return ret;
    }

[Ll][Oo][Oo][Pp][_]      {                        /* loop_ */
    int ret;

    ret = ProcessLoopFromScanner();
    if (ret == LOOP_CIF)
        return ret;
    }

[Ss][Tt][Oo][Pp][_]      {                        /* stop_ */
    ProcessStopFromScanner();
    }

[.]  {                                        /* Unknown */
    int ret;

    ret = ProcessDotFromScanner();
    if (ret == UNKNOWN_CIF)
        return ret;
     }

[?]  {                                        /* Missing */
    int ret;

    ret = ProcessQuestionFromScanner();
    if (ret == MISSING_CIF)
        return ret;
    }

[#][^\n\r]*   {                             /* COMMENT */
    ProcessCommentFromScanner();
    }

[_][^\n\r \t]+  {                /* Item names start with _ end with spaces*/
    int ret;

    ret = ProcessUnderscoreFromScanner();
    if (ret == ITEMNAME_CIF)
        return ret;
   }

[^\n\r \t]+   {  /* unquoted or mixed quote strings strings */
    int ret;

    ret = ProcessBadStringsFromScanner();
    if (ret == ITEMVALUE_CIF)
        return ret;
   }

[\'](.)*[\']([ \t]|{EOL})  { /*  Single Quoted Strings */
    int ret;

    ret = ProcessSQuotedStringsFromScanner();
    if (ret == ITEMVALUE_CIF)
        return ret;
   }


[\"](.)*[\"]([ \t]|{EOL})  {  /* Double Quoted Strings */
    int ret;

    ret = ProcessDQuotedStringsFromScanner();
    if (ret == ITEMVALUE_CIF)
        return ret;
   }


<<EOF>> {
    return (ProcessEofFromScanner());
}

%%

void cif_yy_less(int i)
{
    yyless(i);
}

int yywrap(void)
{
    return (1);
}
