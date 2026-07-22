%{
/* 
  PURPOSE:    DDL 2.1 compliant DIC file lexer ...
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DICParser.h"
#include "DICScannerInt.h"

/*
#undef  YY_DECL
#define YY_DECL int DICParser::yylex()
#undef   DIC_FLEX_LEXER_INCLUDED
#define  DIC_FLEX_LEXER_INCLUDED
#include "DICParser.h"
#include "DICScanner.h"

*/
%}

/* %option noyywrap c++ yyclass="DICScanner" prefix="DIC" */
%option noinput
%option nounput

%%

^;[^\n]*[\n] { 
    int ret;

    ret = ProcessNoneFromDICScanner();
    if (ret == LSITEMVALUE_DIC)
        return ret;
       }



[\t \n]+      {                                    /* white space */ 
    ProcessWhiteSpaceFromDICScanner();
    }

[Dd][Aa][Tt][Aa][_][^\n\t ]*   {                 /* data_ */ 
    int ret;

    ret = ProcessDataFromDICScanner();
    if (ret == DATABLOCK_DIC)
        return ret;
    }


[Ss][Aa][Vv][Ee][_][_][^\n\t ]+ {             /*Save frame begins */
    int ret;

    ret = ProcessItemSaveBeginFromDICScanner();
    if (ret == SAVE_BEGIN_DIC)
        return ret;
   }


[Ss][Aa][Vv][Ee][_][^\n\t ]+ {             /*Save frame begins */
    int ret;

    ret = ProcessCategorySaveBeginFromDICScanner();
    if (ret == SAVE_BEGIN_DIC)
        return ret;
    }

[Ss][Aa][Vv][Ee][_]  {                /*Save frame ends */
    int ret;

    ret = ProcessSaveEndFromDICScanner();
    if (ret == SAVE_END_DIC)
        return ret;
   }


[Ll][Oo][Oo][Pp][_]      {                        /* loop_ */
    int ret;

    ret = ProcessLoopFromDICScanner();
    if (ret == LOOP_DIC)
        return ret;
    }

[Ss][Tt][Oo][Pp][_]      {                        /* stop_ */
    ProcessStopFromDICScanner();
    }

[.]  {                                        /* Unknown */
    int ret;

    ret = ProcessDotFromDICScanner();
    if (ret == UNKNOWN_DIC)
        return ret;
     }

[?]  {                                        /* Missing */
    int ret;

    ret = ProcessQuestionFromDICScanner();
    if (ret == MISSING_DIC)
        return ret;
    }

[#][^\n]*   {                             /* COMMENT */
    ProcessCommentFromDICScanner();
    }

[_][^\n \t]+  {                /* Item names start with _ end with spaces*/
    int ret;

    ret = ProcessItemNameFromDICScanner();
    if (ret == ITEMNAME_DIC)
        return ret;
   }

[^\n \t]+   {  /* unquoted or mixed quote strings strings */
    int ret;

    ret = ProcessUnquotedStringFromDICScanner();
    if (ret == ITEMVALUE_DIC)
        return ret;
   }

[\'](.)*[\'][ \t\n]  { /*  Single Quoted Strings */
    int ret;

    ret = ProcessSQuotedStringFromDICScanner();
    if (ret == ITEMVALUE_DIC)
        return ret;
   }


[\"](.)*[\"][ \t\n]        {  /* Double Quoted Strings */
    int ret;

    ret = ProcessDQuotedStringFromDICScanner();
    if (ret == ITEMVALUE_DIC)
        return ret;
   }


<<EOF>> {
    return (ProcessEofFromDICScanner());
}

%%

void dic_yy_less(int i)
{
    yyless(i);
}

int yywrap(void)
{
    return (1);
}
