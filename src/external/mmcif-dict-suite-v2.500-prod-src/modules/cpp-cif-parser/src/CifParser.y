%{
/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define  LVAL yylval

#include "CifParserInt.h"
/* #include "CifParserBase.h" */
/* #include "CifScanner.h" */
/* #include "CifFileObj.h" */
extern int cifparser_lex();
%}

%union {
	char *cBuf;
}
%token <cBuf> ITEMNAME_CIF
%token <cBuf> ITEMVALUE_CIF
%token <cBuf> LSITEMVALUE_CIF
%token <cBuf> LOOP_CIF
%token <cBuf> DATABLOCK_CIF
%token <cBuf> UNKNOWN_CIF
%token <cBuf> MISSING_CIF
%type  <cBuf> ItemName ItemValue 
%%

Datablocks :     Assignments
             |   Datablocks  Datablock 
             ;

Datablock :    DatablockName Assignments ;


Assignments:  /* empty */ 
             | Assignments error Assignment 
             | Assignments Assignment 
{
  ProcessAssignmentsFromParser();
} 
             ;

Assignment:       ItemName ItemValue  /* instance of a non-looped item name value pair */
{
  int ret;

  ret = ProcessItemValuePairFromParser();
  if (ret == STOP_PARSING)
  {
    YYABORT;
  }

}
  |   ItemNameList ValueList
  ;


ItemNameList:  

  Loop  ItemName           /*  the beginning of a loop_ */
{
  int ret;

  ret = ProcessLoopDeclarationFromParser();
  if (ret == STOP_PARSING)
  {
    YYABORT;
  }
}
 
  | ItemNameList ItemName 
{
  int ret;

  ret = ProcessItemNameListFromParser();
  if (ret == STOP_PARSING)
  {
    YYABORT;
  }

}
  ;


ValueList:  ItemValue
{
  int ret;

  ret = ProcessValueListFromParser();
  if (ret == STOP_PARSING)
  {
    YYABORT;
  }

}
  | ValueList ItemValue
{
  int ret;

  ret = ProcessValueListFromParser();
  if (ret == STOP_PARSING)
  {
    YYABORT;
  }

}
;

ItemName : ITEMNAME_CIF  
{
   Glob_tBufKeyword = $1;
   ProcessItemNameFromParser();
}
;

Loop :    LOOP_CIF
{
  ProcessLoopFromParser();
}
;

ItemValue:   ITEMVALUE_CIF 
{
   Glob_pBufValue = $1;
   ProcessItemValueFromParser();
}
        | LSITEMVALUE_CIF
{
   Glob_pBufValue = $1;
   ProcessLsItemValueFromParser();
}
        | UNKNOWN_CIF 
{
   ProcessUnknownValueFromParser();
}

        | MISSING_CIF
{
   ProcessMissingValueFromParser();
}
;

DatablockName: DATABLOCK_CIF
{
  Glob_dataBlockName = $1;
  ProcessDataBlockNameFromParser();
}
;

%%

