%{
/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define  LVAL yylval

#include "DICParserInt.h"
extern int dicparser_lex();
extern void dicparser_error(const char*);
%}

%union {
	char *cBuf;
}
%token <cBuf> ITEMNAME_DIC
%token <cBuf> ITEMVALUE_DIC
%token <cBuf> LSITEMVALUE_DIC
%token <cBuf> LOOP_DIC
%token <cBuf> DATABLOCK_DIC
%token <cBuf> UNKNOWN_DIC
%token <cBuf> MISSING_DIC
%token <cBuf> SAVE_BEGIN_DIC
%token <cBuf> SAVE_END_DIC
%type  <cBuf> ItemName ItemValue 

%%

Datablocks :     Sections     /*   Assignments */
             |   Datablocks  Datablock 
/*
             |   Datablocks error  Datablock 
{yyclearin;}
*/
             ;

Datablock :  DatablockName  Sections;

Sections:
             | Sections Section

             ;
Section:     Assignment
{}
	     | SaveBegin Assignments SaveEnd
{}
;

Assignments:  /* empty */ 
             | Assignments error Assignment 
{ /* yyclearin; */}
             | Assignments Assignment 
{
    ProcessAssignmentsFromDICParser();
} 
             ;

Assignment:       ItemName ItemValue  /* instance of a non-looped item name value pair */
{
    ProcessOneAssignmentFromDICParser();
}
  |   ItemNameList ValueList
  ;


ItemNameList:  

  Loop  ItemName           /*  the beginning of a loop_ */
{
    ProcessItemNameListLoopFromDICParser();
}
 
  | ItemNameList ItemName 
{
    ProcessItemNameListNameFromDICParser();
}
  ;


ValueList:  ItemValue 
{
    ProcessValueListFromDICParser();
}
  | ValueList ItemValue
{
    ProcessValueListFromDICParser();
}
;

ItemName : ITEMNAME_DIC  
{
   Glob_tBufKeywordSaveDIC = $1;
   ProcessItemNameFromDICParser();
}
;

Loop :    LOOP_DIC
{
    ProcessLoopFromDICParser();
}
;

ItemValue:   ITEMVALUE_DIC 
{
   Glob_pBufValueDIC = $1;
   ProcessItemValueFromDICParser();
}
        | LSITEMVALUE_DIC
{
   Glob_pBufValueDIC = $1;
   ProcessLsItemValueFromDICParser();
}
        | UNKNOWN_DIC 
{
    ProcessUnknownValueFromDICParser();
}

        | MISSING_DIC
{
    ProcessMissingValueFromDICParser();
}
;

SaveBegin: SAVE_BEGIN_DIC
{
    Glob_dataBlockNameDIC = $1;
    ProcessSaveBeginFromDICParser();
}
;
SaveEnd:   SAVE_END_DIC
{
  Glob_pBufValueDIC = $1;
  ProcessSaveEndFromDICParser();
}
;

DatablockName: DATABLOCK_DIC
{
  Glob_dataBlockNameDIC = $1;
  ProcessDataBlockNameFromDICParser();
}
;

%%

