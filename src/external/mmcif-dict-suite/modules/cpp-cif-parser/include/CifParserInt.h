/*!
** \file CifParserInt.h
**
** \brief Header file for bison interfacing to CifParser class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#ifndef CIF_PARSER_INT_H
#define CIF_PARSER_INT_H

#define STOP_PARSING 2

#ifdef __cplusplus
extern "C" {
#endif
void ProcessAssignmentsFromParser();
int ProcessItemValuePairFromParser();
int ProcessLoopDeclarationFromParser();
void ProcessLoopFromParser();
int ProcessItemNameListFromParser();
int ProcessValueListFromParser();
void ProcessItemNameFromParser();
void ProcessItemValueFromParser();
void ProcessLsItemValueFromParser();
void ProcessUnknownValueFromParser();
void ProcessMissingValueFromParser();
void ProcessDataBlockNameFromParser();
void cifparser_error(const char*);
#ifdef __cplusplus
}
#endif
extern char* Glob_tBufKeyword;
extern char* Glob_pBufValue;
extern char* Glob_dataBlockName;

#endif /* CIF_PARSER_BASE_H */
