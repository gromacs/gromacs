/*!
** \file DICParserInt.h
**
** \brief Header file for bison interfacing to DICParser class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#ifndef DIC_PARSER_INT_H
#define DIC_PARSER_INT_H

/*
#include <string>

*/

#ifdef __cplusplus
extern "C" {
#endif
void ProcessAssignmentsFromDICParser();
void ProcessOneAssignmentFromDICParser();
void ProcessItemNameListLoopFromDICParser();
void ProcessItemNameListNameFromDICParser();
void ProcessValueListFromDICParser();
void ProcessItemNameFromDICParser();
void ProcessLoopFromDICParser();
void ProcessItemValueFromDICParser();
void ProcessLsItemValueFromDICParser();
void ProcessUnknownValueFromDICParser();
void ProcessMissingValueFromDICParser();
void ProcessSaveBeginFromDICParser();
void ProcessSaveEndFromDICParser();
void ProcessDataBlockNameFromDICParser();
void dicparser_error(const char*);
#ifdef __cplusplus
}
#endif

extern char* Glob_tBufKeywordSaveDIC;
extern char* Glob_pBufValueDIC;
extern char* Glob_dataBlockNameDIC;

#endif /* DIC_PARSER_INT_H */
