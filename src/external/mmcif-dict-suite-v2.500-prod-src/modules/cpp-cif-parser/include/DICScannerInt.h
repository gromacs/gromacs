/*!
** \file DICScannerInt.h
**
** \brief Header file for flex interfacing to DICScanner class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#ifndef DIC_SCANNER_INT_H
#define DIC_SCANNER_INT_H

#ifdef __cplusplus
extern "C" {
#endif

int ProcessNoneFromDICScanner();
void ProcessWhiteSpaceFromDICScanner();
int ProcessDataFromDICScanner();
int ProcessItemSaveBeginFromDICScanner();
int ProcessCategorySaveBeginFromDICScanner();
int ProcessSaveEndFromDICScanner();
int ProcessLoopFromDICScanner();
void ProcessStopFromDICScanner();
int ProcessDotFromDICScanner();
int ProcessQuestionFromDICScanner();
void ProcessCommentFromDICScanner();
int ProcessItemNameFromDICScanner();
int ProcessUnquotedStringFromDICScanner();
int ProcessSQuotedStringFromDICScanner();
int ProcessDQuotedStringFromDICScanner();
int ProcessEofFromDICScanner();
#ifdef __cplusplus
}
#endif
#endif /* DIC_SCANNER_BASE_H */
