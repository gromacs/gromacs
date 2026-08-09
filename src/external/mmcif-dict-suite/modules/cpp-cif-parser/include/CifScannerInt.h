/*!
** \file CifScannerInt.h
**
** \brief Header file for flex interfacing to CifScanner class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#ifndef CIF_SCANNER_INT_H
#define CIF_SCANNER_INT_H

#ifdef __cplusplus
extern "C" {
#endif

int ProcessNoneFromScanner();
void ProcessWhiteSpaceFromScanner();
int ProcessDataFromScanner();
int ProcessLoopFromScanner();
void ProcessStopFromScanner();
int ProcessDotFromScanner();
int ProcessQuestionFromScanner();
void ProcessCommentFromScanner();
int ProcessUnderscoreFromScanner();
int ProcessBadStringsFromScanner();
int ProcessSQuotedStringsFromScanner();
int ProcessDQuotedStringsFromScanner();
int ProcessEofFromScanner();

#ifdef __cplusplus
}
#endif
#endif /* CIF_SCANNER_BASE_H */
