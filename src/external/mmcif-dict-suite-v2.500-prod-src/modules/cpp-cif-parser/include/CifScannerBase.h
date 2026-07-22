//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifScannerBase.h
**
** \brief Header file for CifScanner class.
*/


/* 
  PURPOSE:    DDL 2.1 compliant CIF file lexer ...
*/


#ifndef CIFSCANNERBASE_H
#define CIFSCANNERBASE_H


/*
#if !defined(FLEX_LEXER_INCLUDED)
#undef    yyFlexLexer
#define   yyFlexLexer CifFlexLexer
#include "FlexLexer.h"
#endif
*/

#include <string>


#include <fstream>
#include <stdio.h>
#include <string.h>

#ifndef  DEBUG
#define DEBUG  0
#endif


/**
** \class CifScanner
**
** \brief Private class that represents a CIF scanner.
*/
class CifScanner // : public CifFlexLexer 
{
 protected:

  std::string *_tBuf;
  int   _isText;
  int   _i, _j, _len;

 protected:
  std::ofstream log;
  std::string errorLog;
  bool _verbose;
  void alt_yymore(void);
  void OpenLog(const std::string& logName, bool verboseLevel);

 public:
  int NDBlineNo;
  CifScanner(std::istream *yyin);
  CifScanner();
  int ProcessNone();
  void ProcessWhiteSpace();
  int ProcessData();
  int ProcessLoopScanner();
  void ProcessStop();
  int ProcessDot();
  int ProcessQuestion();
  void ProcessComment();
  int ProcessUnderscore();
  int ProcessBadStrings();
  int ProcessSQuotedStrings();
  int ProcessDQuotedStrings();
  int ProcessEof();
  void Clear();
  void Reset();
  virtual int yylex();
  virtual ~CifScanner() {Reset();};
};

#endif /* CIFSCANNERBASE_H */
