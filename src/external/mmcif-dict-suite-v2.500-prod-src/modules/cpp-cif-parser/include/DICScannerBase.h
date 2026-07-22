//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file DICScannerBase.h
**
** \brief Header file for DICScanner class.
*/


/* 
  PURPOSE:    DDL 2.1 compliant CIF file lexer ...
*/


#ifndef DICSCANNERBASE_H
#define DICSCANNERBASE_H


/*
#if !defined(DDL_FLEX_LEXER_INCLUDED)
#undef    yyFlexLexer
#define   yyFlexLexer DDLFlexLexer
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
** \class DICScanner
**
** \brief Private class that represents a dictionary scanner.
*/
class DICScanner // : public DDLFlexLexer 
{
 protected:

  std::string *_tBuf;
  int   _isText;
  int   _i, _j, _len;

 protected:
  std::ofstream log;
  bool _verbose;
  void alt_yymore(void);
  void OpenLog(const std::string& logName, bool verboseLevel);

 public:
  int NDBlineNo;
  int isSave;
  DICScanner(std::istream *yyin);
  DICScanner();
  void Clear();
  void Reset();
  int ProcessNone();
  void ProcessWhiteSpace();
  int ProcessData();
  int ProcessItemSaveBegin();
  int ProcessCategorySaveBegin();
  int ProcessSaveEndScanner();
  int ProcessLoopScanner();
  void ProcessStop();
  int ProcessDot();
  int ProcessQuestion();
  void ProcessComment();
  int ProcessItemNameScanner();
  int ProcessUnquotedString();
  int ProcessSQuotedString();
  int ProcessDQuotedString();
  int ProcessEof();
  virtual int yylex();
  virtual ~DICScanner() {Reset();};
};

#endif

