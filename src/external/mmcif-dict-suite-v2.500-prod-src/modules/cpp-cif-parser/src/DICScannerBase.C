//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file DICScannerBase.C
**
** \brief Implementation file for DICScanner class.
*/


/* 
  PURPOSE:    DDL 2.1 compliant CIF file lexer ...
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DICScannerBase.h"
#include "DICScannerInt.h"
#include "DICParserBase.h"
#include "DICParser.h"

extern int dicparser_leng;
extern char* dicparser_text;
extern YYSTYPE dicparser_lval;
extern FILE* dicparser_in;

#define yyleng dicparser_leng
#define yytext dicparser_text
#define yylval dicparser_lval
#define yyin dicparser_in

extern "C" void dic_yy_less(int i);

extern DICParser* DICParserP;


#if 0
int yywrap(void)
{
   return(1);
}
#endif

using std::ios;
using std::endl;

DICScanner::DICScanner() {
  Clear();
  _verbose=false;
   // VLAD: WATCH HERE
  _tBuf = new string;
  _tBuf->clear();
 if (_verbose) log << "Default constructor called" << endl;
}

void DICScanner::OpenLog(const string& logName, bool verboseLevel)
{
  _verbose = verboseLevel;
//  if (_verbose && logName) log.open(logName,ios::out|ios::trunc);
  if (!logName.empty()) log.open(logName.c_str(),ios::out|ios::trunc);
}

/*
DICScanner::DICScanner(istream *in) {
  Clear();
  _verbose=false;
  // VLAD: WATCH HERE
  _tBuf = new string(1025,512);
  _tBuf->Clear();
  yyin=in;
 if (_verbose) log << "DICScanner::DICScanner(istream *in, int verbose) constructor called" << endl;
}
*/

void DICScanner::Clear(void) {
//if (_tBuf) delete _tBuf;
  _tBuf=NULL;
  _isText = false;
  NDBlineNo=1;
}


void DICScanner::Reset(void) {
if (_tBuf) delete _tBuf;
  _tBuf=NULL;
  Clear();
}

int DICScanner::yylex()
{
   return(0);
}

int DICScanner::ProcessNone()
{
#if DEBUG
          log << "LS0: line "<<  NDBlineNo <<  " length " << yyleng << " yytext=" << yytext << endl;
#endif
          NDBlineNo++;
          if (_isText == true) {          /* end of text value */
             for (_i=yyleng-1; _i >= 0; _i--) {
               if ( yytext[_i] == ' ' || yytext[_i] == '\t' ||  yytext[_i] == '\n') {
                  yytext[_i]='\0';
               } else if ( yytext[_i] == ';') {
                    yytext[_i]='\0';
                    break;
               } else
                  break;
             }
             (*_tBuf)+=yytext;
//          _tBuf->InsertAt(strlen(_tBuf->Text())-1,"\0");
          _tBuf->erase(strlen(_tBuf->c_str())-1,1);
//printf("%d   \n",strlen(_tBuf->Text()),_tBuf->Text());
//cout<<strlen(_tBuf->Text());
//_tBuf->Print();
//cout<<endl;
             yylval.cBuf=(char*)_tBuf->c_str();
             _isText = false;
#if DEBUG
          log << "LS1: String[" <<  strlen(yylval.cBuf) << "] " << yylval.cBuf << endl;
#endif
             return(LSITEMVALUE_DIC);
          } else {  /* text value begins */
             _isText = true;
             for (_i=0; _i < yyleng; _i++) {
                 if (yytext[_i] == ';') {  break; }
             }
             _tBuf->clear();
             (*_tBuf) += &yytext[_i+1];

          }

          return(0);
}

void DICScanner::ProcessWhiteSpace()
{
         for (_i=0; _i < yyleng; _i++)
           if (yytext[_i] == '\n') NDBlineNo++;
    if (_isText)
        (*_tBuf) += yytext;

/*      if (_isText)
        (*_tBuf) += yytext;
       else {
         for (_i=0; _i < yyleng; _i++)
           if (yytext[_i] == '\n') NDBlineNo++;
       }*/
}

int DICScanner::ProcessData()
{
      if (_isText)
        (*_tBuf) += yytext;
      else {
        yylval.cBuf=yytext;
        return (DATABLOCK_DIC);
      }

      return(0);
}

int DICScanner::ProcessItemSaveBegin()
{

    if (_isText)
        (*_tBuf) += yytext;
    else {
#if DEBUG
   log<< "SF: line "<<NDBlineNo<<", Starting save frame "<<&yytext[5]<<endl;
#endif
       if (isSave)
          log<< "Syntax error line "<< NDBlineNo<<" with "<< &yytext[5]<<", end of save expected"<<endl;
       else {
          yylval.cBuf=yytext;
          isSave=2;
          return(SAVE_BEGIN_DIC);
       }
    }

    return(0);
}

int DICScanner::ProcessCategorySaveBegin()
{

   if (_isText)
      (*_tBuf) += yytext;
   else {
#if DEBUG
   log<< "SF: line "<<NDBlineNo<<" Starting save frame "<<&yytext[5]<<endl;
#endif
       if (isSave)
          log<< "Syntax error line "<< NDBlineNo<<" with "<< &yytext[5]<<", end of save expected"<<endl;
       else {
          yylval.cBuf=yytext;
          isSave=1;
          return(SAVE_BEGIN_DIC);
          }
       }

    return(0);
}

int DICScanner::ProcessSaveEndScanner()
{
    if (_isText)
        (*_tBuf) += yytext;
    else {
    #if DEBUG
    log<<"SF: line "<<NDBlineNo<<" Ending save frame"<<endl;
    #endif
       if (!isSave)
          log<< "Syntax error line "<< NDBlineNo<<" no open save frame "<<endl;
       else {
                    yylval.cBuf=yytext;
                    isSave=0;
          return(SAVE_END_DIC);
       }
    }

    return(0);
}

int DICScanner::ProcessLoopScanner()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (LOOP_DIC);

      return(0);
}

void DICScanner::ProcessStop()
{
      if (_isText)
        (*_tBuf) += yytext;
}

int DICScanner::ProcessDot()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (UNKNOWN_DIC);

      return(0);
}

int DICScanner::ProcessQuestion()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (MISSING_DIC);

      return(0);
}

void DICScanner::ProcessComment()
{
      if (_isText)
        (*_tBuf) += yytext;
/*      else
        NDBlineNo++;*/
}

int DICScanner::ProcessItemNameScanner()
{
      if (_isText) {
        (*_tBuf) += yytext;
      } else {
        /* If the beginning of text is in buffer yytext */
         for (_i=0; _i<yyleng; _i++)
           if (yytext[_i] == '_') break;
         yylval.cBuf=yytext;
         return(ITEMNAME_DIC);
      }

      return(0);
}

int DICScanner::ProcessUnquotedString()
{
     if (!_isText) {
        _j=0;
        /*
        ** The string is not part of a multiline CIF value, but a CIF value
        ** anywhere else, in a loop or out of the loop.
        */
/*        for (_i=yyleng-1; _i >= 0; _i--) {
            if ( yytext[_i] == '\'' || yytext[_i] == '\"') {
               yytext[_i]='\0';
               break;
            } else
             break;
        }
        for (_i=0; _i < yyleng; _i++) {
             if (yytext[_i] == '\'' || yytext[_i] == '\"') {
                _j++;
                break;
            } else
               break;
        }
*/
        yylval.cBuf=&yytext[_j];
#if DEBUG
        log << "UQ: String " << yylval.cBuf << endl;
#endif

#ifdef REPORT_EMBEDDED_QUOTES
        unsigned int cBufLen = strlen(yylval.cBuf);
        for (unsigned int i = 0; i < cBufLen; ++i)
        {
            if ((yylval.cBuf[i] == '\'') || (yylval.cBuf[i] == '\"'))
            {
                log << "ERROR - Invalid character at line " <<
                  String::IntToString(NDBlineNo) << " in CIF value " <<
                  yylval.cBuf << endl;
            }
        }
#endif // REPORT_EMBEDDED_QUOTES

        return(ITEMVALUE_DIC);
     }
     else {
        /*
        ** The string is part of a multiline CIF value. It is processed as is.
        */
#if DEBUG
          log << "UQx: String " << yytext<< endl;
#endif
        (*_tBuf) += yytext;
     }

     return(0);
}

int DICScanner::ProcessSQuotedString()
{
char * p;
     if (!_isText) {
        p=yytext;
                  p++;
        while ((p=strchr(p,'\''))) {
          p++;
          if ( p[0] == ' ' || p[0] == '\t' || p[0] == '\n') {
             _i=yyleng-strlen(p);
             dic_yy_less(_i);
             p=&yytext[yyleng];
          }
        }
        yylval.cBuf=&yytext[1];
        yylval.cBuf[_i-2]='\0';
#if DEBUG
//        ndb_log_message_text(NDB_MSG_DEBUG,"SQ String:%s",cifp_lval.TempBuffer);
#endif
        return(ITEMVALUE_DIC);
     }
     else {
        if (yytext[yyleng-1] == '\n') dic_yy_less(yyleng-1);
        (*_tBuf) += yytext;

     }

     return(0);
}

int DICScanner::ProcessDQuotedString()
{
char * p;
     if (!_isText) {
        p=yytext;
                  p++;
        while ((p=strchr(p,'\"'))) {
          p++;
          if ( p[0] == ' ' || p[0] == '\t' || p[0] == '\n') {
             _i=yyleng-strlen(p);
             dic_yy_less(_i);
             p=&yytext[yyleng];
          }
        }
        yylval.cBuf=&yytext[1];
        yylval.cBuf[_i-2]='\0';
#if DEBUG
//        ndb_log_message_text(NDB_MSG_DEBUG,"SQ String:%s",cifp_lval.TempBuffer);
#endif
        return(ITEMVALUE_DIC);
     }
     else {
        if (yytext[yyleng-1] == '\n') dic_yy_less(yyleng-1);
        (*_tBuf) += yytext;

     }

     return(0);
}

int DICScanner::ProcessEof()
{
   if (_isText == true) {
      _isText=false;
           log<<"String is not not finish at line "<<NDBlineNo<< endl;
           return(1);
        }
        else
           return(0);
}


int ProcessNoneFromDICScanner()
{
    return (DICParserP->ProcessNone());
}

void ProcessWhiteSpaceFromDICScanner()
{
    DICParserP->ProcessWhiteSpace();
}

int ProcessDataFromDICScanner()
{
    return (DICParserP->ProcessData());
}

int ProcessItemSaveBeginFromDICScanner()
{
    return (DICParserP->ProcessItemSaveBegin());
}

int ProcessCategorySaveBeginFromDICScanner()
{
    return (DICParserP->ProcessCategorySaveBegin());
}

int ProcessSaveEndFromDICScanner()
{
    return (DICParserP->ProcessSaveEndScanner());
}

int ProcessLoopFromDICScanner()
{
    return (DICParserP->ProcessLoopScanner());
}

void ProcessStopFromDICScanner()
{
    DICParserP->ProcessStop();
}

int ProcessDotFromDICScanner()
{
    return (DICParserP->ProcessDot());
}

int ProcessQuestionFromDICScanner()
{
    return (DICParserP->ProcessQuestion());
}

void ProcessCommentFromDICScanner()
{
    DICParserP->ProcessComment();
}

int ProcessItemNameFromDICScanner()
{
    return (DICParserP->ProcessItemNameScanner());
}

int ProcessUnquotedStringFromDICScanner()
{
    return (DICParserP->ProcessUnquotedString());
}

int ProcessSQuotedStringFromDICScanner()
{
    return (DICParserP->ProcessSQuotedString());
}

int ProcessDQuotedStringFromDICScanner()
{
    return (DICParserP->ProcessDQuotedString());
}

int ProcessEofFromDICScanner()
{
    return (DICParserP->ProcessEof());
}

