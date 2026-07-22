//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file DICParserBase.C
**
** \brief Implementation file for DICParser class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#include <stdexcept>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Exceptions.h"
#include "GenString.h"
#include "RcsbFile.h"
#include "CifString.h"
#include "DICScannerBase.h"
#include "DICParserInt.h"
#include "DICParserBase.h"

extern FILE* dicparser_in;

extern "C"
{
  int dicparser_parse();
  void dicparser_restart(FILE*);
}

char* Glob_tBufKeywordSaveDIC;
char* Glob_pBufValueDIC;
char* Glob_dataBlockNameDIC;

DICParser* DICParserP = NULL;

using std::exception;
using std::endl;
#ifdef VLAD_DEBUG
using std::cout;
#endif

DICParser::DICParser(DicFile* fo, CifFile* ddl_in, bool verbose)
{

    if (DICParserP != NULL)
    {
        // Attempting to create a new parser, during the lifetime of
        // an existing parser.
        throw AlreadyExistsException("Cannot create a new parser, since "\
          "one already exists.", "DICParser::DICParser");
    }

    Clear();

    if (fo != NULL)
        _fobj = fo;
    else
        throw EmptyValueException("fo is a NULL pointer",
          "DICParser::DICParser");

    if (ddl_in != NULL)
        ddl = ddl_in;
    else
        throw EmptyValueException("ddl_in is a NULL pointer",
          "DICParser::DICParser");

    _verbose=verbose;

    format=fo->GetFormatTable();
    errorLog.clear();
    listcat.clear(); listcat.push_back("id");
    string BlockName;
    BlockName = ddl->GetFirstBlockName();

    Block& block = ddl->GetBlock(BlockName);
    cattbl = block.GetTablePtr("category");

    const vector<string>& cattblCols = cattbl->GetColumnNames();

    cattbl->SetFlags(cattblCols[0], ISTable::DT_STRING | ISTable::CASE_INSENSE);
    cattbl->SetFlags(cattblCols[1], ISTable::DT_STRING | ISTable::CASE_INSENSE);
    cattbl->SetFlags(cattblCols[2], ISTable::DT_STRING | ISTable::CASE_INSENSE);
    cattbl->CreateIndex("index0",listcat);

    listitem.clear(); listitem.push_back("category_id");
    listitem.push_back("name");
    BlockName = ddl->GetFirstBlockName();
    itemtbl = block.GetTablePtr("item");
    if (!itemtbl->IsColumnPresent("category_id"))
    {
        itemtbl->AddColumn("category_id");
    }

    pdbxitemtbl = block.GetTablePtr("pdbx_item");
    if (pdbxitemtbl != NULL)
    {
        if (!pdbxitemtbl->IsColumnPresent("category_id"))
        {
            pdbxitemtbl->AddColumn("category_id");
        }
    }

    itemtbl->CreateIndex("index0",listitem);
    if (pdbxitemtbl != NULL)
    {
        pdbxitemtbl->CreateIndex("index0",listitem);
    }

    listitem2.clear(); listitem2.push_back("category_id"); listitem2.push_back("mandatory_code");
    itemtbl->CreateIndex("index2",listitem2);
    if (pdbxitemtbl != NULL)
    {
        pdbxitemtbl->CreateIndex("index2",listitem2);
    }

    DICParserP = this;
}

void DICParser::Parse(const string& fileName, string& diagnostics)
{
    diagnostics.clear();

    FILE* dicIn;

    if ((dicIn = fopen(fileName.c_str(), "r")) == NULL )
    {
        diagnostics = "Unable to open file.";
        throw NotFoundException("File \"" + fileName + "\" cannot be opened",
          "DICParser::Parse");
    }

    string logFileName;
    RcsbFile::RelativeFileName(logFileName, fileName);

    logFileName += "-parser.log";

    OpenLog(logFileName, _verbose);

    dicparser_in = dicIn;

    dicparser_restart(dicparser_in);

    int ret = dicparser_parse();
    if (ret != 0)
    {
        int b = 0;
        b++;
    }

    fclose(dicIn);

    AfterParseProcessing();

    if (RcsbFile::IsEmpty(log))
    {
        log.close();
        RcsbFile::Delete(logFileName);
    }
    else
    {
        log.close();
    }

    if (this->errorLog.size() > 0)
    {
        diagnostics = this->errorLog;
    }

}

void dicparser_error(const char *s)
{
    DICParserP->Error(s);
}

DICParser::~DICParser()
{
    Reset();
    DICParserP = NULL;
}



void DICParser::Error(const char *s)
/*
 * Purpose:  yyerror() Print errors in DICScanner log.
 */
{
  errorLog += s;
  errorLog += " near line ";
  errorLog += String::IntToString(NDBlineNo);
  errorLog += '\n';
  /*if (_verbose) */
  log << s << " near line " <<  NDBlineNo << endl;
}

void DICParser::Reset()
{
  if (_tbl && (_curItemNo > 0)) { // write the current table / management of _tbl by _fobj
    CheckDDL();
    if (_nTablesInBlock) {
      Block& block = _fobj->GetBlock(_curDataBlockName);
      block.WriteTable(_tbl);
    } else {
      Block& block = _fobj->GetBlock(_prevDataBlockName);
      block.WriteTable(_tbl);
    }
    _nTablesInBlock++;
    //    delete _tbl;
  }
  _fieldList.clear();
  _fieldListSave.clear();
  if (_prevtbl) delete _prevtbl;
  if (_savetbl != NULL)
  {
      delete _savetbl;
      _savetbl = NULL;
  }
}

void DICParser::Clear()
{
  _nTablesInBlock=0;
  _numDataBlocks=0;
  _afterLoop = false;
  _curItemNo=0;
  _curRow=0;
  _curValueNo=0;
  _fieldListAlloc=100;
  _fieldList.reserve(_fieldListAlloc);
  _pBufValue.clear();
  _fobj=NULL;
  _tbl=NULL;
  isSave = 0;
  _nTablesInBlockSave=0;
  _curItemNoSave=0;
  _saveobj=NULL;
  _savetbl=NULL;
  _prevtbl=NULL;
  _fieldListAllocSave=100;
  _fieldListSave.reserve(_fieldListAllocSave);
  _tBufKeyword.clear();
  _curCategoryName.clear();
  _curDataBlockName.clear();
  _prevDataBlockName.clear();
  _curCategoryNameSave.clear();
  _curDataBlockName = "MISSING_DIC";
  _prevDataBlockName = "MISSING_DIC";
  // jdw xxx
  _curDataBlockNameSave.clear();
  _tmpDataBlockNameSave.clear();
  vector<string> list;
}


void DICParser::ProcessLoopDeclaration(void)
/* ----------------------------------------------------------------------
     Purpose: DICParser::ProcessLoopDeclaration(void)

              Handles initialization for a new loop, by creating a new
              category and adding the current item name to this category.
   ---------------------------------------------------------------------- */
{
  string categoryName;

#if DEBUG
  if (_verbose) log << "Processing loop declaration at line " << NDBlineNo << " value " <<_tBufKeyword << endl;
#endif
  _afterLoop = true;
  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Error in category name at line " << NDBlineNo << " value " << _tBufKeyword << endl;
    errorLog += "Error in category name at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
  }

  Block& block = _fobj->GetBlock(_curDataBlockName);

  if (block.IsTablePresent(categoryName)) {
    if (!isSave) {
      log << "Warning - Duplicate category name " << categoryName << " at line " << NDBlineNo << endl;
/*
      errorLog += "Duplicate category name ";
      errorLog += categoryName;
      errorLog += " at line ";
      errorLog += String::IntToString(NDBlineNo);
      errorLog += '\n';
*/
    }
  /*  else*/ {
      if (_tbl) { // write the current table / management of _tbl by _fobj
  CheckDDL();
        if (_nTablesInBlock) {
          Block& block = _fobj->GetBlock(_curDataBlockName);
          block.WriteTable(_tbl);
        } else {
          Block& block = _fobj->GetBlock(_prevDataBlockName);
          block.WriteTable(_tbl);
        }
        _nTablesInBlock++;
        //      delete _tbl;
      }

      Block& block = _fobj->GetBlock(_curDataBlockName);

      _tbl = block.GetTablePtr(categoryName);

      if (isSave ==0) {
        if (_prevtbl) delete _prevtbl;
        _prevtbl=new ISTable();
        *_prevtbl = *_tbl;
      }

      _curRow = _tbl->GetNumRows();
    ProcessItemNameList();
      _curCategoryName = categoryName;
      _tmpDataBlockNameSave = _curDataBlockNameSave;
    }

  }
else {
    if (_tbl) { // write the current table / management of _tbl by _fobj
  CheckDDL();
      if (_nTablesInBlock) {
        Block& block = _fobj->GetBlock(_curDataBlockName);
        block.WriteTable(_tbl);
        } else {
          Block& block = _fobj->GetBlock(_prevDataBlockName);
          block.WriteTable(_tbl);
        }
      _nTablesInBlock++;
      //      delete _tbl;
    }
    _tbl = new ISTable(categoryName);

    if (isSave ==0) {
      if (_prevtbl) delete _prevtbl;
      _prevtbl=new ISTable();
      *_prevtbl = *_tbl;
    }

    _curRow = 0;
    _curCategoryName = categoryName;
    _tmpDataBlockNameSave = _curDataBlockNameSave;

    ProcessItemNameList();
  }
}


void DICParser::ProcessItemNameList(void)
/* ----------------------------------------------------------------------
   Purpose: DICParser::ProcessItemNameList(void)

            Registers the item keyword for the the current item in the
            current category.  Maintains an index array of "valid" keyword
            names in fieldList[].  This array is used to indirectly
            reference between keywords and values ...
 * ----------------------------------------------------------------------*/
{

  string keywordName;
  string categoryName;

#if DEBUG
  if (_verbose) log << "Processing item name list at line " <<  NDBlineNo << " keyword " <<  _tBufKeyword << endl;
#endif

  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Error in category name at line " << NDBlineNo << " value " << _tBufKeyword << endl;
    errorLog += "Error in category name at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
  }

  if ( _curItemNo > _fieldListAlloc - 1) {
    _fieldListAlloc = _curItemNo + _fieldListAlloc;
    _fieldList.reserve(_fieldListAlloc);
  }

  vector<string> colNames, target, target2;
  colNames.push_back("dbName");
  colNames.push_back("type");
  colNames.push_back("catName");
  target.push_back(_curDataBlockName);
  if (isSave == 0) {
    target.push_back("data");
  }
  if (isSave == 1) {
    target.push_back("category");
  }
  if (isSave == 2) {
    target.push_back("item");
  }
  target.push_back(categoryName);
  unsigned int resIndex = format->FindFirst(target, colNames);
  if (resIndex == format->GetNumRows()) {
    format->AddRow(target);
  }

    if (isSave) {
  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if ((categoryName == _curCategoryName) && !keywordName.empty()) {
#if DEBUG
    if (_verbose) log << "*Line " << NDBlineNo << " keyword is " << keywordName << endl;
#endif
    if (!_tbl->IsColumnPresent(keywordName))
        _tbl->AddColumn(keywordName);

    if (_curItemNo >= (int)_fieldList.size())
      _fieldList.push_back(keywordName);
    else
      _fieldList[_curItemNo] = keywordName;
  }
  else {
    CifString::GetItemFromCifItem(keywordName, _tBufKeyword);

    if (!_tbl->IsColumnPresent(keywordName))
        _tbl->AddColumn(keywordName);

    if (_curItemNo >= (int)_fieldList.size())
      _fieldList.push_back(keywordName);
    else
      _fieldList[_curItemNo] = keywordName;
  }
    }
    else {
      _curCategoryName = categoryName;
  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if ((categoryName == _curCategoryName) && !keywordName.empty()) {
    if (_tbl->IsColumnPresent(keywordName)) {
      log << "Warning - Duplicate item name " <<_tBufKeyword  << " line " << NDBlineNo << endl;
/*
      errorLog += "Duplicate item name ";
      errorLog += _tBufKeyword;
      errorLog += " at line ";
      errorLog += String::IntToString(NDBlineNo);
      errorLog += '\n';
      return;
*/
    }
 /*   else */{
#if DEBUG
    if (_verbose) log << "Line " << NDBlineNo << " keyword is " << keywordName << endl;
#endif
    _tbl->AddColumn(keywordName);
    if (_curItemNo >= (int)_fieldList.size())
      _fieldList.push_back(keywordName);
    else
      _fieldList[_curItemNo] = keywordName;
    }
  } else {
    log << "*Syntax error at line " << NDBlineNo << " at item " << _tBufKeyword << endl;
    errorLog += "Syntax error at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _fieldList[_curItemNo].clear();
  }
    }
  _curItemNo++;

}
void DICParser::ProcessValueList(void)
/* ----------------------------------------------------------------------
     Purpose:  DICParser::ProcessValueList(void)

               Add the current value to the appropriate column in the
               the current row.  Start a new row if necessary.
 * ----------------------------------------------------------------------*/
{
#if DEBUG
  if (_verbose) {
    if (!_pBufValue.empty())
      log << "Processing value at line " << NDBlineNo << " value [" << _pBufValue << "]" << endl;
  }
#endif

  if (!_fieldList[_curValueNo].empty()) {

    if (_curValueNo == 0) {
      vector<string> rowBuf(_tbl->GetNumColumns(), CifString::UnknownValue);
      _curRow++;
      _tbl->AddRow(rowBuf);
    }

    if (!_pBufValue.empty()) {
      if (_pBufValue != CifString::InapplicableValue)
      {
        if (!_fieldList[_curValueNo].empty())
        {
#ifdef VLAD_DELETED
          try
#endif
          {
          _tbl->UpdateCell(_curRow-1, _fieldList[_curValueNo], _pBufValue);
          }
#ifdef VLAD_DELETED
          catch (...)
#endif
          {
          }
        }
      }
      else
      {
        if (!_fieldList[_curValueNo].empty())
        {
#ifdef VLAD_DELETED
          try
#endif
          {
          _tbl->UpdateCell(_curRow-1, _fieldList[_curValueNo], CifString::InapplicableValue);
          }
#ifdef VLAD_DELETED
          catch (...)
#endif
          {
          }
        }
      }
      } else {
        if (!_fieldList[_curValueNo].empty())
           _tbl->UpdateCell(_curRow-1, _fieldList[_curValueNo], CifString::UnknownValue);
      }
  }
  _curValueNo++;

  if (_curValueNo == _curItemNo) {
#if DEBUG
  if (_verbose) {
    log << "Loading row " << _curRow -1 << " with " <<  rowBuf.size() << " elements" << endl;
    for (int i=0; i < rowBuf.size(); i++) {
      log << "Column [" << i << "] value "<< rowBuf[i] << endl;
    }
  }
#endif
    _curValueNo = 0;

  }
  
}





void DICParser::ProcessItemValuePair(void)
/* ----------------------------------------------------------------------
      Purpose: ndb_cif_process_item_name_value_pair()

               Assign the current value to its associated item name.
 * ----------------------------------------------------------------------*/
{
  string categoryName;
  string keywordName;
  _curItemNo  = 1;
  _curValueNo = 0;
#if DEBUG
  if (_verbose)  {
    if (!_pBufValue.empty())
      log << "Processing " << _tBufKeyword << " at " <<  NDBlineNo << " value " << _pBufValue  << endl;
  }
#endif

  if ((_tBufKeyword == "_item.name") && (_pBufValue != _curDataBlockNameSave))
  {
      log << "ERROR - " << "In save frame \"save_" << _curDataBlockNameSave <<
        "\", \"_item.name\" has value \"" << _pBufValue << "\" at line " <<
        NDBlineNo << endl;
  }

  try
  {
    CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  }
  catch (const exception& exc)
  {

  }

  if (categoryName.empty())
  {
    log << "Error in category name at line " << NDBlineNo << " value " << _tBufKeyword << endl;
    errorLog += "Error in category name at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    return;
  }
  vector<string> colNames, target;
  colNames.push_back("dbName");
  colNames.push_back("type");
  colNames.push_back("catName");
  target.push_back(_curDataBlockName);
  if (isSave == 0) {
    target.push_back("data");
  }
  if (isSave == 1) {
    target.push_back("category");
  }
  if (isSave == 2) {
    target.push_back("item");
  }
  target.push_back(categoryName);
  unsigned int resIndex = format->FindFirst(target, colNames);
  if (resIndex == format->GetNumRows()) {
    format->AddRow(target);
  }

  Block& block = _fobj->GetBlock(_curDataBlockName);

  if (block.IsTablePresent(categoryName)) { //  duplicates a persistent table?
    if (!isSave) {
      log << "Warning - Duplicate category name " << categoryName << " at line " << NDBlineNo << endl;
/*
      errorLog += "Duplicate category name ";
      errorLog += categoryName;
      errorLog += " at line ";
      errorLog += String::IntToString(NDBlineNo);
      errorLog += '\n';
      return;
*/
    }
   /* else*/ {
  if ((categoryName != _curCategoryName)|| _afterLoop) {
    if (_tbl) { // write the current table / management of _tbl by _fobj
  CheckDDL();
      if (_nTablesInBlock) {
        Block& block = _fobj->GetBlock(_curDataBlockName);
        block.WriteTable(_tbl);
      } else {
        Block& block = _fobj->GetBlock(_prevDataBlockName);
        block.WriteTable(_tbl);
      }
      _nTablesInBlock++;
      //      delete _tbl;
    }
/*
  if ((categoryName == _curCategoryName)&&_afterLoop) { //  duplicates a persistent table?
    log << "Duplicate category name " << categoryName << " at line " << NDBlineNo << endl;
    errorLog += "Duplicate category name ";
    errorLog += categoryName;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    return;
  }
*/
      Block& block = _fobj->GetBlock(_curDataBlockName);

      _tbl = block.GetTablePtr(categoryName);

      if (isSave ==0) {
        if (_prevtbl) delete _prevtbl;
        _prevtbl=new ISTable();
        *_prevtbl = *_tbl;
      }

      _tbl->AddRow();
      _curCategoryName = categoryName;
      _tmpDataBlockNameSave = _curDataBlockNameSave;
  }
  _afterLoop=false;

  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);

  _curRow = _tbl->GetNumRows();
   if (!_tbl->IsColumnPresent(keywordName)) {
    {
      // Created a block to be able to re-use emptyColumn
    vector<string> emptyColumn;
    for (int i=0; i<_curRow; i++)
      emptyColumn.push_back(CifString::UnknownValue);
    _tbl->AddColumn(keywordName);
    _tbl->FillColumn(keywordName, emptyColumn);
    emptyColumn.clear();
    }
   }

   if (!_pBufValue.empty()) {
     _tbl->UpdateCell(_curRow-1, keywordName, _pBufValue);
   } else {
     _tbl->UpdateCell(_curRow-1, keywordName, CifString::UnknownValue);
   }
  
    }
  }
  else {
  if (categoryName != _curCategoryName) {
    if (_tbl) { // write the current table / management of _tbl by _fobj
  CheckDDL();
      if (_nTablesInBlock) {
        Block& block = _fobj->GetBlock(_curDataBlockName);
        block.WriteTable(_tbl);
      } else {
        Block& block = _fobj->GetBlock(_prevDataBlockName);
        block.WriteTable(_tbl);
      }
      _nTablesInBlock++;
      //      delete _tbl;
    }
    _tbl = new ISTable(categoryName);

    if (isSave ==0) {
      if (_prevtbl) delete _prevtbl;
      _prevtbl=new ISTable();
      *_prevtbl = *_tbl;
    }

    _curCategoryName = categoryName;
    _tmpDataBlockNameSave = _curDataBlockNameSave;

  }

  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if (!keywordName.empty()) {
#if DEBUG
    if (_verbose) log << "Line " << NDBlineNo << " keyword is " << keywordName << endl;
#endif
    _tbl->AddColumn(keywordName);
  }
  else {
    log << "Syntax2 error line " << NDBlineNo << " at item " << _tBufKeyword << endl;
    errorLog += "Syntax2 error at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
  }

  if (!_pBufValue.empty())
  {
    if (_pBufValue != CifString::InapplicableValue)
    {
    if (_tbl->GetNumRows() == 0)
        _tbl->AddRow();
    
    _tbl->UpdateCell(_tbl->GetNumRows() - 1, keywordName, _pBufValue);
    }
    else
    {
    if (_tbl->GetNumRows() == 0)
        _tbl->AddRow();
    
    _tbl->UpdateCell(_tbl->GetNumRows() - 1, keywordName, CifString::InapplicableValue);
    }
  }
  else
  {
    if (_tbl->GetNumRows() == 0)
        _tbl->AddRow();
    
    _tbl->UpdateCell(_tbl->GetNumRows() - 1, keywordName, CifString::UnknownValue);
  }

  }

}


//******************************

void DICParser::ProcessLoopDeclarationSave(void)
/* ----------------------------------------------------------------------
     Purpose: DICParser::ProcessLoopDeclarationSave(void)

              Handles initialization for a new loop, by creating a new
              category and adding the current item name to this category.
   ---------------------------------------------------------------------- */
{
  string categoryName;

  _curItemNoSave = 0;  _curValueNoSave = 0; 
  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Error in category name at line (save frame) " << NDBlineNo << " value " << _tBufKeyword << endl;
    errorLog += "Error in category name at line (save frame) ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
  }

  bool tablePresent = false;

  if (_saveobj->IsBlockPresent(_curDataBlockNameSave))
  {
      Block& block = _saveobj->GetBlock(_curDataBlockNameSave);
      tablePresent = block.IsTablePresent(categoryName);
  }

  if (tablePresent) { 
    log << "Duplicate category name in a save frame " << categoryName << " at line " << NDBlineNo << endl;
    errorLog += "Duplicate category name in a save frame ";
    errorLog += categoryName;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
  } else {
    if (_savetbl) { // write the current table /
      if (_nTablesInBlockSave) {
      if (_saveobj->IsBlockPresent(_curDataBlockNameSave))
      {
      Block& block = _saveobj->GetBlock(_curDataBlockNameSave);
      block.WriteTable(_savetbl);
      }
      } else {
      if (_saveobj->IsBlockPresent(_prevDataBlockNameSave))
      {
      Block& block = _saveobj->GetBlock(_prevDataBlockNameSave);
      block.WriteTable(_savetbl);
      }
      }
      _nTablesInBlockSave++;
      //      delete _savetbl;
      if (_prevtbl) delete _prevtbl;
      _prevtbl=new ISTable();
      *_prevtbl = *_savetbl;
    }
    if (_savetbl != NULL)
    {
        delete _savetbl;
        _savetbl = NULL;
    }
    _savetbl = new ISTable(categoryName);
    _curRowSave = 0;
    _curCategoryNameSave = categoryName;
    ProcessItemNameListSave();
  }
}


void DICParser::ProcessItemNameListSave(void)
/* ----------------------------------------------------------------------
   Purpose: DICParser::ProcessItemNameListSave(void)

            Registers the item keyword for the the current item in the
            current category.  Maintains an index array of "valid" keyword
            names in fieldList[].  This array is used to indirectly
            reference between keywords and values ...
 * ----------------------------------------------------------------------*/
{

  string keywordName;
  string categoryName;

  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Error in category name at line (save frame) " << NDBlineNo << " value " << _tBufKeyword << endl;
    errorLog += "Error in category name at line (save frame) ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
  }

  if ( _curItemNoSave > _fieldListAllocSave - 1) {
    _fieldListAllocSave = _curItemNoSave + _fieldListAllocSave;
    _fieldListSave.reserve(_fieldListAllocSave);
  }

  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if ((categoryName == _curCategoryNameSave) && !keywordName.empty()) {
    _savetbl->AddColumn(keywordName);
    if (_curItemNoSave >= (int)_fieldListSave.size())
      _fieldListSave.push_back(keywordName);
    else
      _fieldListSave[_curItemNoSave] = keywordName;
  } else {
    log << "Syntax error at line (save frame) " << NDBlineNo << " at item " << _tBufKeyword << endl;
    errorLog += "Syntax error at line (save frame) ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _fieldListSave[_curItemNoSave].clear();
  }
  _curItemNoSave++;

}
void DICParser::ProcessValueListSave(void)
/* ----------------------------------------------------------------------
     Purpose:  DICParser::ProcessValueListSave(void)

               Add the current value to the appropriate column in the
               the current row.  Start a new row if necessary.
 * ----------------------------------------------------------------------*/
{

  if (!_fieldListSave[_curValueNoSave].empty()) {

    if (_curValueNoSave == 0) {
      vector<string> rowBufSave(_savetbl->GetNumColumns(), CifString::UnknownValue);
      _curRowSave++;
      _savetbl->AddRow(rowBufSave);
    }
    if (!_pBufValue.empty()) {
         _savetbl->UpdateCell(_curRowSave-1, _fieldListSave[_curValueNoSave],
           _pBufValue);
      } else {
         _savetbl->UpdateCell(_curRowSave-1, _fieldListSave[_curValueNoSave],
           CifString::UnknownValue);
      }
  }
  _curValueNoSave++;

  if (_curValueNoSave == _curItemNoSave) {
    _curValueNoSave = 0;
  }
}





void DICParser::ProcessItemValuePairSave(void)
/* ----------------------------------------------------------------------
      Purpose: ProcessItemValuePairSave(void)

               Assign the current value to its associated item name.
 * ----------------------------------------------------------------------*/
{
  string categoryName;
  string keywordName;

  _curItemNoSave  = 1;
  _curValueNoSave = 0;

  try
  {
      CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  }
  catch (const exception& exc)
  {

  }

  if (categoryName.empty())
  {
    log << "Error in category name at line (save frame) " << NDBlineNo << " value " << _tBufKeyword << endl;
    errorLog += "Error in category name at line (save frame) ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    return;
  }

  bool tablePresent = false;

  if (_saveobj->IsBlockPresent(_curDataBlockNameSave))
  {
      Block& block = _saveobj->GetBlock(_curDataBlockNameSave);
      tablePresent = block.IsTablePresent(categoryName);
  }
  
  if (tablePresent) { //  duplicates a persistent table?
    log << "Duplicate category name in a save frame " << categoryName << " at line " << NDBlineNo << endl;
    errorLog += "Duplicate category name in a save frame ";
    errorLog += categoryName;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    return;
  }

  if (categoryName != _curCategoryNameSave) {

    if (_savetbl) { // write the current table
      if (_nTablesInBlockSave) {
      if (_saveobj->IsBlockPresent(_curDataBlockNameSave))
      {
          Block& block = _saveobj->GetBlock(_curDataBlockNameSave);
          block.WriteTable(_savetbl);
      }
      } else {
      if (_saveobj->IsBlockPresent(_prevDataBlockNameSave))
      {
          Block& block = _saveobj->GetBlock(_prevDataBlockNameSave);
          block.WriteTable(_savetbl);
      }
      }
      _nTablesInBlockSave++;
      //      delete _savetbl;
      if (_prevtbl) delete _prevtbl;
      _prevtbl=new ISTable();
      *_prevtbl = *_savetbl;
    }
    if (_savetbl != NULL)
    {
        delete _savetbl;
        _savetbl = NULL;
    }
    _savetbl = new ISTable(categoryName);
    _curCategoryNameSave = categoryName;
  }
  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if (!keywordName.empty()) {
//    cout<<categoryName<<"   "<<_curCategoryNameSave<<"   "<<endl;
    if (_savetbl==NULL)
      _savetbl = new ISTable(categoryName);
    _savetbl->AddColumn(keywordName);
  }
  else {
    log << "Syntax error line (save frame) " << NDBlineNo << " at item " << _tBufKeyword << endl;
    errorLog += "Syntax error at line (save frame) ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
  }
#ifdef VLAD_DELETED
  try
#endif
  {
  if (!_pBufValue.empty()) {
    if (_savetbl->GetNumRows() == 0)
        _savetbl->AddRow();
 
    _savetbl->UpdateCell(_savetbl->GetNumRows() - 1, keywordName, _pBufValue);

    if (_tBufKeyword == "_category.id")
        if (_pBufValue != _curDataBlockNameSave)
            log << "ERROR - In save frame \"save_" << _curDataBlockNameSave << 
              "\", \"_catgory.id\" has value \"" << _pBufValue <<
              "\" at line " << NDBlineNo << endl;
    } else
    {
    if (_savetbl->GetNumRows() == 0)
        _savetbl->AddRow();

    _savetbl->UpdateCell(_savetbl->GetNumRows() - 1, keywordName, CifString::UnknownValue);
    }
  }
#ifdef VLAD_DELETED
  catch (out_of_range)
#endif
  {
    // VLAD - BUG - Ignore for now but fix this. Reproduce in dict util with
    // make odb
  }

}

void  DICParser::CheckDDL(void)
{
  vector<string> target;
  target.push_back(_curCategoryName);

  unsigned int listOut = cattbl->FindFirst(target, listcat);
  if (listOut == cattbl->GetNumRows())
  {
    log << "Category " << _curCategoryName << " isn't defined; line " <<
      NDBlineNo<< endl;
    errorLog += "Category ";
    errorLog += _curCategoryName;
    errorLog += " isn't defined; line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    return;
  }

  unsigned int numRows2 = _tbl->GetNumRows();
  unsigned int numCols2 = _tbl->GetNumColumns();

  const vector<string>& ColumnNames = _prevtbl->GetColumnNames();
  unsigned int numCols = _prevtbl->GetNumColumns();
  unsigned int numRows = _prevtbl->GetNumRows();

  vector<string> listitem1;
  listitem1.push_back("name");

  bool implicitFound = false;
  bool pdbxImplicitFound = false;

  for (unsigned int i=0; i<numCols; i++)
  {
    string elem1;
    CifString::MakeCifItem(elem1, _prevtbl->GetName(), ColumnNames[i]);

    vector<string> target2;
    target2.push_back(elem1);
    unsigned int ret = itemtbl->FindFirst(target2, listitem1);
    if (ret == itemtbl->GetNumRows())
    {
      log << "Item " << elem1 << " isn't defined; line " << NDBlineNo<< endl;
      errorLog += "Item ";
      errorLog += elem1;
      errorLog += " isn't defined; line ";
      errorLog += String::IntToString(NDBlineNo);
      errorLog += '\n';
    }
    else
    {
        const string& mandCode = (*itemtbl)(ret, "mandatory_code");
        if (String::IsCiEqual(mandCode, "implicit"))
        {
#ifdef VLAD_EXAMINE
            cout << "IMPLICIT FOUND FOR itemName ==" << elem1 << endl;
#endif
            implicitFound = true;
        }
    }

    if (pdbxitemtbl != NULL)
    {
        ret = pdbxitemtbl->FindFirst(target2, listitem1);
        if (ret == pdbxitemtbl->GetNumRows())
        {
          log << "Item " << elem1 << " isn't defined; line " << NDBlineNo<< endl;
          errorLog += "Item ";
          errorLog += elem1;
          errorLog += " isn't defined; line ";
          errorLog += String::IntToString(NDBlineNo);
          errorLog += '\n';
        }
        else
        {
            const string& mandCode = (*pdbxitemtbl)(ret, "mandatory_code");
            if (String::IsCiEqual(mandCode, "implicit"))
            {
#ifdef VLAD_EXAMINE
                cout << "PDBX IMPLICIT FOUND FOR itemName ==" << elem1 << endl;
#endif
                pdbxImplicitFound = true;
            }
        }
    }
  }

  if (String::IsCiEqual(_curCategoryName,"item"))
  {
    if (!(_prevtbl->IsColumnPresent("name")))
    {
      if (numRows == numRows2)
      {
        _tbl->AddColumn("name");
      }

      for (int i=numRows; i>0; i--)
          _tbl->UpdateCell(numRows2-i, "name", _curDataBlockNameSave);
    }

    if (!_prevtbl->IsColumnPresent("category_id"))
    {
      string categoryName;
      if (isSave==2)
        CifString::GetCategoryFromCifItem(categoryName,_curDataBlockNameSave);
      if (isSave==1)
        CifString::GetCategoryFromCifItem(categoryName,_prevDataBlockNameSave);

      if (numRows == numRows2)
      {
        _tbl->AddColumn("category_id");
      }

      for (int i=numRows; i>0; i--)
          _tbl->UpdateCell(numRows2-i, "category_id", categoryName);
    }
  }
  else
  {
    if ((!implicitFound) && (!String::IsCiEqual(_curCategoryName,"pdbx_item")))
    //if (!implicitFound)
    {
      vector<string> target2;
      target2.push_back(_curCategoryName);
      target2.push_back("implicit");
      unsigned int listOut = itemtbl->FindFirst(target2, listitem2);
      if (listOut != itemtbl->GetNumRows())
      {
        const string& itemName = (*itemtbl)(listOut, "name");
        string keywordName;
        CifString::GetItemFromCifItem(keywordName,itemName);

        string elem;
        if (String::IsCiEqual(_curCategoryName,"dictionary") ||
            String::IsCiEqual(_curCategoryName,"category"))
        {
          elem = _curDataBlockName;
        }
        else
        {
          elem = _tmpDataBlockNameSave;
        }

        // Original implementation (the same as in the else part of ifdef)
        if (numCols == numCols2)
        {
          if (!_tbl->IsColumnPresent(keywordName))
              _tbl->AddColumn(keywordName);
        }
        
        for (int i=numRows; i>0; i--)
        {
            _tbl->UpdateCell(numRows2-i, keywordName, elem);
        }
      }
    }
  }


  if (String::IsCiEqual(_curCategoryName,"pdbx_item"))
  {
    if (!(_prevtbl->IsColumnPresent("name")))
    {
      if (numRows == numRows2)
      {
        _tbl->AddColumn("name");
      }

      for (int i=numRows; i>0; i--)
          _tbl->UpdateCell(numRows2-i, "name", _curDataBlockNameSave);
    }

    if (!_prevtbl->IsColumnPresent("category_id"))
    {
      string categoryName;
      if (isSave==2)
        CifString::GetCategoryFromCifItem(categoryName,_curDataBlockNameSave);
      if (isSave==1)
        CifString::GetCategoryFromCifItem(categoryName,_prevDataBlockNameSave);

      if (numRows == numRows2)
      {
        _tbl->AddColumn("category_id");
      }

      for (int i=numRows; i>0; i--)
          _tbl->UpdateCell(numRows2-i, "category_id", categoryName);
    }
  }
  else
  {
    if ((!pdbxImplicitFound) && (!String::IsCiEqual(_curCategoryName,"item")) && (pdbxitemtbl != NULL))
    //if (!pdbxImplicitFound)
    {
      vector<string> target2;
      target2.push_back(_curCategoryName);
      target2.push_back("implicit");
      unsigned int listOut = pdbxitemtbl->FindFirst(target2, listitem2);
      if (listOut != pdbxitemtbl->GetNumRows())
      {
        const string& itemName = (*pdbxitemtbl)(listOut, "name");
        string keywordName;
        CifString::GetItemFromCifItem(keywordName,itemName);

        string elem;
        if (String::IsCiEqual(_curCategoryName,"dictionary") ||
            String::IsCiEqual(_curCategoryName,"category"))
        {
          elem = _curDataBlockName;
        }
        else
        {
          elem = _tmpDataBlockNameSave;
        }

        // Original implementation (the same as in the else part of ifdef)
        if (numCols == numCols2)
        {
          if (!_tbl->IsColumnPresent(keywordName))
              _tbl->AddColumn(keywordName);
        }
        
        for (int i=numRows; i>0; i--)
        {
            _tbl->UpdateCell(numRows2-i, keywordName, elem);
        }
      }
    }
  }
}

void DICParser::ProcessAssignments(void)
{
  if (_curValueNo!=0)
  log <<"Number of data values is not exact multiples of the number of data names   (look above line)"<<NDBlineNo<<endl;
}

void DICParser::ProcessOneAssignment(void)
{
  if (isSave)
    ProcessItemValuePairSave();
  ProcessItemValuePair();
}

void DICParser::ProcessItemNameListLoop(void)
{
 if (isSave)
    ProcessLoopDeclarationSave();
  ProcessLoopDeclaration();
}

void DICParser::ProcessItemNameListName(void)
{
  if (isSave)
    ProcessItemNameListSave();
  ProcessItemNameList();
}

void DICParser::ProcessValueListItem(void)
{
  if (isSave)
    ProcessValueListSave();
  ProcessValueList();
}

void DICParser::ProcessItemName(void)
{
    _tBufKeyword = Glob_tBufKeywordSaveDIC;
}

void DICParser::ProcessLoop(void)
{
    _curItemNo = 0;  _curValueNo = 0;
    _curItemNoSave = 0;  _curValueNoSave = 0;
}

void DICParser::ProcessItemValue(void)
{
   _pBufValue = Glob_pBufValueDIC;
}

void DICParser::ProcessLsItemValue(void)
{
   _pBufValue = Glob_pBufValueDIC;
}

void DICParser::ProcessUnknownValue(void)
{
   _pBufValue = CifString::InapplicableValue;
}

void DICParser::ProcessMissingValue(void)
{
   _pBufValue = CifString::UnknownValue;
}

void DICParser::ProcessSaveBegin(void)
{
  _nTablesInBlockSave=1;
  _saveobj = new CifFile();
  _prevDataBlockNameSave = _curDataBlockNameSave;
  _curDataBlockNameSave = &(Glob_dataBlockNameDIC)[5];

  // Check if the save frame was already indicated
  if (_saveFrames.find(_curDataBlockNameSave) != _saveFrames.end())
  {
      log << "INFO - Duplicate save frame \"" << _curDataBlockNameSave <<
        "\" at line " << NDBlineNo << endl;
  }
  else
  {
      _saveFrames.insert(_curDataBlockNameSave);
  }

  _savetbl=NULL;
  _pBufValue = Glob_dataBlockNameDIC;
}

void DICParser::ProcessSaveEnd(void)
{
  if (_prevtbl) delete _prevtbl;
  _prevtbl=new ISTable();
  *_prevtbl = *_savetbl;
  delete _saveobj;
  if (_savetbl != NULL)
  {
      delete _savetbl;
      _savetbl = NULL;
  }
  _pBufValue = Glob_pBufValueDIC; 
}

void DICParser::ProcessDataBlockName(void)
{
  _numDataBlocks++;
  _nTablesInBlock=0;

  if (&(Glob_dataBlockNameDIC)[5] && (strlen(&(Glob_dataBlockNameDIC)[5])>0)) {
    _curDataBlockName = &(Glob_dataBlockNameDIC)[5];
    _fobj->AddBlock(_curDataBlockName);
  } else {
    _curDataBlockName = "UNNAMED-";
    _curDataBlockName += String::IntToString(_numDataBlocks);
  }
//  if (_prevDataBlockName.size() == 0)   _prevDataBlockName = _curDataBlockName;
  if (_prevDataBlockName == string("MISSING_DIC"))   _prevDataBlockName = _curDataBlockName;

  if (_tbl && (_curItemNo > 0)) { // write the current table / management of _tbl by _fobj
  CheckDDL();
      Block& block = _fobj->GetBlock(_prevDataBlockName);
      block.WriteTable(_tbl);
      _nTablesInBlock++;
    //    delete _tbl;
                _tbl=NULL;
  }
  if (_prevDataBlockName != _curDataBlockName)
         _curCategoryName.clear();

  if (_tmpDataBlockNameSave != _curDataBlockNameSave) {
    _tmpDataBlockNameSave.clear();
  }

#if DEBUG
  if (_verbose) log << " Previous data block is now  " <<  _prevDataBlockName << endl;
  if (_verbose) log << " Current  data block is now  " <<  _curDataBlockName << endl;
#endif
}



void ProcessAssignmentsFromDICParser()
{
    DICParserP->ProcessAssignments();
}

void ProcessOneAssignmentFromDICParser()
{
    DICParserP->ProcessOneAssignment();
}

void ProcessItemNameListLoopFromDICParser()
{
    DICParserP->ProcessItemNameListLoop();
}

void ProcessItemNameListNameFromDICParser()
{
    DICParserP->ProcessItemNameListName();
}

void ProcessValueListFromDICParser()
{
    DICParserP->ProcessValueListItem();
}

void ProcessItemNameFromDICParser()
{
    DICParserP->ProcessItemName();
}

void ProcessLoopFromDICParser()
{
    DICParserP->ProcessLoop();
}

void ProcessItemValueFromDICParser()
{
    DICParserP->ProcessItemValue();
}

void ProcessLsItemValueFromDICParser()
{
    DICParserP->ProcessLsItemValue();
}

void ProcessUnknownValueFromDICParser()
{
    DICParserP->ProcessUnknownValue();
}

void ProcessMissingValueFromDICParser()
{
    DICParserP->ProcessMissingValue();
}

void ProcessSaveBeginFromDICParser()
{
    DICParserP->ProcessSaveBegin();
}

void ProcessSaveEndFromDICParser()
{
    DICParserP->ProcessSaveEnd();
}

void ProcessDataBlockNameFromDICParser()
{
    DICParserP->ProcessDataBlockName();
}


void DICParser::AfterParseProcessing()
{
    InsertImplicitOrdinalItems();
}


void DICParser::InsertImplicitOrdinalItems()
{
    Block& ddlBlock = ddl->GetBlock(ddl->GetFirstBlockName());
    ISTable* itemTableP = ddlBlock.GetTablePtr("item");

    vector<string> searchCols;
    searchCols.push_back("mandatory_code");

    vector<string> searchVals;
    searchVals.push_back("implicit-ordinal");

    vector<unsigned int> foundIndices;
    itemTableP->Search(foundIndices, searchVals, searchCols);

    Block& block = _fobj->GetBlock(_fobj->GetFirstBlockName());

    for (unsigned int foundI = 0; foundI < foundIndices.size(); ++foundI)
    {
        const string& catName = (*itemTableP)(foundIndices[foundI],
          "category_id");

        ISTable* catTableP = block.GetTablePtr(catName);
        if (catTableP == NULL)
        {
            continue;
        }

        const string& itemName = (*itemTableP)(foundIndices[foundI], "name");

        string attrName;
        CifString::GetItemFromCifItem(attrName, itemName); 

        if (!catTableP->IsColumnPresent(attrName))
        {
            catTableP->AddColumn(attrName);
        }

        for (unsigned int rowI = 0; rowI < catTableP->GetNumRows(); ++rowI)
        {
            catTableP->UpdateCell(rowI, attrName,
              String::IntToString(rowI + 1));
        }
    }
} // End of DICParser::InsertImplicitOrdinalItems()

