//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file DICParserBase.h
**
** \brief Header file for DIC Parser class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/


#ifndef DIC_PARSER_BASE_H
#define DIC_PARSER_BASE_H


#include <iostream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DICScannerBase.h"
#include "DICParserInt.h"
#include "CifFileReadDef.h"
#include "ISTable.h"
#include "DicFile.h"

/**
**  \class DICParser
**
**  \brief Public class that respresents a dictionary parser.
**
**  This class represents a dictionary parser. This class utilizes flex/bison
**  for syntax/semantic processing and stores the parsed data (dictionary
**  blocks and tables) in a \e DicFile object.
*/
class DICParser : public DICScanner
{
    public:
        /**
        **  Constructs a dictionary parser.
        **
        **  \param[in] dicFileP - pointer to the \e DicFile object that the
        **    dictionary parser is to use to store the parsed data
        **  \param[in] ddlFileP - pointer to the \e CifFile object that holds
        **    the DDL for the dictionary
        **  \param[in] verbose - optional parameter that indicates whether
        **    parsing logging should be turned on (if true) or off (if false).
        **    If \e verbose is not specified, logging is turned off.
        **
        **  \return Not applicable
        ** 
        **  \pre \e dicFileP must not be NULL
        **  \pre \e ddlFileP must not be NULL
        ** 
        **  \post None
        **
        **  \exception EmptyValueException - if \e dicFileP is NULL
        **  \exception EmptyValueException - if \e ddlFileP is NULL
        */
        DICParser(DicFile* dicFileP, CifFile* ddlFileP, bool verbose = false);

        /**
        **  Parses a dictionary file.
        **
        **  \param[in] fileName - relative or absolute name of the dictionary
        **    file that is to be parsed.
        **  \param[out] diagnostics - parsing result. If empty, parsing
        **    completed with no warnings or errors. If non-empty, there were
        **    parsing warnings and/or parsing errors.
        **
        **  \return None
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        void Parse(const string& fileName, string& diagnostics);

        /**
        **  Destructs a dictionary parser by releasing all the used resources.
        **  
        **  \param: Not applicable
        **  
        **  \return Not applicable
        **  
        **  \pre None
        **  
        **  \post None
        **  
        **  \exception: None
        */
        virtual ~DICParser();

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void Error(const char*);

        /**
        **  Method, not currently part of users public API, and will soon be
        **  re-examined.
        */
        void Clear();

        /**
        **  Method, not currently part of users public API, and will soon be
        **  re-examined.
        */
        void Reset();

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessAssignments(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessOneAssignment(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessItemNameListLoop(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessItemNameListName(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessValueListItem(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessItemName(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessLoop(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessItemValue(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessLsItemValue(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessUnknownValue(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessMissingValue(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessSaveBegin(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessSaveEnd(void);

        /** 
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessDataBlockName(void);

    private:
        DicFile *_fobj;
        ISTable *_tbl;
        int _afterLoop;
        CifFile *_saveobj;
        ISTable *_savetbl;
        ISTable *_prevtbl;
        ISTable * format;
        ISTable * cattbl;
        ISTable * itemtbl;
        ISTable * pdbxitemtbl;
        CifFile *ddl;
        vector<string> listcat, listitem;
        vector<string> listitem2;
        int _nTablesInBlock;
        int _curItemNo, _curValueNo, _numDataBlocks, _fieldListAlloc, _curRow;
        vector<string> _fieldList;
        string _pBufValue;
        string _tBufKeyword;
        string _curCategoryName;
        string _curDataBlockName;
        string _prevDataBlockName;
        int _nTablesInBlockSave;
        int _curItemNoSave, _curValueNoSave;
        int _fieldListAllocSave;
        int _curRowSave;
        vector<string> _fieldListSave;
        string _curCategoryNameSave;
        string _curDataBlockNameSave;
        string _prevDataBlockNameSave;
        string _tmpDataBlockNameSave;
        string errorLog;
        std::set<string> _saveFrames;
        void ProcessLoopDeclaration(void);
        void ProcessItemNameList(void);
        void ProcessValueList(void);
        void ProcessItemValuePair(void);
        void ProcessLoopDeclarationSave(void);
        void ProcessItemNameListSave(void);
        void ProcessValueListSave(void);
        void ProcessItemValuePairSave(void);
        void CheckDDL(void);

        void AfterParseProcessing();

        void InsertImplicitOrdinalItems();
};
 
#endif /* DIC_PARSER_BASE_H */
