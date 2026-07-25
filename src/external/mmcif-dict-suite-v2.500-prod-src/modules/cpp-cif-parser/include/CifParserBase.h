//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifParserBase.h
**
** \brief Header file for CifParser class.
*/


#ifndef CIF_PARSER_BASE_H
#define CIF_PARSER_BASE_H


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ISTable.h"
#include "CifFile.h"
#include "CifScannerBase.h"
#include "CifParserInt.h"
#include "CifFileReadDef.h"


#define  DATA_TAG "data_"


/**
**  \class CifParser
**
**  \brief Public class that respresents a CIF parser.
**
**  This class represents a CIF parser. This class utilizes flex/bison for
**  syntax/semantic processing and stores the parsed data (data blocks and
**  tables) in a \e CifFile object.
*/
class CifParser : public CifScanner
{
    public:
        /**
        **  Constructs a CIF parser.
        **
        **  \param[in] cifFileP - pointer to the \e CifFile object that the
        **    CIF parser is to use to store the parsed data
        **  \param[in] verbose - optional parameter that indicates whether
        **    parsing logging should be turned on (if true) or off (if false).
        **    If \e verbose is not specified, logging is turned off.
        **
        **  \return Not applicable
        **
        **  \pre \e cifFileP must not be NULL
        **
        **  \post None
        **
        **  \exception EmptyValueException - if \e cifFileP is NULL
        */
        CifParser(CifFile* cifFileP, bool verbose = false);

        /**
        **  Method, not currently part of users public API, and will soon be
        **  re-examined.
        */
        CifParser(CifFile* cifFileP, CifFileReadDef readDef,
          bool verbose = false);

        /**
        **  Parses the CIF file.
        **
        **  \param[in] fileName - relative or absolute name of the CIF file
        **    that is to be parsed.
        **  \param[in] parseLogFileName - relative or absolute name of the file
        **    where parsing log is to be stored.
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
        **  \exception NotFoundException - if file with name \e fileName
        **    does not exist
        */
        void Parse(const string& fileName, string& diagnostics,
          const std::string& parseLogFileName = std::string());

        /**
        **  Parses the CIF data in a string.
        **
        **  \param[in] cifString - a string that contains CIF data that is to
        **    be parsed.
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
        void ParseString(const string& cifString, string& diagnostics);

        /**
        **  Destructs a CIF parser by releasing all the used resources.
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
        virtual ~CifParser();

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
        int ProcessLoopDeclaration(void);

        /**
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        int ProcessItemNameList(void);

        /**
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        int ProcessValueList(void);

        /**
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        int ProcessItemValuePair(void);

        /**
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessAssignments(void);

        /**
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessLoop(void);

        /**
        **  Utility method, not part of users public API, and will soon be
        **  removed.
        */
        void ProcessItemName(void);

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
        void ProcessDataBlockName(void);

    private:
        CifFile *_fobj;
        CifFileReadDef _readDef;
        ISTable *_tbl;
        int _afterLoop;
        int _nTablesInBlock;
        int _curItemNo, _curValueNo, _numDataBlocks, _fieldListAlloc, _curRow;
        vector<string> _fieldList;
        string _pBufValue;
        string _tBufKeyword;
        string _curCategoryName;
        string _curDataBlockName;
        string _prevDataBlockName;
        void _ComplexWriteTable();

        int _err, _warn;
};
 
#endif /* CIF_PARSER_BASE_H */
