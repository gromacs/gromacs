//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifFile.h
**
** \brief Header file for CifFile class.
*/


/*
  PURPOSE:    Base class for read/write cif files
*/


#ifndef CIFFILE_H
#define CIFFILE_H


#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "GenString.h"
#include "CifString.h"
#include "TableFile.h"
#include "CifParentChild.h"



/**
**  \class CifFile
**
**  \brief Public class that represents a CIF file, composed of blocks with
**    tables.
**
**  This class represents a CIF file. In addition to inherited methods from
**  \e TableFile class, this class provides methods for writing the data to
**  a text file, along with methods for controlling how the data is written,
**  and a method for verifying the CIF file against dictionary.
*/
class CifFile : public TableFile
{
  public:
    std::string _parsingDiags;
    std::string _checkingDiags;

    static const unsigned int STD_CIF_LINE_LENGTH = 80;

    enum eQuoting
    {
        eSINGLE = 0,
        eDOUBLE
    };

    /**
    **  Constructs a CIF file.
    **
    **  \param[in] fileMode - CIF file mode. Possible values are
    **    read-only, create, update and virtual. Detailed description of
    **    file mode is given in \e TableFile documentation.
    **  \param[in] fileName - relative or absolute name of the file
    **    where object persistency is maintained. If \e fileMode specifies
    **    virtual mode, this parameter is ignored.
    **  \param[in] verbose - optional parameter that indicates whether
    **    logging should be turned on (if true) or off (if false).
    **    If \e verbose is not specified, logging is turned off.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **  \param[in] maxLineLength - optional parameter that indicates the
    **    maximum number of written characters in one line in the written
    **    text file. If not specified, \e STD_CIF_LINE_LENGTH is used.
    **  \param[in] nullValue - optional parameter that indicates the
    **    character that is to be used to denote unknown value in the written
    **    CIF file. If not specified, \e CifString::UnknownValue is used.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    CifFile(const eFileMode fileMode, const std::string& fileName,
      const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Constructs a CIF file in virtual mode.
    **
    **  \param[in] verbose - optional parameter that indicates whether
    **    logging should be turned on (if true) or off (if false).
    **    If \e verbose is not specified, logging is turned off.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **  \param[in] maxLineLength - optional parameter that indicates the
    **    maximum number of written characters in one line in the written
    **    text file. If not specified, \e STD_CIF_LINE_LENGTH is used.
    **  \param[in] nullValue - optional parameter that indicates the
    **    character that is to be used to denote unknown value in the written
    **    CIF file. If not specified, \e CifString::UnknownValue is used.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    CifFile(const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /* Fake method to export the above method to other languages. */
    CifFile(const bool fake, const bool verbose,
      const unsigned int intCaseSense = 0,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Destructs a CIF file, by releasing all consumed resources.
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
    ~CifFile();

    /**
    **  Sets file name of a file that was the source of the object data.
    **
    **  \param srcFileName - The name of the source data file.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetSrcFileName(const std::string& srcFileName);


    /**
    **  Retrieves source file name.
    **
    **  \param: None
    **
    **  \return - source file name
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    const std::string& GetSrcFileName();


    /**
    **  Retrieves logging option.
    **
    **  \param: None
    **
    **  \return true - if logging is turned on
    **  \return false - if logging is turned off
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline bool GetVerbose();

    /**
    **  Sets smart printing option. Smart printing is used to beautify the
    **  output of a written text file.
    **
    **  \param smartPrint - smart printing. If false, smart printing is
    **    disabled. If true, smart printing is enabled. If not specified,
    **    smart printing is enabled.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline void SetSmartPrint(bool smartPrint = true);


    /**
    **  Retrieves smart printing option.
    **
    **  \param: None
    **
    **  \return true - if smart printing is enabled
    **  \return false - if smart printing is disabled
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline bool IsSmartPrint();

    /**
    **  Sets quoting option. This option is used in order to
    **  select the type of quoting to be used in the written text file.
    **
    **  \param quoting - type of quoting. If \e eSINGLE, single quotes are
    **    used. If \e eDOUBLE, double quotes are used.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetQuoting(eQuoting quoting);

    /**
    **  Retrieves quoting option.
    **
    **  \param: None
    **
    **  \return eSINGLE - if single quotes are used
    **  \return eDOUBLE - if double quotes are used
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    unsigned int GetQuoting();

    /**
    **  This method is used in order to control how single row categories are
    **  written: in form of a "loop_" construct or as an item-value pair.
    **
    **  \param catName - category name
    **  \param looping - category looping option. If false and the
    **    category is a single row category, that category will not be
    **    written with "loop_" construct. Otherwise, if true, single row
    **    category will be written with "loop_" construct.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetLooping(const std::string& catName, bool looping = false);

    /**
    **  Retrieves looping option of a category.
    **
    **  \param catName - category name
    **
    **  \return - category looping option, as described in SetLooping() method.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool GetLooping(const std::string& catName);

    /**
    **  Writes the data out to a text file.
    **
    **  \param[in] cifFileName - relative or absolute name of the text file
    **    to which the data from \e CifFile object is to be written to.
    **  \param[in] sortTables - optional parameter that indicates whether
    **    written tables should be sorted (if true) or not sorted (if false).
    **    If \e sortTables is not specified, tables are not sorted prior to
    **    writing them.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Write(const std::string& cifFileName, const bool sortTables = false,
      const bool writeEmptyTables = false);

    /**
    **  Writes the data out to a text file.
    **
    **  \param[in] cifFileName - relative or absolute name of the text file
    **    to which the data from \e CifFile object is to be written to.
    **  \param[in] tableOrder - vector of table names that indicates the
    **    order of written tables.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Write(const std::string& cifFileName,
      const std::vector<std::string>& tableOrder,
      const bool writeEmptyTables = false);

    /**
    **  Writes the data out to an output stream.
    **
    **  \param[in] outStream - a reference to the output stream
    **  \param[in] sortTables - optional parameter that indicates whether
    **    written tables should be sorted (if true) or not sorted (if false).
    **    If \e sortTables is not specified, tables are not sorted prior to
    **    writing them.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Write(std::ostream& outStream, const bool sortTables = false,
      const bool writeEmptyTables = false);

    /**
    **  Writes the data out to a text file in NMR-STAR format.
    **
    **  \param[in] nmrStarFileName - relative or absolute name of the text file
    **    to which the data from \e CifFile object is to be written to.
    **  \param[in] globalBlockName - the name of the global NMR-STAR block.
    **  \param[in] sortTables - optional parameter that indicates whether
    **    written tables should be sorted (if true) or not sorted (if false).
    **    If \e sortTables is not specified, tables are not sorted prior to
    **    writing them.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteNmrStar(const std::string& nmrStarFileName,
      const std::string& globalBlockName,  const bool sortTables = false,
      const bool writeEmptyTables = false);

    /**
    **  Checks a CIF file (all blocks in it) against the dictionary.
    **
    **  \param[in] dicRef - reference to a dictionary file. The check is
    **    done against the first block in the dictionary file.
    **  \param[in] diagFileName - relative or absolute name of the file,
    **    where diagnostic information is stored.
    **  \param[in] extraDictChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, dictionary checks. If not
    **    specified, those checks are not performed.
    **  \param[in] extraCifChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, CIF checks. If not specified,
    **    those checks are not performed.
    **
    **  \return 0 - if all checks passed
    **  \return different than 0 - if checks failed
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    int DataChecking(CifFile& dicRef, const std::string& diagFileName,
      const bool extraDictChecks = false, const bool extraCifChecks = false);

    /**
    **  Checks a block of CIF file against the specified reference block.
    **
    **  \param[in] block - reference to a block that is to be checked
    **  \param[in] refBlock - reference to a reference block against which
    **    \e block is to be checked
    **  \param[out] diagBuf - diagnostics buffer that holds checking results
    **  \param[in] extraDictChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, checks. If not specified, those
    **    checks are not performed.
    **  \param[in] extraCifChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, CIF checks. If not specified,
    **    those checks are not performed.
    **
    **  \return 0 - if all checks passed
    **  \return different than 0 - if checks failed
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    int DataChecking(Block& block, Block& refBlock, std::ostringstream& buf,
      const bool extraDictChecks = false, const bool extraCifChecks = false);

    /**
    **  Sets enumerations checking option for case-insensitive types.
    **
    **  \param caseSense - case sensitivity of enumeration values check. If
    **    false, enumeration values of case-insensitive types will be checked
    **    as case-insensitive. If true, enumeration values of case-insensitive
    **    types will be checked as case-sensitive.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetEnumCheck(bool caseSense = false);

    /**
    **  Retrieves enumerations checking option for case-insensitive types.
    **
    **  \param: None
    **
    **  \return true - if case-sensitive enumeration check is enabled
    **  \return false - if case-insensitive enumeration check is enabled
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool GetEnumCheck();

    /**
    **  Gets parsing diagnostics.
    **
    **  \param: None
    **
    **  \return - reference to parsing diagnostics
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    const std::string& GetParsingDiags();


    /**
    **  Finds indices of rows that contain all CIF null values. A CIF null
    **  value is defined as a "?" or "".
    **
    **  \param[out] nullRowsIndices - vector of null rows indices.
    **  \param[in] isTable - table reference
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void FindCifNullRows(std::vector<unsigned int>& nullRowsIndices,
      const ISTable& isTable);

    void GetAttributeValue(std::string& attribVal, const std::string& blockId,
      const std::string& category, const std::string& attribute);
    void GetAttributeValueIf(std::string& attribVal, const std::string& blockId,
      const std::string& category, const std::string& attributeA,
      const std::string& attributeB, const std::string& valB);
    bool IsAttributeValueDefined(const std::string& blockId,
      const std::string& category, const std::string& attribute);

    void SetAttributeValue(const std::string& blockId,
      const std::string& category,
      const std::string& attribute, const std::string& value,
      const bool create = false);
    void SetAttributeValueIf(const std::string& blockId,
      const std::string& category, const std::string& attributeA,
      const std::string& valA,
      const std::string& attributeB, const std::string& valB,
       const bool create = false);
    void SetAttributeValueIfNull(const std::string& blockId,
      const std::string& category, const std::string& attribute,
      const std::string& value);

    void GetAttributeValues(std::vector<std::string>& strings,
      const std::string& blockId,
      const std::string& category, const std::string& attribute);
    void GetAttributeValuesIf(std::vector<std::string>& strings,
      const std::string& blockId, const std::string& category,
      const std::string& attributeA,
      const std::string& attributeB, const std::string& valB);

    void SetAttributeValues(const std::string& blockId,
      const std::string& category, const std::string& attribute,
      const std::vector<std::string>& values);

#ifdef VLAD_TO_CIF_FILE_NOT_USED
    void del_attribute_value_where(CifFile *fobj, const char *blockId,
      const char *category, const char *attributeB, const char *valB);
#endif // VLAD_TO_CIF_FILE_NOT_USED not defined

    int CheckCategories(Block& block, Block& refBlock, std::ostringstream& log);
    void CheckCategoryKey(Block& block, std::ostringstream& log);
    void CheckItemsTable(Block& block, std::ostringstream& log);
    int CheckItems(Block& block, Block& refBlock, std::ostringstream& log);


  protected:
    static const unsigned int STD_PRINT_SPACING = 3;
    static const unsigned int SMART_PRINT_SPACING = 1;
    static const unsigned int HEADER_SPACING = 40;

    enum eIdentType
    {
        eNONE = 0,
        eLEFT,
        eRIGHT
    };

    std::string _beginDataKeyword;
    std::string _endDataKeyword;

    std::string _beginLoopKeyword;
    std::string _endLoopKeyword;

    unsigned int _maxCifLineLength;
    std::string _nullValue;
    bool _verbose;
    bool _smartPrint;
    std::string _quotes;
    std::map<std::string, bool> _looping;
    bool _enumCaseSense;

    int _IsQuotableText(const std::string& itemValue);
    eIdentType _FindPrintType(const std::vector<std::string>& values);

    void _PrintItemIdent(std::ostream& cifo, unsigned int& linePos);
    void _PrintItemName(std::ostream& cifo, const std::string& category,
      const std::string& itemName, unsigned int& linePos);
    void _PrintPostItemSeparator(std::ostream& cifo, unsigned int& linePos,
      const bool ident = false, const unsigned int numSpaces = 1);

    int _PrintItemValue(std::ostream& cifo, const std::string& itemValue,
      unsigned int& linePos, const eIdentType identType = eNONE,
      const unsigned int width = 0);

    int _PrintItemNameInHeader(std::ostream& cifo, const std::string& itemValue,
      unsigned int& linePos, const eIdentType identType = eNONE,
      const unsigned int width = 0);

    void _PrintHeaderedItems(std::ostream& cifo,
      const std::vector<std::string>& colNames,
      const std::vector<unsigned int>& colWidths,
      const std::vector<eIdentType> colPrintType);

    void Write(std::ostream& cifo, const std::vector<std::string>& catOrder,
      const bool writeEmptyTables = false);

    void Write(std::ostream& cifo, std::vector<unsigned int>& tables,
      const bool writeEmptyTables = false);


  private:
    std::string _srcFileName;

    bool _extraDictChecks;
    bool _extraCifChecks;

    void Init();

    bool IsCatDefinedInRef(const std::string& catName, ISTable& catTable);
    bool IsItemDefinedInRef(const std::string& catName,
      const std::string& itemName, ISTable& refItemTable);
    bool IsImplicitNatureKey(const string& catName, const string& attribName,
      ISTable& itemTable);
    void GetImplNatureKeysWithMissingValues(vector<string>& implKeyItems,
      const string& catName, ISTable& catTable, ISTable& refItemTable);
    void GetImplNatureKeys(vector<string>& implNatureKeys,
      const string& catName, ISTable& refItemTable);
    bool AreSomeValuesInColumnEmpty(ISTable& table, const string& colName);
    void FixMissingValuesOfImplNatureKeys(ISTable& catTable,
      const vector<string>& implNatureKeys, ISTable& refItemDefaultTable,
      std::ostringstream& log);
    void GetItemDefaultValue(string& defValue, const string& implNatKey,
      ISTable& refItemDefaultTable);
    void CheckKeyItems(const std::string& blockName, ISTable& catTable,
      ISTable& keyTable, std::ostringstream& log);
    void CheckKeyValues(const std::vector<std::string>& keyItems,
      ISTable& catTable, std::ostringstream& log);

    void GetKeyAttributes(std::vector<std::string>& keyAttributes,
      const std::string& catTableName, ISTable& catKeyTable);
    void CheckKeyItems(const std::string& blockName, ISTable& catTable,
      const std::vector<std::string>& keyAttributes, ISTable& itemTable,
      ISTable* itemDefaultTableP, std::ostringstream& log);

    void CheckMandatoryItems(const std::string& blockName, ISTable& catTable,
      ISTable& refItemTable, const std::vector<std::string>& keyItems,
      std::ostringstream& log);

    void CheckAndRectifyItemTypeCode(Block& block, std::ostringstream& log);
    void RectifyItemTypeCode(std::string& retItemTypeCode,
      std::ostringstream& log, Block& block, CifParentChild& cifParentChild,
      const std::string& cifItemName);

    int CheckRegExpRangeEnum(Block& block, ISTable& catTable,
      const std::string& attribName, ISTable& itemTypeTable,
      ISTable& itemTypeListTable, ISTable* itemRangeTableP,
      ISTable* itemEnumTableP, ISTable& parChildTable, ISTable* itemAliasesP,
      std::ostringstream& log);

    int CheckCellRange(const std::string& cell, const std::string& typeCode,
      const std::vector<std::string>& minlist,
      const std::vector<std::string>& maxlist);

    int CheckCellEnum(const std::string& cell, const std::string& typeCode,
      const std::string& primCode, const std::vector<std::string>& enumlist);

    int CheckCellFloatRange(const std::string& cell,
      const std::vector<std::string>& minlist,
      const std::vector<std::string>& maxlist);

    int CheckCellIntRange(const std::string& cell,
      const std::vector<std::string>& minlist,
      const std::vector<std::string>& maxlist);

    int CheckCellFloatEnum(const std::string& cell,
      const std::vector<std::string>& enumlist);

    int CheckCellIntEnum(const std::string& cell,
      const std::vector<std::string>& enumlist);

    int CheckCellOtherEnum(const std::string& cell, const std::string& primCode,
      const std::vector<std::string>& enumlist);

    void GetItemTypeCode(std::string& typeCode, const std::string& cifItemName,
      ISTable& itemTypeTable);

    void ConvertEscapedString(const std::string& inString,
      std::string& outString);
};


inline bool CifFile::GetVerbose()
{
    return(_verbose);
}


inline void CifFile::SetSmartPrint(bool smartPrint)
{
    _smartPrint = smartPrint;
}


inline bool CifFile::IsSmartPrint()
{
    return(_smartPrint);
}


#endif
