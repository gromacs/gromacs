//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file TableFile.h
**
** \brief Header file for Block and TableFile classes.
*/


#ifndef TABLEFILE_H
#define TABLEFILE_H


#include <vector>
#include <set>

#include "rcsb_types.h"
#include "mapped_ptr_vector.h"
#include "mapped_ptr_vector.C"
#include "GenString.h"
#include "ISTable.h"
#include "Serializer.h"


/**
**  \class Block
**
**  \brief Public class that represents a data block, that contains tables.
**
**  This class represents a data block, that can come from DDL,
**  dictionary or CIF files. Data block is a container of tables.
**  This class provides methods for construction and destruction, tables
**  manipulation (addition, retrieval, deleting, writing), data blocks
**  comparison.
*/
class Block
{
  public:
    mapped_ptr_vector<ISTable, StringLess> _tables;

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    **
    **  Constructs a data block.
    **
    **  \param[in] name - the name of the data block
    **  \param[in] serP - pointer to the File Navigator object
    **  \param[in] fileMode - optional parameter that indicates data block
    **    mode. Possible values are read-only, create, update and virtual.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    Block(const string& name, Serializer* serP,
      const eFileMode fileMode = READ_MODE, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    **
    **  Destructs a data block.
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
    ~Block();

    /**
    **  Compares a data block to another data block.
    **
    **  \param[in] inBlock - reference to input data block
    **
    **  \return vector of pairs, where the first value in a pair is a
    **    table name and the second value in a pair is one of the following
    **    indicators of table differences: \n \n
    **  eMISSING - table exists only in the input block and not in this
    **    block \n
    **  eEXTRA - table exists only in this block and not in the input block \n
    **  eCASE_SENSE - table exists in both blocks, but with different column
    **    name case sensitivity \n
    **  eMORE_COLS - table exists in both blocks, but the table in this block
    **    has more columns than the table in the input block \n
    **  eLESS_COLS - table exists in both blocks, but the table in this block
    **    has less columns than the table in the input block \n
    **  eCOL_NAMES - table exists in both blocks, but tables have different
    **    column names \n
    **  eMORE_ROWS - table exists in both blocks, but the table in this block
    **    has more rows than the table in the input block \n
    **  eLESS_ROWS - table exists in both blocks, but the table in this block
    **    has less rows than the table in the input block \n
    **  eCELLS - table exists in both blocks, but tables have different
    **    content \n
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    vector<pair<string, ISTable::eTableDiff> > operator==(Block& inBlock);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    **
    **  Sets the name of a data block.
    **
    **  \param[in] name - the name of the data block
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline void SetName(const string& name);

    /**
    **  Retrieves data block name.
    **
    **  \param: None
    **
    **  \return Constant reference to a string that contains data block name.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline const string& GetName() const;

    /**
    **  Adds a table to the block. If a table with the specified name
    **  already exists, it will be overwritten.
    **
    **  \param[in] name - optional parameter that indicates the name of the
    **    table to be added
    **  \param[in] colCaseSense - optional parameter that indicates case
    **    sensitivity of column names. Possible values are case sensitive and
    **    case in-sensitive. If not specified, a table with case sensitive
    **    column names is constructed.
    **
    **  \return Reference to the table
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    ISTable& AddTable(const std::string& name = string(),
      const Char::eCompareType colCaseSense = Char::eCASE_SENSITIVE);


    /**
    **  Changes the name of a table in the data block.
    **
    **  \param[in] oldName - the name of the table which is to be renamed
    **  \param[in] newName - the new table name
    **
    **  \return None
    **
    **  \pre \e oldName must be non-empty
    **  \pre Table with name \e oldName must be present
    **  \pre \e newName must be non-empty
    **  \pre Table with name \e newName must not be present
    **  \pre Block must be in create or update mode
    **
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e oldName is empty
    **  \exception NotFoundException - if table with name \e oldName does
    **    not exist
    **  \exception EmptyValueException - if \e newName is empty
    **  \exception AlreadyExistsException - if table with name \e newName
    **    already exists
    **  \exception FileModeException - if block is not in create or
    **    update mode
    */
    void RenameTable(const string& oldName, const string& newName);

    /**
    **  Retrieves names of all tables in a data block.
    **
    **  \param[out] tableNames - retrieved table names
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void GetTableNames(vector<string>& tableNames);

    /**
    **  Checks for table presence in the data block.
    **
    **  \param[in] tableName - table name
    **
    **  \return true - if table exists
    **  \return false - if table does not exist
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool IsTablePresent(const string& tableName);

    /**
    **  Retrieves a table reference.
    **
    **  \param[in] tableName - table name
    **
    **  \return Reference to the table, if table was found
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception NotFoundException - if table with name \e tableName
    **    does not exist
    */
    ISTable& GetTable(const string& tableName);

    /**
    **  Retrieves a pointer to the table.
    **
    **  \param[in] tableName - table name
    **
    **  \return Pointer to the table, if table was found
    **  \return NULL, if table was not found
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    ISTable* GetTablePtr(const string& tableName);

    /**
    **  Deletes a table from a data block.
    **
    **  \param[in] tableName - table name
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void DeleteTable(const string& tableName);

    /**
    **  Writes a table to the data block. In this context, writing means
    **  adding it (if it does not already exist) or updating it (if it
    **  already exists).
    **
    **  \param[in] isTable - reference to the table
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteTable(ISTable& isTable);

    /**
    **  Writes a table to the data block. In this context, writing means
    **  adding it (if it does not already exist) or updating it (if it
    **  already exists).
    **
    **  \param[in] isTableP - pointer to the table
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteTable(ISTable* isTableP);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    */
    void Print();

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    *
    *  JDW rename this method to resolve llvm ambiguity issues with public method signature -
    */
    void _AddTable(const string& name, const int indexInFile = 0,
      ISTable* isTableP = NULL);

  private:
    string _name;
    eFileMode _fileMode;
    Serializer* _ser;

    Block(const Block& t);
    Block& operator=(const Block& inBlock);

    ISTable* _GetTablePtr(const unsigned int tableIndex);


};


/**
**  \class TableFile
**
**  \brief Public class that represents a file composed of blocks with tables.
**
**  This class represents an ordered container of data blocks. Data blocks can
**  come from DDL, dictionary or CIF files, where each data block is a
**  container of tables. This class provides methods for construction and
**  destruction, data blocks manipulation (addition, retrieval, renaming.).
**  The class does in-memory management of data blocks, as well as
**  serialization and de-serialization to and from a file. The class supports
**  the following file modes: read-only, create, update and virtual. In
**  read-only mode, blocks and tables can only be read (from an existing table
**  file that has been previously serialized to a file) and cannot be
**  modified. Create mode is used to create a table file from scratch and
**  add/update blocks and tables in it and serialize it to a file. Update mode
**  is used to update an existing table file (that has been previously
**  serialized to a file). Virtual mode only provides in-memory management of
**  data blocks and is used when object persistency is not needed. Hence, all
**  modes except virtual mode provide association between in-memory data
**  blocks and persistent data blocks stored in a file.
*/
class TableFile
{
  public:
    enum eStatusInd
    {
        eCLEAR_STATUS = 0x0000,
        eDUPLICATE_BLOCKS = 0x0001,
        eUNNAMED_BLOCKS = 0x0002
    };

    /**
    **  Constructs a table file.
    **
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names in blocks. Possible values are case
    **    sensitive and case in-sensitive. If not specified, case sensitive
    **    table names are assumed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table file has virtual file mode.
    **
    **  \exception: None
    */
    TableFile(const Char::eCompareType caseSense = Char::eCASE_SENSITIVE);

    /**
    **  Constructs a table file.
    **
    **  \param[in] fileMode - table file mode. Possible values are
    **    read-only, create, update and virtual.
    **  \param[in] fileName - relative or absolute name of the file
    **    where object persistency is maintained. If \e fileMode specifies
    **    virtual mode, this parameter is ignored.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names in blocks. Possible values are case
    **    sensitive and case in-sensitive. If not specified, case sensitive
    **    table names are assumed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    TableFile(const eFileMode fileMode, const string& fileName,
      const Char::eCompareType caseSense = Char::eCASE_SENSITIVE);

    /**
    **  Destructs a table file, by first flushing all the modified tables in
    **  data blocks (for create mode or update mode) and then releasing all
    **  in-memory objects.
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
    virtual ~TableFile();

    /**
    **  Retrieves the name of the file that persistently holds data blocks
    **  and their tables.
    **
    **  \param: None
    **
    **  \return String that contains the file name.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline string GetFileName(void);

    /**
    **  Retrieves table file mode.
    **
    **  \param: None
    **
    **  \return READ_MODE - if read-only mode
    **  \return CREATE_MODE - if create mode
    **  \return UPDATE_MODE - if update mode
    **  \return VIRTUAL_MODE - if virtual mode
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline eFileMode GetFileMode(void);

    /**
    **  Retrieves case sensitivity of table names in blocks.
    **
    **  \param: None
    **
    **  \return eCASE_SENSITIVE - if case sensitive
    **  \return eCASE_INSENSITIVE - if case in-sensitive
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline Char::eCompareType GetCaseSensitivity(void);

    /**
    **  Retrieves table file status in form of one or more flags.
    **
    **  \param: None
    **
    **  \return One or more of these flags: \n
    **    eCLEAR_STATUS - no flag value indicates that there are no flags set \n
    **    eDUPLICATE_BLOCKS - flag that indicates existence of blocks with
    **    the same name, which are internally stored with different names \n
    **    eUNNAMED_BLOCKS - flag that indicates existence of blocks with
    **    empty names
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetStatusInd(void);

    /**
    **  Retrieves the number of data blocks in the table file.
    **
    **  \param: None
    **
    **  \return The number of data blocks in the table file.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetNumBlocks();

    /**
    **  Retrieves data block names.
    **
    **  \param[out] blockNames - retrieved data block names
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void GetBlockNames(vector<string>& blockNames);

    /**
    **  Retrieves the name of the first data block.
    **
    **  \param: None
    **
    **  \return String that contains the name of the first data block.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    string GetFirstBlockName();

    /**
    **  Checks for data block existence.
    **
    **  \param[in] blockName - the name of the data block
    **
    **  \return true - if data block exists
    **  \return false - if data block does not exist
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool IsBlockPresent(const string& blockName);

    /**
    **  Adds a block to the table file. If a block with the specified name
    **  already exists, table file stores it under different internal name,
    **  that is obtained by appending a "#" symbol and the current block
    **  count. After writing blocks, with these kinds of block names,
    **  to an ASCII file, "#" symbol becomes a comment and the text after
    **  it is ignored. This enables the preservation of all duplicate blocks
    **  in the written file.
    **
    **  \param[in] blockName - the name of the data block
    **
    **  \return String that contains the internally assigned data block name.
    **    This value is different from \e blockName, if data block with
    **    the name \e blockName, already exists when this method is invoked.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    string AddBlock(const string& blockName);

    /**
    **  Retrieves a reference to the data block in the table file.
    **
    **  \param[in] blockName - the name of the data block
    **
    **  \return Reference to the data block in the table file.
    **
    **  \pre Data block with name \e blockName must be present
    **
    **  \post None
    **
    **  \exception NotFoundException - if data block with name \e blockName
    **    does not exist
    */
    Block& GetBlock(const string& blockName);

    /**
    **  Changes the data block name.
    **
    **  \param[in] oldBlockName - the name of the data block which is to
    **    be renamed
    **  \param[in] newBlockName - the new data block name
    **
    **  \return String that contains the internally assigned data block name.
    **    This value is different from \e newBlockName, if data block with
    **    the name \e newBlockName, already exists when this method is invoked.
    **
    **  \pre Table file must have at least one data block.
    **  \pre Data block with name \e oldBlockName must be present
    **
    **  \post None
    **
    **  \exception EmptyContainerException - if table file has no data blocks
    **  \exception NotFoundException - if data block with name \e oldBlockName
    **    does not exist
    */
    string RenameBlock(const string& oldBlockName, const string& newBlockName);

    /**
    **  Changes the name of the first data block in table file.
    **
    **  \param[in] newBlockName - the new data block name
    **
    **  \return String that contains the internally assigned data block name.
    **    This value is different from \e newBlockName, if data block with
    **    the name \e newBlockName, already exists when this method is invoked.
    **
    **  \pre Table file must have at least one data block.
    **
    **  \post None
    **
    **  \exception EmptyContainerException - if table file has no data blocks
    */
    inline string RenameFirstBlock(const string& newBlockName);

    /**
    **  Writes only the new or modified tables in data blocks to the
    **  associated persistent storage file (specified at table file
    **  construction time).
    **
    **  \param: None
    **
    **  \return None
    **
    **  \pre Table file must be in create or update mode
    **
    **  \post None
    **
    **  \exception FileModeException - if table file is not in create or
    **    update mode
    */
    void Flush();

    /**
    **  Writes all the data blocks and their tables in the specified file.
    **  The inteded purpose of this method is to write to a file different
    **  than the one specified at construction time.
    **
    **  \param[in] fileName - relative or absolute name of the file
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Serialize(const string& fileName);

    /**
    **  Flushes the table file (if in create or update mode) and closes
    **  the associated persistent storage file.
    **
    **  \param: None
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Close();

  protected:
    string _fileName;

    eFileMode _fileMode;

    // Indicates case sensitivity of identifiers
    Char::eCompareType _caseSense;

    // Indicates the current status of the object
    unsigned int _statusInd;

    mapped_ptr_vector<Block, StringLess> _blocks;

    Serializer* _f;

    void _SetStatusInd(const string& blockName);

    void _AddBlock(const string& blockName, Serializer* serP);

    void _GetNumTablesInBlocks(vector<UInt32>& numTablesInBlocks);

    ISTable* _GetTablePtr(const unsigned int blockIndex,
      const unsigned int tableIndex);
    void _GetAllTables();

    unsigned int GetTotalNumTables();
    void GetTableNames(vector<string>& tableNames);

    void GetTablesIndices(vector<unsigned int>& tablesIndices);
    void GetSortedTablesIndices(vector<unsigned int>& tablesIndices);

    void _ReadFileIndex();
    void _ReadFileIndexVersion0();
    void _ReadFileIndexVersion1();
    void _WriteFileIndex(Serializer* serP,
      const vector<unsigned int>& tableLocs);

  private:
    static const string _version;
    void Init();
    void Open(const string& fileName, const eFileMode fileMode);
    unsigned int GetBlockIndexFromTableId(const string& tableId);
    string GetTableNameFromTableId(const string& tableId);
    string MakeInternalBlockName(const string& blockName,
      const unsigned int blockIndex);
    void PrintHeaderInfo();
};


inline void Block::SetName(const string& name)
{
    _name = name;
}


inline const string& Block::GetName() const
{
    return _name;
}


inline string TableFile::GetFileName(void)
{
    return _fileName;
}


inline eFileMode TableFile::GetFileMode(void)
{
    return _fileMode;
}


inline Char::eCompareType TableFile::GetCaseSensitivity(void)
{
    return(_caseSense);
}


inline unsigned int TableFile::GetStatusInd(void)
{
    return _statusInd;
}


inline unsigned int TableFile::GetNumBlocks()
{
    return _blocks.size();
}


inline string TableFile::RenameFirstBlock(const string& newBlockName)
{
    return(RenameBlock(GetFirstBlockName(), newBlockName));
}


#endif // TABLEFILE_H
