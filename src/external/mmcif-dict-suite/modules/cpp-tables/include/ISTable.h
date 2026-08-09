//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file ISTable.h
**
** \brief Header file for ISTable class.
*/


#ifndef ISTABLE_H
#define ISTABLE_H


#include <float.h>

#include <string>
#include <vector>
#include <map>

#include "mapped_vector.h"
#include "mapped_vector.C"
#include "GenString.h"
#include "ITTable.h"
#include "Serializer.h"


typedef std::multimap<std::string, unsigned int, StringLess> tIndex;


/**
**  \class ISTable
**
**  \brief Public class that respresents a two-dimensional table of strings.
**
**  This class represents a two-dimensional table of cells. Each cell holds
**  a text string. The table is identified by its name. Rows are horizontal
**  table entities identified by row indices, which are unsigned integers
**  ranging from zero to the number of rows minus one. Columns are vertical
**  table entities identified by (non-empty) column names. Column names can
**  be case-sensitive (by default) or case-insensitive (customizable during
**  construction). The class provides methods for table construction and
**  destruction, assignment operator, equal operator, column and row based
**  methods for addition, insertion, retrieval, update, deletion, cell
**  based methods for update and retrieval, column search methods and table
**  printing. Table cells are internally stored as vectors of text strings,
**  where vectors can represent either columns (by default) or rows
**  (customizable during construction).
*/
class ISTable
{
  public:
    typedef ITTable::eOrientation eOrientation;

    static const eOrientation eCOLUMN_WISE = ITTable::eCOLUMN_WISE;
    static const eOrientation eROW_WISE = ITTable::eROW_WISE;

    enum eTableDiff
    {
        eNONE = 0,
        eCASE_SENSE,
        eMORE_COLS,
        eLESS_COLS,
        eCOL_NAMES,
        eMORE_ROWS,
        eLESS_ROWS,
        eCELLS,
        // Used only in block diff to indicate missing table in first block
        eMISSING,
        // Used only in block diff to indicate extra table in first block
        eEXTRA
    };

    typedef ITTable::eSearchType eSearchType;

    static const eSearchType eEQUAL = ITTable::eEQUAL;
    static const eSearchType eLESS_THAN = ITTable::eLESS_THAN;
    static const eSearchType eLESS_THAN_OR_EQUAL = ITTable::eLESS_THAN_OR_EQUAL;
    static const eSearchType eGREATER_THAN = ITTable::eGREATER_THAN;
    static const eSearchType eGREATER_THAN_OR_EQUAL =
      ITTable::eGREATER_THAN_OR_EQUAL;

#ifdef VLAD_SECOND_ITTABLE
    enum eSearchType
    {
        eEQUAL = 0,
        eLESS_THAN,
        eLESS_THAN_OR_EQUAL,
        eGREATER_THAN,
        eGREATER_THAN_OR_EQUAL
    };
#endif

    typedef ITTable::eSearchDir eSearchDir;

    static const eSearchDir eFORWARD = ITTable::eFORWARD;
    static const eSearchDir eBACKWARD = ITTable::eBACKWARD;

#ifdef VLAD_SECOND_ITTABLE
    enum eSearchDir
    {
        eFORWARD = 0,
        eBACKWARD
    };
#endif

    static const unsigned char DT_STRING_VAL = 1;
    static const unsigned char DT_INTEGER_VAL = 2;
    // static const unsigned char DT_DOUBLE_VAL = 3;

    // Sets string comparison case sensitive
    static const unsigned char CASE_SENSE = 0x00;
    // Sets string comparison case insensitive
    static const unsigned char CASE_INSENSE = 0x01;
    // Sets string comparison to be sensitive to whitespace
    static const unsigned char W_SPACE_SENSE = 0x00;
    // Sets string comparison to ignore repeating whitspace.
    // Also ignores leading and trailing whitespace
    static const unsigned char W_SPACE_INSENSE = 0x02;
    // string datatype
    static const unsigned char DT_STRING  = DT_STRING_VAL  << 4;
    // integer datatype
    static const unsigned char DT_INTEGER = DT_INTEGER_VAL << 4;
    // VLAD FEATURE NOT WORKING double is not working, maybe integer. check it      // double datatype
    // static const unsigned char DT_DOUBLE  = DT_DOUBLE_VAL  << 4;

    /**
    **  Constructs a table.
    **
    **  \param[in] colCaseSense - optional parameter that indicates case
    **    sensitivity of column names. Possible values are case sensitive and
    **    case in-sensitive. If not specified, a table with case sensitive
    **    column names is constructed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table has 0 rows and 0 columns.
    **  \post Constructed table is nameless (its name is an empty string).
    **
    **  \exception: None
    */
    ISTable(const Char::eCompareType colCaseSense = Char::eCASE_SENSITIVE);

    /**
    **  Constructs a table.
    **
    **  \param[in] orient - table orientation. Possible values are
    **    row-wise orientation (vectors of strings represent table rows) and
    **    column-wise orientation (vectors of strings represent table columns)
    **  \param[in] colCaseSense - optional parameter that indicates case
    **    sensitivity of column names. Possible values are case sensitive and
    **    case in-sensitive. If not specified, a table with case sensitive
    **    column names is constructed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table has 0 rows and 0 columns.
    **  \post Constructed table is nameless (its name is an empty string).
    **
    **  \exception: None
    */
    ISTable(eOrientation orient, const Char::eCompareType
      colCaseSense = Char::eCASE_SENSITIVE);

    /**
    **  Constructs a table.
    **
    **  \param[in] name - the name of the table to be constructed
    **  \param[in] colCaseSense - optional parameter that indicates case
    **    sensitivity of column names. Possible values are case sensitive and
    **    case in-sensitive. If not specified, a table with case sensitive
    **    column names is constructed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table has 0 rows and 0 columns.
    **
    **  \exception: None
    */
    ISTable(const std::string& name,
      const Char::eCompareType colCaseSense = Char::eCASE_SENSITIVE);

    /**
    **  Constructs a table.
    **
    **  \param[in] name - the name of the table to be constructed
    **  \param[in] orient - table orientation. Possible values are
    **    row-wise orientation (vectors of strings represent table rows) and
    **    column-wise orientation (vectors of strings represent table columns)
    **  \param[in] colCaseSense - optional parameter that indicates case
    **    sensitivity of column names. Possible values are case sensitive and
    **    case in-sensitive. If not specified, a table with case sensitive
    **    column names is constructed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table has 0 rows and 0 columns.
    **
    **  \exception: None
    */
    ISTable(const std::string& name, eOrientation orient,
      const Char::eCompareType colCaseSense = Char::eCASE_SENSITIVE);

    /**
    **  Constructs a table by copying from another table (copy constructor).
    **
    **  \param[in] inTable - reference to a table that will be copied to
    **    the newly constructed table
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table is a clone (has the same name, internal
    **    cells orientation, case sensitivity, column names, content)
    **    as the table referenced by \e inTable.
    **
    **  \exception: None
    */
    ISTable(const ISTable& inTable);

    /**
    **  Destructs a table, by releasing all the used resources.
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
    ~ISTable();

    /**
    **  Copies a table to another table (assignment operator).
    **
    **  \param[in] inTable - reference to the source table
    **
    **  \return Reference to the destination table
    **
    **  \pre None
    **
    **  \post Constructed table is a clone (has the same name, internal
    **    cells orientation, case sensitivity, column names, content)
    **    as the table referenced by \e inTable.
    **
    **  \exception: None
    */
    ISTable& operator=(const ISTable& inTable);

    /**
    **  Compares a table to another table, ignoring the table name.
    **
    **  \param[in] inTable - reference to input table
    **
    **  \return eNONE - if tables are identical
    **  \return eCASE_SENSE - if tables have different column name case
    **    sensitivity
    **  \return eMORE_COLS - if this table has more columns than input table
    **  \return eLESS_COLS - if this table has less columns than input table
    **  \return eCOL_NAMES - if tables have different column names
    **  \return eMORE_ROWS - if this table has more rows than input table
    **  \return eLESS_ROWS - if this table has less rows than input table
    **  \return eCELLS - if tables have different content
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    eTableDiff operator==(ISTable& inTable);

    /**
    **  Retrieves table name.
    **
    **  \param: None
    **
    **  \return Constant reference to a string that contains table name.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline const std::string& GetName() const;

    /**
    **  Changes the table name.
    **
    **  \param[in] name - the new name of the table
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetName(const std::string& name);

    /**
    **  Retrieves the number of columns in the table.
    **
    **  \param: None
    **
    **  \return The number of columns in the table.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetNumColumns() const;

    /**
    **  Retrieves column names.
    **
    **  \param[out] colNames - retrieved column names
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    const std::vector<std::string>& GetColumnNames() const;

    /**
    **  Checks for column existence.
    **
    **  \param[in] colName - the name of the column
    **
    **  \return true - if column exists
    **  \return false - if column does not exist
    **
    **  \pre \e colName must be non-empty
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    */
    bool IsColumnPresent(const std::string& colName);

    /**
    **  Adds a column to the end of the table.
    **
    **  \param[in] colName - the name of the column to be added
    **  \param[in] col - optional parameter that contains the values which
    **    are to be used to fill in the newly added column. If \e col is
    **    specified, filling starts at row index 0 and continues until size
    **    of \e col. If \e col is not specified, the newly added column is
    **    filled with empty values, where filling starts at row index 0 and
    **    ends at row index "number of rows - 1".
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must not be present
    **  \pre If \e col is specified, the size of \e col must be less than or
    **    equal to the number of rows.
    **
    **  \post If table is empty (0 rows) and \e col is specified, the number of
    **    rows is set to the size of \e col. Otherwise, the number of rows is
    **    unchanged.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception AlreadyExistsException - if column with name \e colName
    **    already exists
    **  \exception out_of_range - if size of \e col is greater than the number
    **    of rows
    */
    void AddColumn(const std::string& colName,
      const std::vector<std::string>& col = std::vector<std::string>());

    /**
    **  Inserts a new column at the specified existing column and shifts,
    **    to the right by one, the specified existing column and all columns
    **    after it.
    **
    **  \param[in] colName - the name of the column to be inserted
    **  \param[in] atColName - the name of the column at which the
    **    new column is to be inserted
    **  \param[in] col - optional parameter that contains the values which
    **    are to be used to fill in the newly inserted column. If \e col is
    **    specified, filling starts at row index 0 and continues until size
    **    of \e col. If \e col is not specified, the newly inserted column is
    **    filled with empty values, where filling starts at row index 0 and
    **    ends at row index "number of rows - 1".
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must not be present
    **  \pre \e atColName must be non-empty
    **  \pre Column with name \e atColName must be present
    **  \pre If \e col is specified, the size of \e col must be less than or
    **    equal to the number of rows.
    **  \pre The column which comes, in order, before the column with name
    **    \e atColName, must be non-empty. This is to prevent creation of
    **    non-rectangular tables.
    **
    **  \post If table is empty (0 rows) and \e col is specified, the number of
    **    rows is set to the size of \e col. Otherwise, the number of rows is
    **    unchanged.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception AlreadyExistsException - if column with name \e colName
    **    already exists
    **  \exception EmptyValueException - if \e atColName is empty
    **  \exception NotFoundException - if column with name \e atColName
    **    does not exist
    **  \exception out_of_range - if size of \e col is greater than the number
    **    of rows
    **  \exception out_of_range - if column, which comes, in order, before the
    **    column with name \e atColName, is empty.
    */
    void InsertColumn(const std::string& colName,
      const std::string& afColName, const std::vector<std::string>& col =
      std::vector<std::string>());

    /**
    **  Fills a column with values.
    **
    **  \param[in] colName - the name of the column to be filled
    **  \param[in] col - contains the values which are to be used for filling.
    **    Filling starts at row index 0 and continues until size of \e col.
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **  \pre The size of \e col must be less than or equal to the number of
    **    rows.
    **  \pre The column which comes, in order, before the column with name
    **    \e colName, must be non-empty. This is to prevent creation of
    **    non-rectangular tables.
    **
    **  \post If table is empty (0 rows), the number of rows is set to the
    **    size of \e col. Otherwise, the number of rows is unchanged.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    **  \exception out_of_range - if size of \e col is greater than the number
    **    of rows
    **  \exception out_of_range - if column, which comes, in order, before the
    **    column with name \e colName, is empty.
    */
    void FillColumn(const std::string& colName,
      const std::vector<std::string>& col);

    /**
    **  Retrieves column values.
    **
    **  \param[out] col - retrieved column values
    **  \param[in] colName - the name of the column which content is to be
    **    retrieved.
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    void GetColumn(std::vector<std::string>& col, const std::string& colName);

    /**
    **  Retrieves column values in the specified row range.
    **
    **  \param[out] col - retrieved values
    **  \param[in] colName - the name of the column which content is to be
    **    retrieved.
    **  \param[in] fromRowIndex - the row index of the first cell in the
    **    column to be retrieved.
    **  \param[in] toRowIndex - the row index of the last cell in the column
    **    to be retrieved.
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **  \pre \e fromRowIndex must be less than or equal to the column length
    **  \pre \e toRowIndex must be less than or equal to the column length
    **  \pre \e fromRowIndex must be less than or equal to the \e toRowIndex
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    **  \exception out_of_range - if \e fromRowIndex is greater than the column
    **    length
    **  \exception out_of_range - if \e toRowIndex is greater than the column
    **    length
    **  \exception out_of_range - if \e fromRowIndex is greater than
    **    \e toRowIndex
    */
    void GetColumn(std::vector<std::string>& col, const std::string& colName,
      const unsigned int fromRowIndex, unsigned int toRowIndex);

    /**
    **  Retrieves column values in the specified rows.
    **
    **  \param[out] col - retrieved values
    **  \param[in] colName - the name of the column which content is to be
    **    retrieved
    **  \param[in] rowIndices - row indices of column cells to be retrieved
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **  \pre Row indices in \e rowIndices must be less than or equal to the
    **    column length
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    **  \exception out_of_range - if at least one row index in \e rowIndices is
    **    greater than the column length
    */
    void GetColumn(std::vector<std::string>& col, const std::string& colName,
      const std::vector<unsigned int>& rowIndex);

    /**
    **  Changes the column name.
    **
    **  \param[in] oldColName - the name of the column which is to be renamed
    **  \param[in] newColName - the new column name
    **
    **  \return None
    **
    **  \pre \e oldColName must be non-empty
    **  \pre Column with name \e oldColName must be present
    **  \pre \e newColName must be non-empty
    **  \pre Column with name \e newColName must not be present
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e oldColName is empty
    **  \exception NotFoundException - if column with name \e oldColName
    **    does not exist
    **  \exception EmptyValueException - if \e newColName is empty
    **  \exception AlreadyExistsException - if column with name \e newColName
    **    already exists
    */
    void RenameColumn(const std::string& oldColName,
      const std::string& newColName);

    /**
    **  Sets all cells in the column to empty string.
    **
    **  \param[in] colName - the name of the column
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post Column length is unchanged.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    void ClearColumn(const std::string& colName);

    /**
    **  Deletes a column from the table.
    **
    **  \param[in] colName - the name of the column
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post The number of table columns is reduced by one.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    void DeleteColumn(const std::string& colName);

    /**
    **  Retrieves the number of rows in the table.
    **
    **  \param: None
    **
    **  \return The number of rows in the table.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetNumRows() const;

    /**
    **  Adds a row to the end of the table. For an empty table (0 rows),
    **  the number of inserted cells is equal to the number of table columns.
    **  For a non-empty table, the number of inserted cells is equal to the
    **  number of non-empty columns (this is in order to prevent creation of
    **  non-rectangular tables). The newly added row is, optionally, filled
    **  with values, starting at the first column.
    **
    **  \param[in] row - optional parameter that contains the values which
    **    are to be used to fill in the newly added row. If \e row is
    **    specified, filling starts at the first column and continues until
    **    size of \e row is reached. If \e row is not specified and table is
    **    empty, the newly inserted row is filled with empty values,
    **    where filling starts at the first column and continues until the
    **    number of columns is reached.  If \e row is not specified and table
    **    is not empty, the newly inserted row is filled with empty values,
    **    where filling starts at the first column and continues until the
    **    number of non-empty columns is reached.
    **
    **  \return The new number of rows after the row addition.
    **
    **  \pre Table must have at least one column, which can be empty.
    **  \pre If table is not empty and \e row is specified, the size of \e row
    **    must be less than or equal to the number of non-empty columns.
    **    This is in order to prevent creation of non-rectangular tables.
    **  \pre If table is empty and \e row is specified, the size of \e row must
    **    be less than or equal to the number of columns.
    **
    **  \post The number of rows is increased by one.
    **
    **  \exception EmptyContainerException - if table has no columns.
    **  \exception out_of_range - if table is not empty and size of \e row is
    **    greater than the number of non-empty columns.
    **  \exception out_of_range - if table is empty and size of \e row is
    **    greater than the number of columns.
    */
    unsigned int AddRow(const std::vector<std::string>& row =
      std::vector<std::string>());

    /**
    **  Inserts a row at the specified row index and shifts, down by one,
    **  the old row and all other rows below it. For an empty table (0 rows),
    **  the number of inserted cells is equal to the number of table columns.
    **  For a non-empty table, the number of inserted cells is equal to the
    **  number of non-empty columns (this is in order to prevent creation of
    **  non-rectangular tables). The newly inserted row is optionally filled
    **  with values, starting at the first column.
    **
    **  \param[in] atRowIndex - index of the row at which the new row is to be
    **    inserted. Note: If \e atRowIndex is equal to the number of rows, the
    **    operation of this method is equivalent to AddRow().
    **  \param[in] row - optional parameter that contains the values which
    **    are to be used to fill in the newly inserted row. If \e row is
    **    specified, filling starts at the first column and continues until
    **    size of \e row is reached. If \e row is not specified and table is
    **    empty, the newly inserted row is filled with empty values,
    **    where filling starts at the first column and continues until the
    **    number of columns is reached.  If \e row is not specified and table
    **    is not empty, the newly inserted row is filled with empty values,
    **    where filling starts at the first column and continues until the
    **    number of non-empty columns is reached.
    **
    **  \return The new number of rows after the row insertion.
    **
    **  \pre Table must have at least one column, which can be empty.
    **  \pre \e atRowIndex must be greater than or equal to 0 and less than
    **    or equal to the number of table rows.
    **  \pre If table is not empty and \e row is specified, the size of \e row
    **    must be less than or equal to the number of non-empty columns.
    **    This is in order to prevent creation of non-rectangular tables.
    **  \pre If table is empty and \e row is specified, the size of \e row must
    **    be less than or equal to the number of columns.
    **
    **  \post The number of rows is increased by one.
    **  \post Row indices, of the rows below the inserted row, are invalidated
    **    by being increased by one.
    **
    **  \exception EmptyContainerException - if table has no columns.
    **  \exception out_of_range - if \e atRowIndex is greater than the number
    **    of table rows.
    **  \exception out_of_range - if table is not empty and size of \e row is
    **    greater than the number of non-empty columns.
    **  \exception out_of_range - if table is empty and size of \e row is
    **    greater than the number of columns.
    */
    unsigned int InsertRow(const unsigned int atRowIndex,
      const std::vector<std::string>& row = std::vector<std::string>());

    /**
    **  Fills, with values, a row at the specified row index, starting at the
    **  the first column.
    **
    **  \param[in] rowIndex - index of the row that is to be filled.
    **  \param[in] row - values which are to be used to fill in the row.
    **    Filling starts at the first column and continues until size of
    **    \e row is reached.
    **
    **  \return None
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **  \pre The size of \e row must be less than or equal to the number of
    **    non-empty columns. This is in order to prevent creation of
    **    non-rectangular tables.
    **
    **  \post None
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal to
    **    the number of table rows.
    **  \exception out_of_range - if size of \e row is greater than the number
    **    of non-empty columns.
    */
    void FillRow(const unsigned int rowIndex,
      const std::vector<std::string>& row);

    /**
    **  Retrieves row values.
    **
    **  \param[out] row - retrieved row values
    **  \param[in] rowIndex - index of the row which values are to be
    **    retrieved.
    **  \param[in] fromColName - optional parameter which specifies the
    **    location of the first cell to be retrieved. If not specified
    **    the first column cell is used.
    **  \param[in] toColName - optional parameter which specifies the
    **    location of the last cell to be retrieved. If not specified
    **    the last non-empty-column cell is used.
    **
    **  \return None
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **  \pre If \e fromColName is specified, the column with name
    **    \e fromColName must be present and must be non-empty
    **  \pre If \e toColName is specified, the column with name
    **    \e toColName must be present and must be non-empty
    **  \pre If \e fromColName is different than \e toColName, it must come
    **    prior to it in the column order.
    **
    **  \post None
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal to
    **    the number of table rows.
    **  \exception NotFoundException - If \e fromColName is specified and
    **    column with name \e fromColName does not exist
    **  \exception NotFoundException - If \e toColName is specified and
    **    column with name \e toColName does not exist
    **  \exception out_of_range - If \e fromColName is specified and
    **    column with name \e fromColName exists but is empty
    **  \exception out_of_range - If \e toColName is specified and
    **    column with name \e toColName exists but is empty
    **  \exception out_of_range - if \e fromColName is different than
    **    \e toColName and it comes after it in the column order.
    */
    void GetRow(std::vector<std::string>& row, const unsigned int rowIndex,
      const std::string& fromColName = std::string(),
      const std::string& toColName = std::string());

    /**
    **  Retrieves a constant reference to a row of values the table.
    **
    **  \param[in] rowIndex - index of a row to which a reference is to be
    **    retrieved.
    **
    **  \return Constant reference to the row of values in the table.
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **
    **  \post None
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal to
    **    the number of table rows.
    */
    const std::vector<std::string>& GetRow(const unsigned int rowIndex);

    /**
    **  Sets all cells in the row to empty string.
    **
    **  \param[in] rowIndex - index of the row that is to be cleared.
    **
    **  \return None
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **
    **  \post Number of table rows is unchanged.
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal to
    *     the number of table rows.
    */
    void ClearRow(const unsigned int rowIndex);

    /**
    **  Deletes a row with the specified index and shifts, up by one,
    **  all other rows below it.
    **
    **  \param[in] rowIndex - index of the row that is to be deleted.
    **
    **  \return None
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **
    **  \post Number of table rows is reduced by one.
    **  \post Row indices of the rows which are below the deleted row are
    **    invalidated by being reduced by one.
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal to
    **    the number of table rows.
    */
    void DeleteRow(const unsigned int rowIndex);

    /**
    **  Deletes rows with specified indices.
    **
    **  \param[in] rows - indices of rows that are to be deleted.
    **
    **  \return None
    **
    **  \pre indices in \e rows must be greater than or equal to 0 and less
    **    than the number of table rows.
    **
    **  \post Number of table rows is reduced by the size of \e rows.
    **  \post Row indices of the remaining rows are invalidated by being
    **    appropriatelly adjusted.
    **
    **  \exception out_of_range - if any row index in \e rows is greater than
    **    or equal to the number of table rows.
    */
    void DeleteRows(const std::vector<unsigned int>& rows);

    /**
    **  Retrieves the row index of the last row in the table.
    **
    **  \param: None
    **
    **  \return The index of the last row in the table.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetLastRowIndex();

    /**
    **  Updates a cell in the table.
    **
    **  \param[in] rowIndex - row index of the cell that is to be updated
    **  \param[in] colName - the name of the column of the cell that is to be
    **    updated
    **  \param[in] value - the new value
    **
    **  \return None
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post None
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal
    **    to the number of table rows.
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    void UpdateCell(const unsigned int rowIndex, const std::string& colName,
      const std::string& value);

    /**
    **  Retrieves a constant reference to the cell in the table.
    **
    **  \param[in] rowIndex - row index of the cell
    **  \param[in] colName - the name of the column of the cell
    **
    **  \return Constant reference to the cell in the table.
    **
    **  \pre \e rowIndex must be greater than or equal to 0 and less than
    **    the number of table rows.
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post None
    **
    **  \exception out_of_range - if \e rowIndex is greater than or equal
    **    to the number of table rows.
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    const std::string& operator()(const unsigned int rowIndex,
      const std::string& colName) const;

    /**
    **  Sets column flags that are only used in column search. These flags
    **  control how cell values in a column are interpreted at the time of
    **  search. They can be interpreted as strings or integers, as
    **   case-sensitive or case-insensitive strings, as space-ignoring
    **  or space-non-ignoring strings. Multiple flags can be specified using
    **  "|" operator when invoking this method.
    **
    **  \param[in] colName - the name of the column
    **  \param[in] flags - column search flags. It can have any or multiple
    **    "or"-ed values of: DT_STRING, DT_INTEGER, CASE_SENSE, CASE_INSENSE,
    **    W_SPACE_SENSE, W_SPACE_INSENSE.
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    void SetFlags(const std::string& colName, const unsigned char flags);

    /**
    **  Retrieves data type flag of a column.
    **
    **  \param[in] colName - the name of the column
    **
    **  \return DT_STRING_VAL - if data type of a column is string.
    **  \return DT_INTEGER_VAL - if data type of a column is integer.
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    unsigned char GetDataType(const std::string& colName);

    /**
    **  Searches the columns for the first occurrence of target values and
    **  returns the row index, where the match was found. The performed search
    **  is a forward search (starts at row index 0) and search criteria is
    **  value equality. If match was not found, the number of table rows
    **  is returned.
    **
    **  \param[in] targets - values that are to be searched for
    **  \param[in] colNames - column names that are to be searched
    **  \param[in] indexName - optional parameter not used and will be soon
    **    be removed
    **
    **  \return the first row index, where the match was found
    **  \return the number of rows, if the match was not found
    **
    **  \pre Each column name in \e colNames must be non-empty
    **  \pre Each column name in \e colNames must be present
    **  \pre \e colNames and \e targets must have the same size
    **
    **  \post None
    **
    **  \exception EmptyValueException - if one or more column names in
    **    \e colNames is empty
    **  \exception NotFoundException - if one or more column names in
    **    \e colNames does not exist
    **  \exception out_of_range - if \e colNames and \e targets have
    **    different sizes
    */
    unsigned int FindFirst(const std::vector<std::string>& targets,
      const std::vector<std::string>& colNames,
      const std::string& indexName = std::string());

    /**
    **  Searches one column for all occurrences of target value and
    **  returns row indices, where the match was found.
    **
    **  \param[out] res - vector of row indices, where the match was found
    **  \param[in] target - value that is to be searched for
    **  \param[in] colName - column name that is to be searched
    **  \param[in] searchType - optional parameter that specifies search
    **    criteria: equality, less than, less than or equal, greater than,
    **    greater than or equal. These are the text strings search criteria.
    **    If not specified, the search criteria is equality.
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    */
    void Search(std::vector<unsigned int>& res, const std::string& target,
      const std::string& colName, const unsigned int fromRowIndex = 0,
      const eSearchDir searchDir = ISTable::eFORWARD,
      const eSearchType searchType = ISTable::eEQUAL);

    /**
    **  Searches the columns for all occurrences of target values and
    **  returns row indices, where the match was found.
    **
    **  \param[out] res - vector of row indices, where the match was found
    **  \param[in] targets - values that are to be searched for
    **  \param[in] colNames - column names that are to be searched
    **  \param[in] searchType - optional parameter that specifies search
    **    criteria: equality, less than, less than or equal, greater than,
    **    greater than or equal. These are the text strings search criteria.
    **    If not specified, the search criteria is equality.
    **  \param[in] indexName - optional parameter not used and will be soon
    **    be removed.
    **
    **  \return None
    **
    **  \pre Each column name in \e colNames must be non-empty
    **  \pre Each column name in \e colNames must be present
    **  \pre \e colNames and \e targets must have the same size
    **
    **  \post None
    **
    **  \exception EmptyValueException - if one or more column names in
    **    \e colNames is empty
    **  \exception NotFoundException - if one or more column names in
    **    \e colNames does not exist
    **  \exception out_of_range - if \e colNames and \e targets have
    **    different sizes
    */
    void Search(std::vector<unsigned int>& res,
      const std::vector<std::string>& targets,
      const std::vector<std::string>& colNames,
      const unsigned int fromRowIndex = 0,
      const eSearchDir searchDir = ISTable::eFORWARD,
      const eSearchType searchType = ISTable::eEQUAL,
      const std::string& indexName = std::string());

    /**
    **  Finds duplicate rows and, optionally, deletes them.
    **
    **  \param[out] duplRows - vector of pairs of indices, where each pair
    **    indicates a row and its duplicate row
    **  \param[in] colNames - column names that are of inerest in determining
    **    duplicate rows. Note that determination of duplicate rows is not
    **    done based on all values in a row, but based on the cell values
    **    in the columns specified in this parameter.
    **  \param[in] keepDuplRows - indicates whether duplicate rows should be
    **    kept in the table (if true) or deleted (if false).
    **  \param[in] searchDir - optional parameter which specifies search
    **    direction. This parameter is only relevant when duplicate rows are
    **    deleted. If \e searchDir specifies forward search, the deleted
    **    duplicate rows will have bigger index than the original row.
    **    If \e searchDir specifies backward search, the deleted duplicate
    **    rows will have smaller index than the original row.
    **
    **  \return None
    **
    **  \pre Each column name in \e colNames must be non-empty
    **  \pre Each column name in \e colNames must be present
    **
    **  \post If deletion of duplicate rows is requested, the number of
    **    table rows will be reduced by the number of duplicate rows.
    **
    **  \exception EmptyValueException - if one or more column names in
    **    \e colNames is empty
    **  \exception NotFoundException - if one or more column names in
    **    \e colNames does not exist
    */
    void FindDuplicateRows(std::vector<std::pair<unsigned int,
      unsigned int> >& duplRows, const std::vector<std::string>& colNames,
      const bool keepDuplRows, const eSearchDir searchDir = ISTable::eFORWARD);

    /**
    **  Retrieves case sensitivity of column names.
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
    inline Char::eCompareType GetColCaseSense() const;

    /**
    **  Utility method, not part of users public API.
    */
    inline void SetModified(const bool modified);

    /**
    **  Utility method, not part of users public API.
    */
    inline bool GetModified();

    /**
    **  Utility method, not part of users public API.
    */
    void SetSerializer(Serializer* ser);

    /**
    **  Utility method, not part of users public API.
    */
    int WriteObject(Serializer* ser, int& size);

    /**
    **  Utility method, not part of users public API.
    */
    int GetObject(UInt32 index, Serializer* ser);

    /**
    **  Utility method, not part of users public API.
    */
    void Read(unsigned int indexInFile);

    /**
    **  Utility method, not part of users public API.
    */
    int Write();

    /**
    **  Utility method, not part of users public API.
    */
    // typeOfMerge is 0 for overwrite, 1 for overlap
    static ISTable* Merge(ISTable& firstTable, ISTable& secondTable,
      unsigned int typeOfMerge = 0);

    /**
    **  Utility method, not part of users public API.
    */
    bool PrintDiff(ISTable& inTable);

    /**
    **  Utility method, not part of users public API.
    */
    inline bool IndexExists(const std::string& indexName);

    /**
    **  Utility method, not part of users public API.
    */
    void CreateIndex(const std::string& indexName,
      const std::vector<std::string>& colNames,
      const unsigned int unique = 0);

    /**
    **  Utility method, not part of users public API.
    */
    void UpdateIndex(const std::string& indexName, const unsigned int rowIndex);

    /**
    **  Utility method, not part of users public API.
    */
    void RebuildIndex(const std::string& indexName);

    /**
    **  Utility method, not part of users public API.
    */
    void RebuildIndices();

    /**
    **  Utility method, not part of users public API.
    */
    void DeleteIndex(const std::string& indexName);

    /**
    **  Utility method, not part of users public API.
    */
    inline unsigned int GetNumIndices();

    /**
    **  Utility method, not part of users public API.
    */
    void CreateKey(const std::vector<std::string>& colNames);

    /**
    **  Utility method, not part of users public API.
    */
    void DeleteKey();

    /**
    **  Utility method, not part of users public API.
    */
    static void SetUnion(const std::vector<unsigned int>& a,
      const std::vector<unsigned int>& b, std::vector<unsigned int>& ret);

    /**
    **  Utility method, not part of users public API.
    */
    static void SetIntersect(const std::vector<unsigned int>& a,
      const std::vector<unsigned int>& b, std::vector<unsigned int>& ret);

    /**
    **  Utility method, not part of users public API.
    */
    void GetColumnsIndices(std::vector<unsigned int>& colIndices,
      const std::vector<std::string>& colNames);

    /**
    **  Utility method, not part of users public API.
    */
    void GetColumn(std::vector<std::string>& col, const std::string& colName,
      const std::string& indexName);

  private:
    static const unsigned int MAX_NUM_ITTABLE_ROWS = 1000;

    // number of digit DBL_MIN_10_EXP, letter e is not included in size
    static const unsigned int EXPONENT      =  4;
    static const unsigned int MAX_PRECISION = DBL_DIG;
    //???DBL_MANT_DIG;
    static const unsigned int MANTISSA       =  MAX_PRECISION + 2;
    static const unsigned int INT_LIMIT      = 11;

    // datatype mask
    static const unsigned char DT_MASK        = 15 << 4;
    // string comparison sensitivity mask
    static const unsigned char SC_MASK        = 0x01;
    // white space sensitivity mask
    static const unsigned char WS_MASK        = 0x02;
    static const unsigned char LAST_DT_VALUE  = 3;
    static const unsigned int  DEFAULT_PRECISION = MAX_PRECISION;
    static const unsigned char DEFAULT_OPTIONS;

    static const std::string _version;

    std::string _name;

    std::vector<ITTable> _ittables;

    ITTable::eOrientation _orient;

    Char::eCompareType _colCaseSense;

    mapped_vector<std::string, StringLess> _colNames;

    std::vector<unsigned int> _precision;
    std::vector<unsigned char> _compare_opts;

    std::vector<std::string> _indexNames;
    std::vector<std::vector<unsigned int> > _listsOfColumns;
    std::vector<unsigned int> _unique;

    Serializer* _ser;

    bool _modified; // Indicates whether table has been modified

    unsigned int _numRows;

    mutable unsigned int _rowIndexCache;
    mutable std::pair<unsigned int, unsigned int> _rowLocCache;

    void InsertColumn(const std::string& colName, const unsigned int atColIndex,
      const std::vector<std::string>& col = std::vector<std::string>());
    void CreateColumn(const std::string& colName, const unsigned int atColIndex,
      const std::vector<std::string>& col = std::vector<std::string>());
    int UpdateCell(const std::string& cell, const unsigned int colIndex,
      const unsigned int rowIndex);
    const std::string& operator()(const unsigned int rowIndex,
      const unsigned int colIndex) const;
    int SetFlags(const unsigned char newOpts, const unsigned int colIndex);
    void FindDuplicateRows(const std::vector<unsigned int>& colIndices,
      std::vector<std::pair<unsigned int, unsigned int> >& duplRows,
      const unsigned int keep, const eSearchDir searchDir = eFORWARD);
    void VerifyColumnsIndices(const std::vector<unsigned int>& colIndices);
    bool AreListsOfColumnsValid(const std::vector<unsigned int>& colIndices);
    void CreateIndex(const std::string& indexName,
      const std::vector<unsigned int>& colIndices,
      const unsigned int unique = 0);
    void CreateKey(const std::vector<unsigned int>& colIndices);
    unsigned int FindFirst(const std::vector<std::string>& targets,
      const std::vector<unsigned int>& colIndices,
      const std::string& indexName = std::string());
    void Search(std::vector<unsigned int>& res,
      const std::vector<std::string>& targets,
      const std::vector<unsigned int>& colIndices,
      const unsigned int fromRowIndex = 0,
      const eSearchDir searchDir = eFORWARD,
      const eSearchType searchType = eEQUAL,
      const std::string& indexName = std::string());

    void Init();
    void Clear();

    Char::eCompareType
      GetCompareType(const std::vector<unsigned int>& colIndices);

    std::string CellValue(const unsigned int colIndex,
      const unsigned int rowIndex);
    std::string ConvertString(const std::string& value,
      const unsigned int colIndex);
    std::string MultiStringsValue(const std::vector<std::string>& values,
      const std::vector<unsigned int>& colIndices);
    std::string SubRowValue(const std::vector<unsigned int>& colIndices,
      const unsigned int rowIndex);
    std::string AggregateRow(const std::vector<unsigned int>& colIndices,
      const unsigned int rowIndex);

    inline void AppendToAndDelimit(std::string& to,
      const std::string& appending);

    void ValidateOptions(unsigned int colIndex);

    std::string CreateInternalIndexName(const unsigned int indexIndex);
    void UpdateIndex(const unsigned int indexIndex,
      const unsigned int rowIndex);
    void RebuildIndex(const unsigned int indexIndex);
    void ClearIndex(const unsigned int indexIndex);
    void DeleteIndex(const unsigned int indexIndex);

    int FindIndex(const std::string& indexName);
    int FindIndex(const std::vector<unsigned int>& colIndices);

    void UpdateIndices(const unsigned int rowIndex);
    void ClearIndices();

    bool IsColumnInIndex(const unsigned int indexIndex,
      const unsigned int colIndex);

    int FindKeyIndex();

    void UpdateColListOnColInsert(const unsigned int colIndex);
    void UpdateColListOnColDelete(const unsigned int colIndex);
    void UpdateColListOnCellUpdate(const unsigned int rowIndex,
      const unsigned int colIndex);

    unsigned int FindFirst(const std::vector<std::string>& targets,
      const std::vector<unsigned int>& colIndices,
      const unsigned int indexIndex);

    int WriteObjectV9(Serializer*, int& size);

    int GetObjectV9(UInt32 index, Serializer*);
    int GetObjectV8(UInt32 index, Serializer*);
    int GetObjectV7(UInt32 index, Serializer*);
    int GetObjectV6(UInt32 index, Serializer*);
    int GetObjectV3(UInt32 index, Serializer*);
    int GetObjectV2(UInt32 index, Serializer*);
    int GetObjectV1(UInt32 index, Serializer*);
    int GetObjectV1_1(UInt32 index, Serializer*);

    void ConvertToInt(const std::string& a, std::string& ret);
    void ConvertDouble(const std::string& a, std::string& ret);
    void ConvertToLowerNoWhiteSpace(const std::string& a, std::string& ret);

    void GetRowLocation(std::pair<unsigned int, unsigned int>& rowLoc,
      const unsigned int rowIndex) const;
    void CacheRowLocation(const unsigned int rowIndex) const;

    void CreateSubtables(const unsigned int numRows);
    void CreateSubtableColumns(const unsigned int colIndex,
      const std::vector<std::string>& col);
    void CreateColumn(const unsigned int atColIndex,
      const std::vector<std::string>& col);

    void Print(const std::string& indexName);

    unsigned int GetColumnIndex(const std::string& colName) const;

};


std::ostream& operator<<(std::ostream& out, const ISTable& isTable);


inline unsigned int ISTable::GetLastRowIndex()
{

    return(GetNumRows() - 1);

}


inline unsigned int ISTable::GetNumIndices()
{

    return(_indexNames.size());

}


inline bool ISTable::IndexExists(const std::string& indexName)
{

    int ret = FindIndex(indexName);

    if (ret == -1)
    {
        return(false);
    }
    else
    {
        return(true);
    }

}


inline void ISTable::AppendToAndDelimit(std::string& to,
  const std::string& appending)
{

    to += appending;
    // VLAD HARDCODED CONST
    to += " ";

}


inline void ISTable::SetModified(const bool modified)
{
    _modified = modified;
}


inline bool ISTable::GetModified()
{
    return _modified;
}


inline const std::string& ISTable::GetName() const
{
    return(_name);
}


inline unsigned int ISTable::GetNumRows() const
{
    return(_numRows);
}


inline unsigned int ISTable::GetNumColumns() const
{
    return(_colNames.size());
}


inline Char::eCompareType ISTable::GetColCaseSense() const
{
    return(_colCaseSense);
}


#endif // ISTABLE_H
