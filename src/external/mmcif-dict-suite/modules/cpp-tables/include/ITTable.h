//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file ITTable.h
**
** \brief Header file for ITTable class.
*/


#ifndef ITTABLE_H
#define ITTABLE_H


#include <float.h>

#include <string>
#include <vector>
#include <map>

#include "GenString.h"
#include "TTable.h"
#include "Serializer.h"


typedef std::multimap<std::string, unsigned int, StringLess> tIndex;


/**
**  \class ITTable
** 
**  \brief Private class that respresents a two-dimensional table of strings.
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
class ITTable
{
  public:
    enum eOrientation
    {
        eCOLUMN_WISE = 0,
        eROW_WISE
    };

    enum eSearchType
    {
        eEQUAL = 0,
        eLESS_THAN,
        eLESS_THAN_OR_EQUAL,
        eGREATER_THAN,
        eGREATER_THAN_OR_EQUAL
    };

    enum eSearchDir
    {
        eFORWARD = 0,
        eBACKWARD
    };

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
    **  \return Not applicable
    **  
    **  \pre None
    **  
    **  \post Constructed table has 0 rows and 0 columns.
    **  \post Constructed table is nameless (its name is an empty string).
    **
    **  \exception: None
    */
    ITTable();

    /** 
    **  Constructs a table.
    **  
    **  \param[in] orient - table orientation. Possible values are
    **    row-wise orientation (vectors of strings represent table rows) and
    **    column-wise orientation (vectors of strings represent table columns)
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
    ITTable(eOrientation orient);

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
    ITTable(const ITTable& inTable);

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
    ~ITTable();
 
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
    ITTable& operator=(const ITTable& inTable);

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
      const std::string& atColName, const std::vector<std::string>& col =
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
    **  Appends one cell to the specified column. 
    **
    **  \param[in] colName - the name of the column to which the cell is to
    **   be appended
    **  \param[in] cell - contains the value that is to be appended
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **  \pre The column, which comes, in order, before the column with name
    **    \e colName, must be non-empty. This is to prevent creation of
    **    non-rectangular tables.
    **
    **  \post The number of rows is increased by one. 
    **  \post Cells in other columns of the, newly appended row, are set to
    **    empty values.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    **  \exception out_of_range - if column, which comes, in order, before the
    **    column with name \e colName, is empty.
    */
    void AppendToColumn(const std::string& colName, const std::string& cell);

    /**
    **  Appends cells to the specified column. 
    **
    **  \param[in] colName - the name of the column to which the cells are to
    **   be appended
    **  \param[in] col - contains the values which are to be appended.
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Column with name \e colName must be present
    **  \pre The column, which comes, in order, before the column with name
    **    \e colName, must be non-empty. This is to prevent creation of
    **    non-rectangular tables.
    **
    **  \post The number of rows is increased by the size of \e col. 
    **  \post Cells in other columns of the, newly appended rows, are set to
    **    empty values.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if column with name \e colName
    **    does not exist
    **  \exception out_of_range - if column, which comes, in order, before the
    **    column with name \e colName, is empty.
    */
    void AppendToColumn(const std::string& colName,
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
    **  Retrieves a reference to the cell in the table.
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
#ifdef VLAD_SECOND_ITTABLE
    unsigned int FindFirst(const std::vector<std::string>& targets,
      const std::vector<std::string>& colNames,
      const std::string& indexName = std::string());
#endif

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
#ifdef VLAD_SECOND_ITTABLE
    void Search(std::vector<unsigned int>& res, const std::string& target,
      const std::string& colName, const eSearchType searchType = eEQUAL);
#endif

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
#ifdef VLAD_SECOND_ITTABLE
    void Search(std::vector<unsigned int>& res,
      const std::vector<std::string>& targets,
      const std::vector<std::string>& colNames,
      const eSearchType searchType = eEQUAL,
      const std::string& indexName = std::string());
#endif
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
      unsigned int> >& duplRows,
      const std::vector<std::string>& colNames, const bool keepDuplRows,
      const eSearchDir searchDir = eFORWARD);

    /* \todo Figure this out. */
    void ValidateOptions(unsigned int colIndex);
    void UpdateIndex(const unsigned int indexIndex,
      const unsigned int rowIndex);
    void InsertIndexEntry(const unsigned int indexIndex,
      const unsigned int rowIndex);
    void DeleteIndexEntry(const unsigned int indexIndex,
      const unsigned int rowIndex);
    void VerifyColumnsIndices(const std::vector<unsigned int>& colIndices);
    int FindIndex(const std::vector<unsigned int>& colIndices);
    int SetFlags(const unsigned char newOpts, const unsigned int colIndex);
    unsigned int FindFirst(const std::vector<std::string>& targets,
      const std::vector<unsigned int>& colIndices, const unsigned int indexIndex);

    void Search(std::vector<unsigned int>& res,
      const std::vector<std::string>& targets,
      const std::vector<unsigned int>& colIndices,
      const unsigned int indexIndex,
      const eSearchType searchType = eEQUAL);
    void DeleteIndex(const unsigned int indexIndex);
    void Search(std::vector<unsigned int>& res, const std::string& target,
      const unsigned int colIndex, const eSearchType searchType = eEQUAL);
    void FindDuplicateRows(const std::vector<unsigned int>& colIndices,
      std::vector<std::pair<unsigned int, unsigned int> >& duplRows,
      const bool keep, const eSearchDir searchDir = eFORWARD);
    void RebuildIndex(const unsigned int indexIndex);
    void InsertColumn(const unsigned int colIndex,
      const std::vector<std::string>& col = std::vector<std::string>());
    void InsertColumn(const unsigned int colIndex,
      std::vector<std::string>::const_iterator colBeg,
      std::vector<std::string>::const_iterator colEnd);
    void Clear();
    const std::string& operator()(const unsigned int rowIndex,
      const unsigned int colIndex) const;
    int UpdateCell(const std::string& cell, const unsigned int colIndex,
      const unsigned int rowIndex);
    void FillColumn(const std::vector<std::string>& col,
      const unsigned int colIndex);
    void FillColumn(const unsigned int colIndex,
      std::vector<std::string>::const_iterator colBeg,
      std::vector<std::string>::const_iterator colEnd);
    void AppendToColumn(const unsigned int colIndex,
      const std::vector<std::string>& col);
    void AppendToColumn(const unsigned int colIndex, const std::string& cell);
    void CreateColumn(const unsigned int atColIndex,
      const std::vector<std::string>& col = std::vector<std::string>());
    void CreateColumn(const unsigned int atColIndex,
      std::vector<std::string>::const_iterator colBeg,
      std::vector<std::string>::const_iterator colEnd);
    void GetColumn(std::vector<std::string>& col, const unsigned int colIndex,
      const unsigned int fromRowIndex, unsigned int toRowIndex);
    void GetColumn(std::vector<std::string>& col, const unsigned int colIndex,
      const std::vector<unsigned int>& rowIndex);
    void ClearColumn(const unsigned int colIndex);
    void DeleteColumn(const unsigned int colIndex);
    void GetColumn(std::vector<std::string>& col, const unsigned int colIndex);
    void GetRow(std::vector<std::string>& row, const unsigned int rowIndex,
      const unsigned int fromColIndex, unsigned int toColIndex);
    const std::vector<std::string>& GetRow(const unsigned int rowIndex);
    eOrientation GetOrientation();
    void CreateIndex(const std::vector<unsigned int>& colIndices,
      const unsigned int unique = 0);
    
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
    int Read(unsigned int indexInFile, Serializer* ser);

    /**
    **  Utility method, not part of users public API.
    */
    int Write(Serializer* ser, int& size);

    /**
    **  Utility method, not part of users public API.
    */
    void RebuildIndices();

    void InsertEntry(const unsigned int rowIndex);
    void DeleteEntry(const unsigned int rowIndex);

    /**
    **  Utility method, not part of users public API.
    */
    inline unsigned int GetNumIndices();

    /**
    **  Utility method, not part of users public API.
    */
    void GetColumn(std::vector<std::string>& col, const unsigned int colIndex,
      const unsigned int indexIndex);

  private:
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

    // static const string _version;

    TTable _ttable;

    eOrientation _orient;

    Serializer* _ser;

    std::vector<unsigned char> _compare_opts;

    std::vector<std::vector<unsigned int> > _listsOfColumns;
    std::vector<unsigned int> _unique;
    std::vector<tIndex> _indices;

    bool AreListsOfColumnsValid(const std::vector<unsigned int>& colIndices);
    void CreateKey(const std::vector<unsigned int>& colIndices);

    void Init();

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

    std::string CreateInternalIndexName(const unsigned int indexIndex);
    void ClearIndex(const unsigned int indexIndex);

    void UpdateIndices(const unsigned int rowIndex);
    void ClearIndices();

    bool IsColumnInIndex(const unsigned int indexIndex,
      const unsigned int colIndex);

    int FindKeyIndex();

    void UpdateColListOnColInsert(const unsigned int colIndex);
    void UpdateColListOnColDelete(const unsigned int colIndex);
    void UpdateIndicesOnCellUpdate(const unsigned int rowIndex,
      const unsigned int colIndex);

    void ConvertToInt(const std::string& a, std::string& ret);
    void ConvertDouble(const std::string& a, std::string& ret);
    void ConvertToLowerNoWhiteSpace(const std::string& a, std::string& ret);

    void Print(unsigned int indexIndex);
};


std::ostream& operator<<(std::ostream& out, const ITTable& isTable);


inline unsigned int ITTable::GetLastRowIndex()
{

    return(GetNumRows() - 1);

}


inline unsigned int ITTable::GetNumIndices()
{

    return(_listsOfColumns.size());

}


inline void ITTable::AppendToAndDelimit(std::string& to,
  const std::string& appending)
{

    to += appending;
    // VLAD HARDCODED CONST
    to += " ";

}


inline unsigned int ITTable::GetNumColumns() const
{
    if (_orient == eCOLUMN_WISE)
        return(_ttable.GetNumTuples());
    else
        return(_ttable.GetNumColumns());
}

inline unsigned int ITTable::GetNumRows() const
{
    if (_orient == eCOLUMN_WISE)
        return(_ttable.GetNumColumns());
    else
        return(_ttable.GetNumTuples());
}

#endif // ITTABLE_H
