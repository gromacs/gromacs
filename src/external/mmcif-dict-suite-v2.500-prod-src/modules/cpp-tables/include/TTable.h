//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file TTable.h
**
** \brief Header file for TTable class.
*/


#ifndef TTABLE_H
#define TTABLE_H


#include <string>
#include <vector>

#include "TableError.h"
#include "Serializer.h"

#define TTABLE_COLUMN_DELETE_AS_REMOVE


/**
** \class TTable
**
** \brief Private class that represents a table of tuples.
**
** This class represents a two-dimensional table of cells. Each cell is
** represented by a text string. Tuples are horizontal table entities
** identified by tuple indices, which are unsigned integers ranging from zero
** to the number of tuples minus one. Tuples are vertical table entities
** identified by tuple indices. The class provides methods for table
** construction and destruction, assignment operator, tuple and column based
** methods for addition, insertion, retrieval, update, deletion, cell based
** methods for update and retrieval and table printing.
*/
class TTable
{
  public:
    /**
    **  Constructs a tuple table.
    **
    **  \b Parameters: None
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table has 0 columns and 0 tuples.
    **
    **  \b Exceptions: None
    */
    TTable();

    /**
    **  Constructs a tuple table by copying from another table
    **    (copy constructor).
    **
    **  \param[in] inTable - reference to a table that will be copied to
    **    the newly constructed table
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table has the same content as the table
    **    referenced by \e inTable.
    **
    **  \b Exceptions: None
    */
    TTable(const TTable& inTable);

    /**
    **  Destructs a table.
    **
    **  \b Parameters: None
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \b Exceptions: None
    */
    virtual ~TTable();

    /**
    **  Copies a tuple table to another table (assignment operator).
    **
    **  \param[in] inTable - reference to the source table
    **
    **  \return Reference to the destination table
    **
    **  \pre None
    **
    **  \post Constructed table has the same content as the table
    **    referenced by \e inTable.
    **
    **  \b Exceptions: None
    */
    TTable& operator=(const TTable& inTable);

    /**
    **  Deletes all the content from the table.
    **
    **  \b Parameters: None
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post Table has 0 columns and 0 tuples.
    **
    **  \b Exceptions: None
    */
    void Clear();

    /**
    **  Retrieves the number of tuples in the table.
    **
    **  \b Parameters: None
    **
    **  \return The number of tuples in the table.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \b Exceptions: None
    */
    inline unsigned int GetNumTuples() const;

    /* \todo Re-use much of the comment of InsertTuple() */
    unsigned int AddTuple(const std::vector<std::string>& tuple =
      std::vector<std::string>());

    /**
    **  Inserts a new tuple at the specified tuple index and shifts,
    **    to the right by one, the specified existing tuple and all tuples
    **    after it.
    **
    **  \param[in] tupleIndex - the index of the tuple at which the new tuple
    **    is to be inserted
    **  \param[in] tuple - optional parameter that contains the values which
    **    are to be used to fill in the newly inserted tuple. If \e tuple is
    **    specified, filling starts at column index 0 and continues until size
    **    of \e tuple. If \e tuple is not specified, the newly inserted tuple is
    **    filled with empty values, where filling starts at column index 0 and
    **    ends at column index "number of columns - 1".
    **
    **  \return None
    **
    **  \pre \e tupleIndex must be greater than 0 and less than or equal to
    **    the number of tuples
    **
    **  \pre If \e tuple is specified, the size of \e tuple must be less than or
    **    equal to the number of columns.
    **  \pre The tuple which comes, in order, before the tuple with name
    **    \e atColName, must be non-empty. This is to prevent creation of
    **    non-rectangular tables.
    **
    **  \post If table is empty (0 columns) and \e tuple is specified,
    **    the number of columns is set to the size of \e tuple. Otherwise, the
    **    number of columns is unchanged.
    **
    **  \exception out_of_range - if \e tupleIndex is greater than than the
    **    number of tuples
    **  \exception out_of_range - if size of \e tuple is greater than the number
    **    of columns
    **  \exception out_of_range - if tuple, which comes, in order, before the
    **    tuple with name \e atColName, is empty.
    */
    void InsertTuple(const unsigned int tupleIndex,
      const std::vector<std::string>& tuple = std::vector<std::string>());

    void InsertTuple(const unsigned int tupleIndex,
      std::vector<std::string>::const_iterator tupleBeg,
      std::vector<std::string>::const_iterator tupleEnd);

    /**
    **  Inserts a new tuple at the specified tuple index and shifts,
    **    to the right by one, the specified existing tuple and all tuples
    **    after it.
    **
    **  \param[in] tupleIndex - the index of the tuple at which the new tuple
    **    is to be inserted
    **  \param[in] tuple - contains the values which are to be used to fill in
    **    the newly inserted tuple. Filling starts at column index 0 and
    **    continues until size of \e tuple.
    **
    **  \return None
    **
    **  \pre \e tupleIndex must be greater than 0 and less than the number
    **    of tuples
    **
    **  \pre The size of \e tuple must be less than or equal to the number
    **    of columns.
    **  \pre The tuple which comes, in order, before the tuple with name
    **    \e atColName, must be non-empty. This is to prevent creation of
    **    non-rectangular tables.
    **
    **  \post If table is empty (0 columns) and \e tuple is specified, the
    **    number of columns is set to the size of \e tuple. Otherwise, the
    **    number of columns is unchanged.
    **
    **  \exception out_of_range - if \e tupleIndex is greater than than the
    **    number of tuples
    **  \exception out_of_range - if size of \e tuple is greater than the number
    **    of columns
    **  \exception out_of_range - if tuple, which comes, in order, before the
    **    tuple with name \e atColName, is empty.
    */
    void FillTuple(const unsigned int tupleIndex,
      const std::vector<std::string>& tuple,
      const unsigned int fromColIndex = 0);

    void GetTuple(std::vector<std::string>& tuple,
      const unsigned int tupleIndex,
      const unsigned int fromColIndex, unsigned int toColIndex);

    const std::vector<std::string>& GetTuple(const unsigned int tupleIndex);

    /**
    **  Sets all cells in the tuple to empty string.
    **
    **  \param[in] colName - the name of the tuple
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Tuple with name \e colName must be present
    **
    **  \post Tuple length is unchanged.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if tuple with name \e colName
    **    does not exist
    */
    void ClearTuple(const unsigned int tupleIndex);

    /**
    **  Deletes a tuple from the table.
    **
    **  \param[in] colName - the name of the tuple
    **
    **  \return None
    **
    **  \pre \e colName must be non-empty
    **  \pre Tuple with name \e colName must be present
    **
    **  \post The number of table tuples is reduced by one.
    **
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if tuple with name \e colName
    **    does not exist
    */
    void DeleteTuple(const unsigned int tupleIndex);

    /**
    **  Retrieves the number of columns in the table.
    **
    **  \b Parameters: None
    **
    **  \return The number of columns in the table.
    **
    **  \pre None
    **
    **  \post None
    */
    inline unsigned int GetNumColumns() const;

    /**
    **  Adds a new column to the bottom end of the table. For an empty table,
    **  the number of inserted cells is equal to the number of table tuples.
    **  For a non-empty table, the number of inserted cells is equal to the
    **  number of non-empty tuples (this is in order to prevent creation of
    **  non-rectangular tables). The newly added column is, optionally,
    **  filled with values, starting at the first tuple.
    **
    **  \param[in] column - optional parameter that contains the values which
    **    are to be used to fill in the newly added column. Filling starts at
    **    the first tuple and continues until size of \e column.
    **
    **  \return The new number of columns after the column addition.
    **
    **  \pre Table must have at least one tuple, which can be empty.
    **  \pre If table is not empty and \e column is specified, the size of
    **    \e column must be less than or equal to the number of non-empty
    **    tuples. This is in order to prevent creation of non-rectangular
    **    tables.
    **  \pre If table is empty and \e column is specified, the size of
    **    \e column must be less than or equal to the number of tuples.
    **
    **  \post The number of columns is increased by one. 
    **
    **  \exception EmptyContainerException - if table has no tuples.
    **  \exception out_of_range - if table is not empty and size of \e column is
    **    greater than the number of non-empty tuples.
    **  \exception out_of_range - if table is empty and size of \e column is
    **    greater than the number of tuples.
    */
    unsigned int AddColumn(const std::vector<std::string>& col =
      std::vector<std::string>());

    /**
    **  Inserts a new column at the specified column index and shifts, down by
    **  one, the old column with the specified column index and all other
    **  columns below it.
    **  For an empty table, the number of inserted cells is equal to the
    **  number of table tuples. For a non-empty table, the number of
    **  inserted cells is equal to the number of non-empty tuples (this is
    **  in order to prevent creation of non-rectangular tables). The newly
    **  inserted column is optionally filled with values, starting at the
    **  first tuple.
    **
    **  \param[in] atColIndex - index of the column at which the new column is
    **    to be inserted. Note: If \e atColIndex is equal to the number of
    **    columns, the operation of this method is equivalent to AddRow().
    **  \param[in] column - optional parameter that contains the values which
    **    are to be used to fill in the newly inserted column. Filling starts at
    **    the first tuple and continues until size of \e column.
    **
    **  \return The new number of columns after the column insertion.
    **
    **  \pre Table must have at least one tuple, which can be empty.
    **  \pre \e atColIndex must be less than or equal to the number of table
    **    columns.
    **  \pre If table is not empty and \e column is specified, the size of
    **    \e column must be less than or equal to the number of non-empty
    **    tuples. This is in order to prevent creation of non-rectangular
    **    tables.
    **  \pre If table is empty and \e column is specified, the size of
    **    \e column must be less than or equal to the number of tuples.
    **
    **  \post The number of columns is increased by one. 
    **  \post Row indices of the columns which are below the inserted column are
    **    invalidated by being increased by 1.
    **
    **  \exception EmptyContainerException - if table has no tuples.
    **  \exception out_of_range - if \e atColIndex is greater than the number
    **    of table columns.
    **  \exception out_of_range - if table is not empty and size of \e column is
    **    greater than the number of non-empty tuples.
    **  \exception out_of_range - if table is empty and size of \e column is
    **    greater than the number of tuples.
    */
    unsigned int InsertColumn(const unsigned int atColIndex,
      const std::vector<std::string>& col = std::vector<std::string>());

    void InsertColumn(const unsigned int atColIndex,
      std::vector<std::string>::const_iterator colBeg,
      std::vector<std::string>::const_iterator colEnd);

    /**
    **  Fills, with values, a column at the specified column index, starting
    **  at the the first tuple.
    **
    **  \param[in] colIndex - index of the column that is to be filled.
    **  \param[in] column - values which are to be used to fill in the column.
    **    Filling starts at the first tuple and continues until size of
    **    \e column.
    **
    **  \return None
    **
    **  \pre \e colIndex must be greater than 0 and less than the number of
    **    table columns.
    **  \pre The size of \e column must be less than or equal to the number of
    **    non-empty tuples. This is in order to prevent creation of
    **    non-rectangular tables.
    **
    **  \post None 
    **
    **  \exception out_of_range - if \e colIndex is greater than or equal to
    **    the number of table columns.
    **  \exception out_of_range - if size of \e column is greater than the
    **    number of non-empty tuples.
    */
    void FillColumn(const unsigned int colIndex,
      const std::vector<std::string>& col,
      const unsigned int fromTupleIndex = 0);

    void FillColumn(const unsigned int colIndex,
      std::vector<std::string>::const_iterator colBeg,
      std::vector<std::string>::const_iterator colEnd,
      const unsigned int fromTupleIndex = 0);

    /**
    **  Retrieves the values in the specified column. 
    **
    **  \param[out] column - retrieved column values
    **  \param[in] colIndex - index of the column which values are to be
    **    retrieved.
    **  \param[in] fromColName - optional parameter which specifies the
    **    column location of the first cell to be retrieved. If not specified
    **    the first tuple cell is used.
    **  \param[in] toColName - optional parameter which specifies the
    **    column location of the last cell to be retrieved. If not specified
    **    the last non-empty-tuple cell is used.
    **
    **  \return None
    **
    **  \pre \e colIndex must be greater than 0 and less than the number of
    **    table columns.
    **  \pre If \e fromColName is specified, the tuple with name
    **    \e fromColName must be present and must be non-empty
    **  \pre If \e toColName is specified, the tuple with name
    **    \e toColName must be present and must be non-empty
    **  \pre If \e fromColName is different than \e toColName, it must come
    **    prior to it in the tuple order.
    **
    **  \post None
    **
    **  \exception out_of_range - if \e colIndex is less than 0 or greater
    **    than or equal to the number of table columns.
    **  \exception NotFoundException - If \e fromColName is specified and
    **    tuple with name \e fromColName does not exist
    **  \exception NotFoundException - If \e toColName is specified and
    **    tuple with name \e toColName does not exist
    **  \exception out_of_range - If \e fromColName is specified and
    **    tuple with name \e fromColName exists but is empty
    **  \exception out_of_range - If \e toColName is specified and
    **    tuple with name \e toColName exists but is empty
    **  \exception out_of_range - if \e fromColName is different than
    **    \e toColName and it comes after it in the tuple order.
    */
    void GetColumn(std::vector<std::string>& col, const unsigned int colIndex,
      const unsigned int fromTupleIndex, unsigned int toTupleIndex);

    /**
    **  Sets all cells in the column to empty string.
    **
    **  \param[in] colIndex - index of the column that is to be cleared.
    **
    **  \return None
    **
    **  \pre \e colIndex must be greater than 0 and less than the number of
    **    table columns.
    **
    **  \post None
    **
    **  \exception out_of_range - if \e colIndex is less than 0 or greater
    **    than or equal to the number of table columns.
    */
    void ClearColumn(const unsigned int colIndex);

    /**
    **  Deletes a column with the specified column index.
    **
    **  \param[in] colIndex - index of the column that is to be deleted.
    **
    **  \return None
    **
    **  \pre \e colIndex must be greater than 0 and less than the number of
    **    table columns.
    **
    **  \post Number of table columns is reduced by 1.
    **  \post Row indices of the columns which are below the deleted column are
    **    invalidated by being reduced by 1.
    **
    **  \exception out_of_range - if \e colIndex is less than 0 or greater
    **    than or equal to the number of table columns.
    */
    void DeleteColumn(const unsigned int colIndex);

    /**
    **  Updates a cell in the table.
    **
    **  \param[in] colIndex - column index of the cell that is to be updated.
    **  \param[in] colName - the name of the tuple
    **
    **  \return None
    **
    **  \pre \e colIndex must be greater than 0 and less than the number of
    **    table columns.
    **  \pre \e colName must be non-empty
    **  \pre Tuple with name \e colName must be present
    **
    **  \post None
    **
    **  \exception out_of_range - if \e colIndex is less than 0 or greater
    **    than or equal to the number of table columns.
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if tuple with name \e colName
    **    does not exist
    */
    std::string& operator()(const unsigned int tupleIndex,
      const unsigned int colIndex);

    /**
    **  Retrieves a reference to the cell in the table.
    **
    **  \param[in] colIndex - column index of the cell that is to be updated.
    **  \param[in] colName - the name of the tuple
    **
    **  \return Constant reference to a cell in the table.
    **
    **  \pre \e colIndex must be greater than 0 and less than the number of
    **    table columns.
    **  \pre \e colName must be non-empty
    **  \pre Tuple with name \e colName must be present
    **
    **  \post None
    **
    **  \exception out_of_range - if \e colIndex is less than 0 or greater
    **    than or equal to the number of table columns.
    **  \exception EmptyValueException - if \e colName is empty
    **  \exception NotFoundException - if tuple with name \e colName
    **    does not exist
    */
    const std::string& operator()(const unsigned int tupleIndex,
      const unsigned int colIndex) const;

    int Write(Serializer* ser, unsigned int& size);
    int Read(UInt32 index, Serializer* ser);

  private:
    unsigned int _numCols;

    std::vector<std::vector<std::string>*> _tuples;

    inline unsigned int IntColIndex(const unsigned int colIndex) const;

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    unsigned int _numDelCols;
    std::vector<unsigned int> _delColMap;

    void EnlargeColMap(const unsigned int numCols);
    void ReduceColMap(const unsigned int numCols);
    void MarkColDeleted(const unsigned int colIndex);
    void UnMarkColDeleted(const unsigned int colIndex);
#endif

};


std::ostream& operator<<(std::ostream& out, const TTable& sTable);


inline unsigned int TTable::GetNumTuples() const
{
    return(_tuples.size());
}


inline unsigned int TTable::GetNumColumns() const
{
#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    return(_numCols - _numDelCols);
#else
    return(_numCols);
#endif
}

inline unsigned int TTable::IntColIndex(const unsigned int colIndex) const
{

    // Returns the TTable internal column index
#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    return(_delColMap[colIndex]);
#else
    return(colIndex);
#endif
}

#endif // TTABLE_H
