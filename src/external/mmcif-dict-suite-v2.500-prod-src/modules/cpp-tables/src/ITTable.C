//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file ITTable.C
**
** \brief Implementation file for ITTable class.
*/


#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "Exceptions.h"
#include "GenString.h"
#include "ITTable.h"


using std::out_of_range;
using std::find;
using std::sort;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::make_pair;
using std::ostringstream;
using std::ostream;
using std::setw;
using std::endl;


const unsigned char ITTable::DEFAULT_OPTIONS = (ITTable::DT_STRING_VAL << 4);


ITTable::ITTable() : _ttable(), _orient(eROW_WISE) 
{

    Init();

}


ITTable::ITTable(eOrientation orient) : _ttable(), _orient(orient) 
{

    Init();

}


ITTable::ITTable(const ITTable& inTable)
{

    _ttable = inTable._ttable;

    _orient = inTable._orient;

    _compare_opts = inTable._compare_opts;

    _listsOfColumns = inTable._listsOfColumns;
    _unique = inTable._unique;

    _ser = inTable._ser;

    _indices = inTable._indices;

}


ITTable::~ITTable()
{

    Clear();

}


ITTable& ITTable::operator=(const ITTable& inTable)
{

    if (this != &inTable)
    {
        Clear();

        _ttable = inTable._ttable;

        _orient = inTable._orient;

        _compare_opts = inTable._compare_opts;

        _listsOfColumns = inTable._listsOfColumns;
        _unique = inTable._unique;

        _ser = inTable._ser;

        _indices = inTable._indices;
    }

    return(*this);

}


ITTable::eOrientation ITTable::GetOrientation()
{

    return(_orient);

}


void ITTable::AppendToColumn(const unsigned int colIndex, const string& cell)
{

        if (colIndex >= _compare_opts.size())
        {
            _compare_opts.insert(_compare_opts.end(),
                colIndex - _compare_opts.size() + 1, DEFAULT_OPTIONS);
        }

        // VLAD: INDEX CORRUPTION
        // This is not safe after index/indices is/are created
        // indices will be corrupted

        unsigned int numColumns;
        if (_orient == eCOLUMN_WISE)
            numColumns = _ttable.GetNumTuples();
        else
            numColumns = _ttable.GetNumColumns();

        if (colIndex < numColumns)
        {
            if (_orient == eCOLUMN_WISE)
            {
                _ttable.AddColumn();
                _ttable(colIndex, GetNumRows() - 1) = cell;
            }
            else
            {
                _ttable.AddTuple();
                _ttable(GetNumRows() - 1, colIndex) = cell;
            }

        }
        else
        {
            vector<string> tmpCol;
            tmpCol.push_back(cell);

            if (_orient == eCOLUMN_WISE)
                _ttable.InsertTuple(colIndex, tmpCol);
            else
            {
                _ttable.InsertColumn(colIndex, tmpCol);
            }
        }

        try
        {
            InsertEntry(GetNumRows() - 1);
        }
        catch (AlreadyExistsException& exc)
        {
            if (_orient == eCOLUMN_WISE)
                _ttable.DeleteTuple(colIndex);
            else
                _ttable.DeleteColumn(colIndex);

            throw;
        }

}


void ITTable::GetColumn(vector<string>& col, const unsigned int colIndex)
{

    col.clear();


    if (_orient == eCOLUMN_WISE)
        _ttable.GetTuple(col, colIndex, 0, _ttable.GetNumColumns());
    else
        _ttable.GetColumn(col, colIndex, 0, _ttable.GetNumTuples());

}


void ITTable::GetColumn(vector<string>& col, const unsigned int colIndex,
  const unsigned int fromRowIndex, unsigned int toRowIndex)
{

    col.clear();


    for (unsigned int rowI = fromRowIndex; rowI <= toRowIndex; ++rowI)
    {
        col.push_back(operator()(rowI, colIndex));
    }

}


void ITTable::GetColumn(vector<string>& col, const unsigned int colIndex,
  const vector<unsigned int>& rowIndex)
{


    col.clear();

    for (unsigned int index = 0; index < rowIndex.size(); ++index)
    {
        col.push_back(operator()(rowIndex[index], colIndex));
    }

}


void ITTable::ClearColumn(const unsigned int colIndex)
{

    vector<string> saveColumn;

    // VLAD - IMPROVE - if no indices present this is not needed !!!
    GetColumn(saveColumn, colIndex);

    if (_orient == eCOLUMN_WISE)
        _ttable.ClearTuple(colIndex);
    else
        _ttable.ClearColumn(colIndex);

    try
    {
        for (unsigned int rowI = 0; rowI < GetNumRows(); ++rowI)
            UpdateIndicesOnCellUpdate(rowI, colIndex);
    }
    catch (AlreadyExistsException& exc)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.FillTuple(colIndex, saveColumn);
        else
            _ttable.FillColumn(colIndex, saveColumn);

        throw;
    }

}


void ITTable::DeleteColumn(const unsigned int colIndex)
{

    if (_orient == eCOLUMN_WISE)
        _ttable.DeleteTuple(colIndex);
    else
        _ttable.DeleteColumn(colIndex);

    _compare_opts.erase(_compare_opts.begin() + colIndex);

    // Deletes all indices build on the column colIndex
    unsigned int numIndices = _indices.size();

    for (unsigned int i = 0; i < numIndices; i++)
    {
        if (IsColumnInIndex(i, colIndex))
        {
            DeleteIndex(i);
        }
    }

    UpdateColListOnColDelete(colIndex);

}


unsigned int ITTable::AddRow(const vector<string>& row)
{

    return(InsertRow(GetNumRows(), row));

}


unsigned int ITTable::InsertRow(const unsigned int rowIndex,
  const vector<string>& row)
{

    if ((!row.empty()) && (row.size() > GetNumColumns()))
    {
        // Wrong row index.
        throw out_of_range("Invalid row size in ITTable::InsertRow");
    }

    unsigned int oldNumRows = GetNumRows();

    if (oldNumRows == 0)
    {
        for (unsigned int colI = 0; colI < GetNumColumns(); ++colI)
        {
            if (colI < row.size())
                AppendToColumn(colI, row[colI]);
            else
                AppendToColumn(colI, string());
        }
    }
    else
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.InsertColumn(rowIndex, row);
        else
            _ttable.InsertTuple(rowIndex, row);

        // VLAD-PERFORMANCE. THIS CAN BE IMPROVED BY ONLY UPDATING
        // THE ROWS BELOW THE INSERTED ONE !!!
        try
        {
            InsertEntry(rowIndex);
        }
        catch (AlreadyExistsException& exc)
        {
            if (_orient == eCOLUMN_WISE)
                _ttable.DeleteColumn(rowIndex);
            else
                _ttable.DeleteTuple(rowIndex);

            throw;
        }
    }

    return(GetNumRows());

}

void ITTable::FillRow(const unsigned int rowIndex, const vector<string>& row)
{

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::FillRow");
    }

    vector<string> oldRow;
    GetRow(oldRow, rowIndex, 0, GetNumColumns() - 1);

    if (_orient == eCOLUMN_WISE)
        _ttable.FillColumn(rowIndex, row);
    else
        _ttable.FillTuple(rowIndex, row);

    try
    {
        UpdateIndices(rowIndex);
    }
    catch (AlreadyExistsException& exc)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.FillColumn(rowIndex, oldRow);
        else
            _ttable.FillTuple(rowIndex, oldRow);

        UpdateIndices(rowIndex);

        throw;
    }

}


void ITTable::GetRow(vector<string>& row, const unsigned int rowIndex,
  const unsigned int fromColIndex, unsigned int toColIndex)
{

    row.clear();

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::GetRow");
    }

    if (_orient == eCOLUMN_WISE)
        _ttable.GetColumn(row, rowIndex, fromColIndex, toColIndex + 1);
    else
        _ttable.GetTuple(row, rowIndex, fromColIndex, toColIndex + 1);

}


const vector<string>& ITTable::GetRow(const unsigned int rowIndex)
{

    if (_orient != eROW_WISE)
    {
        // Wrong table orientation.
        throw InvalidStateException("Cannot get row reference on "\
          "column-wise oriented table", "ITTable::GetRow");
    }

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::GetRow");
    }

    return(_ttable.GetTuple(rowIndex));

}


void ITTable::ClearRow(const unsigned int rowIndex)
{

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::ClearRow");
    }

    vector<string> oldRow;
    GetRow(oldRow, rowIndex, 0, GetNumColumns() - 1);

    if (_orient == eCOLUMN_WISE)
        _ttable.ClearColumn(rowIndex);
    else
        _ttable.ClearTuple(rowIndex);

    try
    {
        UpdateIndices(rowIndex);
    }
    catch (AlreadyExistsException& exc)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.FillColumn(rowIndex, oldRow);
        else
            _ttable.FillTuple(rowIndex, oldRow);

        UpdateIndices(rowIndex);

        throw;
    }

}


void ITTable::DeleteRow(const unsigned int rowIndex)
{

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::DeleteRow");
    }

    if (_orient == eCOLUMN_WISE)
        _ttable.DeleteColumn(rowIndex);
    else
        _ttable.DeleteTuple(rowIndex);

    DeleteEntry(rowIndex);

}


void ITTable::DeleteRows(const vector<unsigned int>& rows)
{

    // VLAD ERROR CHECKING ALL THE INDICES in rows
    // VLAD IMPROVE WITH REVERSE ITERATORS or more meaningfull index loops

    for (unsigned int rowIndex = 0; rowIndex < rows.size(); ++rowIndex)
    {
        DeleteRow(rows[rows.size() - 1 - rowIndex]);
    }

}


void ITTable::FindDuplicateRows(const vector<unsigned int>& colIndices,
  vector<pair<unsigned int, unsigned int> >& duplRows, const bool keep,
  const eSearchDir searchDir)
{

    duplRows.clear();

    // VLAD ERROR HANDLING - CHECK INDICES
    // for keep = 1, reports duplicate rows, but keeps them in table
    // for keep = 0, reports duplicate rows and deletes them from table
    // searchDir can be forward or backward

    unsigned int realRowIndex;
    tIndex::iterator pos;
    string value;

    Char::eCompareType compareType = GetCompareType(colIndices);

    tIndex tmpindex(compareType);

    for (unsigned int rowIndex = 0; rowIndex < GetNumRows(); rowIndex++)
    {
        if (searchDir == eFORWARD)
        {
            realRowIndex = rowIndex;
        }
        else
        {
            realRowIndex = GetNumRows() - rowIndex - 1;
        }

        value = AggregateRow(colIndices, realRowIndex);

        pos = tmpindex.find(value);
        if (pos != tmpindex.end())
        {
            duplRows.push_back(make_pair(realRowIndex, pos->second));
        }
        else
        {
            tIndex::value_type valType(value, realRowIndex);
            tmpindex.insert(valType);
        }

        value.clear();
    }

    if (!keep)
    {
        vector<unsigned int> duplRowIndices;
        for (unsigned int rowI = 0; rowI < duplRows.size(); rowI++)
        {
            duplRowIndices.push_back(duplRows[rowI].first);
        }

        DeleteRows(duplRowIndices);
    }

}


bool ITTable::AreListsOfColumnsValid(const vector<unsigned int>& colIndices)
{

    for (unsigned int index = 0; index < colIndices.size(); index++)
    {
        if (colIndices[index] >= GetNumColumns())
        {
            return(false);
        }
    }

    return(true);

}


void ITTable::CreateIndex(const vector<unsigned int>& colIndices,
  const unsigned int unique)
{

    if (colIndices.empty())
    {
        throw EmptyValueException("Empty column indices",
          "ITTable::CreateIndex");
    }

    VerifyColumnsIndices(colIndices);

    Char::eCompareType compareType = GetCompareType(colIndices);

    // Create one or another type of multimap with/without
    // case-insensitivity
    tIndex index(compareType);

    _listsOfColumns.push_back(colIndices);
    _indices.push_back(index);
    _unique.push_back(unique);

    unsigned int lastIndex = _indices.size() - 1;

    for (unsigned int i = 0; i < GetNumRows(); i++)
        UpdateIndex(lastIndex, i);

}


void ITTable::InsertEntry(const unsigned int rowIndex)
{

    unsigned int numIndices = _indices.size();

    for (unsigned int i = 0; i < numIndices; i++)
    {
        if (_unique[i] == 1)
        {
            // Found key index. Cannot change the table. Throw exception.
            throw AlreadyExistsException("Attempting to change the table "\
              "that has a key search index", "ITTable::InsertEntry");
        }
    }

    for (unsigned int i = 0; i < numIndices; i++)
    {
        if (_unique[i] == 0)
            InsertIndexEntry(i, rowIndex);
    }

}


void ITTable::DeleteEntry(const unsigned int rowIndex)
{

    unsigned int numIndices = _indices.size();

    for (unsigned int i = 0; i < numIndices; i++)
    {
        DeleteIndexEntry(i, rowIndex);
    }

}


void ITTable::InsertIndexEntry(const unsigned int indexIndex,
  const unsigned int rowIndex)
{

#ifdef OLD_DESCRIPTION_FIX_IT
    // indexIndex is the index of an index
    // rowIndex is the row index
    // This function preforms the following:
    //   if unique flag is set to 0, table rows with duplicate values
    //     are allowed
    //   if unique flag is set to 1, table must not have duplicate values
    //     and rows with repeated values are deleted.
    //
    // More information on unique flag being set to 1. If a row with rowIndex
    // already exists and the value of the new row with rowIndex already
    // exists, delete the row with rowIndex from the table and from the index.
    // VLAD: SIMPLIFY THE ABOVE unique == 1 processing.
#endif

    // New row has been inserted at index rowIndex
    // See if new values in a row are unique for unique = 1. If not throw
    // exception and exit.

    // If entry can go in, find the old rowIndex and insert new values
    // Find all indices with larger row index than rowIndex and increase
    // the value by one


    // VLAD ERROR HANDLING
    if (rowIndex > GetNumRows() - 1)
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::InsertIndexEntry");
    }

    if (!AreListsOfColumnsValid(_listsOfColumns[indexIndex]))
    {
        // VLAD - ERROR HANDLING - A LOT OF PLACES ATTEMPT TO CREATE INDICES
        // ON EMPTY COLUMNS
        return;
    }

    string value = SubRowValue(_listsOfColumns[indexIndex], rowIndex);

    // See if new values in a row are unique for unique = 1. If not throw
    // exception and exit.
    if (_unique[indexIndex] == 1)
    {
        pair<tIndex::iterator, tIndex::iterator> range;

        range = _indices[indexIndex].equal_range(value);

        if (range.first != range.second)
        {
            // Found a value. Values must be unique. Throw exception.
            throw AlreadyExistsException("Duplicate value on key index",
              "ITTable::InsertIndexEntry");
        }
    }

    // Find all indices with larger row index than rowIndex and increase
    // the value by one
    // VLAD IMPROVE. This and below for loop can be made into a method
    // IndexSearchForRow()
    for (tIndex::iterator pos = _indices[indexIndex].begin();
      pos != _indices[indexIndex].end(); ++pos)
    {
        if (pos->second >= rowIndex)
        {
            ++(pos->second);
        }
    }

    tIndex::value_type valType(value, rowIndex);

    _indices[indexIndex].insert(valType);

}


void ITTable::DeleteIndexEntry(const unsigned int indexIndex,
  const unsigned int rowIndex)
{

#ifdef OLD_DESCRIPTION_FIX_IT
    // indexIndex is the index of an index
    // rowIndex is the row index
    // This function preforms the following:
    //   if unique flag is set to 0, table rows with duplicate values
    //     are allowed
    //   if unique flag is set to 1, table must not have duplicate values
    //     and rows with repeated values are deleted.
    //
    // More information on unique flag being set to 1. If a row with rowIndex
    // already exists and the value of the new row with rowIndex already
    // exists, delete the row with rowIndex from the table and from the index.
    // VLAD: SIMPLIFY THE ABOVE unique == 1 processing.
#endif

    // Row has been deleted at index rowIndex

    // Delete entry with that row index
    // Find all indices with larger row index than rowIndex and decrease
    // the value by one

    // VLAD ERROR HANDLING
    if (rowIndex >= GetNumRows() + 1)
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::InsertIndexEntry");
    }

    if (!AreListsOfColumnsValid(_listsOfColumns[indexIndex]))
    {
        // VLAD - ERROR HANDLING - A LOT OF PLACES ATTEMPT TO CREATE INDICES
        // ON EMPTY COLUMNS
        return;
    }

    for (tIndex::iterator pos = _indices[indexIndex].begin();
      pos != _indices[indexIndex].end(); ++pos)
    {
        if (pos->second == rowIndex)
        {
            _indices[indexIndex].erase(pos);
            break;
        }
    }

    // Find all indices with larger row index than rowIndex and increase
    // the value by one
    // VLAD IMPROVE. This and below for loop can be made into a method
    // IndexSearchForRow()
    for (tIndex::iterator pos = _indices[indexIndex].begin();
      pos != _indices[indexIndex].end(); ++pos)
    {
        if (pos->second >= rowIndex)
        {
            --(pos->second);
        }
    }

}


void ITTable::UpdateIndex(const unsigned int indexIndex,
  const unsigned int rowIndex)
{

    // indexIndex is the index of an index
    // rowIndex is the row index
    // This function preforms the following:
    //   if unique flag is set to 0, table rows with duplicate values
    //     are allowed
    //   if unique flag is set to 1, table must not have duplicate values
    //     and rows with repeated values are deleted.
    //
    // More information on unique flag being set to 1. If a row with rowIndex
    // already exists and the value of the new row with rowIndex already
    // exists, delete the row with rowIndex from the table and from the index.
    // VLAD: SIMPLIFY THE ABOVE unique == 1 processing.

    // VLAD ERROR HANDLING
    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::UpdateIndex");
    }

    if (!AreListsOfColumnsValid(_listsOfColumns[indexIndex]))
    {
        // VLAD - ERROR HANDLING - A LOT OF PLACES ATTEMPT TO CREATE INDICES
        // ON EMPTY COLUMNS
        return;
    }

    string value = SubRowValue(_listsOfColumns[indexIndex], rowIndex);

    // See if new values in a row are unique for unique = 1. If not throw
    // exception and exit.
    if (_unique[indexIndex] == 1)
    {
        pair<tIndex::iterator, tIndex::iterator> range;

        range = _indices[indexIndex].equal_range(value);

        if (range.first != range.second)
        {
            // Found a value. Values must be unique. Throw exception.
            throw AlreadyExistsException("Duplicate value on key index",
              "ITTable::InsertIndexEntry");
        }
    }

    // VLAD IMPROVE. This and below for loop can be made into a method
    // IndexSearchForRow()

    // See if a pair with the row having rowIndex already exists. If it
    // does delete it first.
    for (tIndex::iterator pos = _indices[indexIndex].begin();
      pos != _indices[indexIndex].end(); ++pos)
    {
        if (pos->second == rowIndex)
        {
            _indices[indexIndex].erase(pos);
            break;
        }
    }

    // Duplicate values are allowed, just store value/index pair.
    tIndex::value_type valType(value, rowIndex);
    _indices[indexIndex].insert(valType);

}


void ITTable::ClearIndex(const unsigned int indexIndex)
{

    _listsOfColumns[indexIndex].clear();
    _indices[indexIndex].clear();

}


void ITTable::DeleteIndex(const unsigned int indexIndex)
{

    // VLAD: THIS MAY NOT BE NEEDED. DEPENDS ON HOW STL CONTAINERS
    // FREE RESOURCES IN ERASE CALLS
    ClearIndex(indexIndex);

    _listsOfColumns.erase(_listsOfColumns.begin() + indexIndex);
    _indices.erase(_indices.begin() + indexIndex);
    _unique.erase(_unique.begin() + indexIndex);

}


int ITTable::FindIndex(const vector<unsigned int>& colIndices)
{

  int equal;

  if (colIndices.empty())
  {
      return(-1);
  }

  if (_indices.empty())
  {
      return(-1);
  }

  Char::eCompareType compareType = GetCompareType(colIndices);

  unsigned int numIndices = _indices.size();

  for (unsigned int i = 0; i < numIndices; i++)
  {
      if (_listsOfColumns[i].size() != colIndices.size())
      {
          continue;
      }

      // VLAD INDEX SIMPLIFY CompareColIndices()
      // Assume that all columns are equal
      equal = 0;
      for (unsigned int j = 0; j < _listsOfColumns[i].size(); j++)
      {
          if (_listsOfColumns[i][j] != colIndices[j])
          {
              // Found one non-equal column. Set the flag and brake the loop.
              equal = 1;
              break;
          }
      }

      if (equal == 0)
      {
          if ((_indices[i].key_comp()).GetCompareType() == compareType)
          {
              // The index is matched only if it also has the same
              // compare type
              return(i);
          }
      }
  }

  return(-1);

}


void ITTable::UpdateIndices(const unsigned int rowIndex)
{

    unsigned int numIndices = _indices.size();

    for (unsigned int i = 0; i < numIndices; i++)
    {
        if (_unique[i] == 1)
            // Found key index. Cannot change the table. Throw exception.
            throw AlreadyExistsException("Attempting to change the table "\
              "that has a key search index", "ITTable::UpdateIndices");
    }

    for (unsigned int i = 0; i < numIndices; i++)
    {
        if (_unique[i] == 0)
            UpdateIndex(i, rowIndex);
    }

}


void ITTable::ClearIndices()
{

    for(unsigned int i = 0; i < _listsOfColumns.size(); i++)
    {
        ClearIndex(i);
    }

}


bool ITTable::IsColumnInIndex(const unsigned int indexIndex,
  const unsigned int colIndex)
{

    vector<unsigned int>::iterator pos =
      find(_listsOfColumns[indexIndex].begin(),
      _listsOfColumns[indexIndex].end(), colIndex);

    if (pos != _listsOfColumns[indexIndex].end())
    {
        return(true);
    }
    else
    {
        return(false);
    }

}


string ITTable::CellValue(const unsigned int colIndex,
  const unsigned int rowIndex)
{

    return(ConvertString(operator()(rowIndex, colIndex), colIndex));

}


string ITTable::ConvertString(const string& value,
  const unsigned int colIndex)
{
  string ret;
  unsigned char temp = _compare_opts[colIndex];
  switch (( temp & DT_MASK) >> 4) {
  case DT_INTEGER_VAL: {
      ConvertToInt(value, ret);
      break;
  }
#ifdef VLAD_DELETED
  case DT_DOUBLE_VAL: {
      ConvertDouble(value, ret);
      break;
  }
#endif
  case DT_STRING_VAL: {
      if ( temp & SC_MASK ) {
        if ( temp & WS_MASK ) {
          ConvertToLowerNoWhiteSpace(value, ret);
	  break;
	}
        else {
          String::LowerCase(value, ret);
	  break;
	}
      } else if (temp & WS_MASK) {
	String::RemoveWhiteSpace(value, ret);
	break;
      } else {
	ret=value;
	break;
      }
  }
  default: {
      ret=value;
      break;
  }
  }
  return ret;
}


string ITTable::MultiStringsValue(const vector<string>& values,
    const vector<unsigned int>& colIndices)
{

    string value;

    for (unsigned int index = 0; index < values.size(); index++)
    {
        AppendToAndDelimit(value,
          ConvertString(values[index], colIndices[index]));
    }

    return(value);

}


string ITTable::SubRowValue(const vector<unsigned int>& colIndices,
  const unsigned int rowIndex)
{

    string value;

    for (unsigned int index = 0; index < colIndices.size(); index++)
    {
        AppendToAndDelimit(value, CellValue(colIndices[index], rowIndex));
    }

    return(value);

}


string ITTable::AggregateRow(const vector<unsigned int>& colIndices,
  const unsigned int rowIndex)
{

    string value;

    for (unsigned int index = 0; index < colIndices.size(); index++)
    {
        AppendToAndDelimit(value, operator()(rowIndex, colIndices[index]));
    }

    return(value);

}


void ITTable::Init()
{

    _ser = NULL;

}


void ITTable::Clear()
{

    _ttable.Clear();

    _listsOfColumns.clear();
    _indices.clear();
    _unique.clear();
    _compare_opts.clear();

}


void ITTable::ValidateOptions(unsigned int colIndex)
{

    if (colIndex >= _compare_opts.size())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ITTable::ValidateOptions");
    }

    unsigned char a = 0, b = 0;

    b = _compare_opts[colIndex];

    a = b & DT_MASK;

    if (((a >> 4) > LAST_DT_VALUE) || (a == 0))
    {
        a = DEFAULT_OPTIONS;
    }
    else
        a = b;

    a &= ~(3 << 2); // not necessary, but clean up 3rd and 4th bit

    _compare_opts[colIndex] = a;

}


string ITTable::CreateInternalIndexName(const unsigned int indexIndex)
{

    // VLAD HARCODED CONST
    string name("index_");

    name += String::IntToString(indexIndex);

    return(name);

}


void ITTable::Search(vector<unsigned int>& res, const string& target,
  const unsigned int colIndex, const eSearchType searchType)
{

    res.clear();

    vector<string> targets;
    targets.push_back(target);

    vector<unsigned int> colIds;

    colIds.push_back(colIndex);

    Search(res, targets, colIds, searchType);

}


void ITTable::Search(vector<unsigned int>& res, const vector<string>& targets,
  const vector<unsigned int>& colIds, const unsigned int indexNum,
  const eSearchType searchType)
{

    res.clear();

    unsigned int len = targets.size();
    if ((len > colIds.size()) || (len < 1))
    {
        // To index less than from index.
        throw out_of_range("Invalid search target range in ITTable::Search");
    }

    VerifyColumnsIndices(colIds);

    string value = MultiStringsValue(targets, colIds);

    /*
    ** This code searches the container for all key/value pairs, where key
    ** equals "value" and stores all the values in the ret.
    */

    pair<tIndex::iterator, tIndex::iterator> range;

    range = _indices[indexNum].equal_range(value);

    tIndex::iterator start, afterLast;

    switch (searchType)
    {
        case eEQUAL:
           start = range.first;
           afterLast = range.second;
           break;
        case eLESS_THAN:
           start = _indices[indexNum].begin();
           afterLast = range.first;
           break;
        case eLESS_THAN_OR_EQUAL:
           start = _indices[indexNum].begin();
           afterLast = range.second;
           break;
        case eGREATER_THAN:
           start = range.second;
           afterLast = _indices[indexNum].end();
           break;
        case eGREATER_THAN_OR_EQUAL:
           start = range.first;
           afterLast = _indices[indexNum].end();
           break;
        default:
           // VLAD ERROR CHECKING: DO THIS CHECK PRIOR TO SEARCH
           return;
           break;
    }

    for (tIndex::iterator pos = start; pos != afterLast; ++pos)
    {
        res.push_back(pos->second);
    }

    if (searchType == eEQUAL)
        sort(res.begin(), res.end());

}


unsigned int ITTable::FindFirst(const vector<string>& targets,
  const vector<unsigned int>& colIds, const unsigned int indexNum)
{

    unsigned int ret = GetNumRows();

    if (colIds.empty())
    {
        throw EmptyValueException("Empty column indices",
          "ITTable::FindFirst");
    }

    unsigned int len = targets.size();
    if ((len < colIds.size()) || (len < 1))
    {
        // To index less than from index.
        throw out_of_range("Invalid search target range in ITTable::FindFirst");
    }

    VerifyColumnsIndices(colIds);

    string value = MultiStringsValue(targets, colIds);

    // Assume that it is not found

    tIndex::iterator pos = _indices[indexNum].find(value);
    if (pos != _indices[indexNum].end())
    {
        // Found
        ret = pos->second;
    }

    return(ret);

}


ostream& operator<<(ostream& out, const ITTable& itTable)
{

    unsigned int numRows = itTable.GetNumRows();

    for (unsigned int j = 0; j < numRows; j++)
    {
        out << setw(5) << j;

        for (unsigned int i = 0; i < itTable.GetNumColumns(); i++)
        {
            // VLAD HARDCODED CONST
            out << setw(10) << itTable(j, i) << " ";
        }

        out << endl;
    }

    return(out);

}


void ITTable::SetSerializer(Serializer* ser)
{

    _ser = ser;

}


int ITTable::WriteObject(Serializer* ser, int& size)
{

    SInt32 ret = 0;

    unsigned int unused = 0;
    _ttable.Write(ser, unused);

    unused = ser->WriteUInt32(_orient);

    return ret;

}


int ITTable::Read(unsigned int indexInFile, Serializer* ser)
{
    return(GetObject(indexInFile, ser));
}


int ITTable::Write(Serializer* ser, int& size)
{
    return(WriteObject(ser, size));
}


int ITTable::GetObject(UInt32 index, Serializer* ser)
{

    index = _ttable.Read(index, ser);

    _orient = (eOrientation)ser->ReadUInt32(index);
    index++;

    if (_compare_opts.size() != GetNumColumns())
    {
        _compare_opts.insert(_compare_opts.end(),
          GetNumColumns() - _compare_opts.size(), DEFAULT_OPTIONS);
    }

    return(index);

}


void ITTable::ConvertToInt(const string& a, string& ret)
{

  ret.clear();

  const unsigned int len = a.size();

  if (len > INT_LIMIT)
  {
    return;
  }

  int intValue = 0;

  if (len != 0)
  {
      intValue = String::StringToInt(a);
  }

  if (intValue < 0)
  {
      const string maxNum(len - 1, '9');

      ret.push_back('-');
      ret.append(INT_LIMIT - len, '9');

      // maxNum - (-intValue)
      ret += String::IntToString(String::StringToInt(maxNum) + intValue);
  }
  else
  {
      ret.assign(INT_LIMIT - len, '0');
      ret += a;
  }

}


void ITTable::ConvertDouble(const string& a, string& ret){
/* This method convert string representing a double value into
   another string. when we wont to compare two double value
   we can compare two string value (converted value) and
   result will be the same
*/

  const int MAX_DECIMAL_DIGIT = 9;
  string m;  // Mantissa before decimal point
  string r;  // Mantissa after decimal point (remainder)
  string e;  // Exponent, without letter "e"
  string aa; // Exponential representation of the input double
  // Maximum remainder complement value

  int eint;  // Exponent value
  double adouble; // Input value


  // Convert the original double representation to mantissa/exponent
  // representation.
  adouble = String::StringToDouble(a);

  ostringstream outStringStream;
  outStringStream << std::scientific << std::showpoint <<
    std::setprecision(MAX_PRECISION) << adouble;

  aa = outStringStream.str();

  // Extract mantissa, remainder and exponent
  string::size_type dotCharIndex = aa.find('.');
  if (dotCharIndex == string::npos)
  {
    return;
  }

  m.append(aa, 0, dotCharIndex);

  string::size_type expCharIndex = aa.find('e');
  if (expCharIndex == string::npos)
  {
    // VLAD BUG: Fix this
    return;
  }

  r.append(aa, dotCharIndex + 1, expCharIndex - dotCharIndex - 1);
  e.append(aa, expCharIndex + 1, aa.size() - expCharIndex - 1);

  eint = String::StringToInt(e);

  if (adouble < 0)
  {
    ret.push_back('-');
  }
  else
  {
    ret.push_back('0');
  }

  if (eint < 0)
  {
    const string maxExpNum(e.size() - 1, '9');

    ret.push_back('-');

    // Complement the exponent
    ret.append(EXPONENT - e.size(), '9');

    // maxExpNum - (-eint)
    ret += String::IntToString(String::StringToInt(maxExpNum) + eint);
  }
  else
  {
    ret.push_back('0');
    ret.append(EXPONENT - e.size(), '0');
    ret.append(e, 1, e.size() - 1);
  }

  if (adouble < 0)
  {
    // Complement mantissa
    // MAX_DECIMAL_DIGIT - (-mantissa before dot)
    ret += String::IntToString(MAX_DECIMAL_DIGIT + String::StringToInt(m));

    ret.push_back('.');

    // Complement each digit in the remainer
    for (unsigned int index = 0; index < r.size(); ++index)
    {
      ret += String::IntToString(MAX_DECIMAL_DIGIT -
        String::StringToInt(r.substr(index, 1)));
    }
  }
  else
  {
    ret += m;
    ret += '.';
    ret += r;
  }

}


void ITTable::ConvertToLowerNoWhiteSpace(const string& a, string& ret)
{
    ret.clear();

    String::RemoveWhiteSpace(a, ret);
    String::LowerCase(ret);
}


Char::eCompareType
  ITTable::GetCompareType(const vector<unsigned int>& colIndices)
{

    // Initially set it to case sensitive type.
    Char::eCompareType compareType = Char::eCASE_SENSITIVE;

    // In setting the compare type, only the options of the column
    // at index 0 are considered. It is assumed that all the rest of
    // the columns have the same options as the first column.

    unsigned char opts = _compare_opts[colIndices[0]];
    switch ((opts & DT_MASK) >> 4)
    {
        case DT_STRING_VAL:
            if (opts & SC_MASK)
            {
                // Case in-sensitive column.
                compareType = Char::eCASE_INSENSITIVE;
            }
            break;
        default:
            break;
    }

    return(compareType);

}


void ITTable::UpdateColListOnColInsert(const unsigned int colIndex)
{

    // Updating lists of columns. If a column in the list has column
    // index >= colIndex then increment it by 1
    for (unsigned int i = 0; i < _listsOfColumns.size(); i++)
    {
        for (unsigned int j = 0; j < _listsOfColumns[i].size(); j++)
        {
            if (_listsOfColumns[i][j] >= colIndex)
                _listsOfColumns[i][j]++;
        }
    }

}


void ITTable::UpdateColListOnColDelete(const unsigned int colIndex)
{

    // Updating lists of columns. If a column in the list has column
    // index >= colIndex then decrement it by 1
    for (unsigned int i = 0; i < _listsOfColumns.size(); i++)
    {
        for (unsigned int j = 0; j < _listsOfColumns[i].size(); j++)
        {
            if (_listsOfColumns[i][j] >= colIndex)
                _listsOfColumns[i][j]--;
        }
    }

}


void ITTable::UpdateIndicesOnCellUpdate(const unsigned int rowIndex,
  const unsigned int colIndex)
{

    for (unsigned int i = 0; i < _indices.size(); i++)
    {
        if (_unique[i] == 1)
        {
            // Found key index. Cannot change the table. Throw exception.
            throw AlreadyExistsException("Attempting to change the table "\
              "that has a key search index",
              "ITTable::UpdateIndicesOnCellUpdate");
        }
    }

    for (unsigned int i = 0; i < _indices.size(); i++)
    {
        if (_unique[i] == 0)
        {
            for (unsigned int k = 0; k < _listsOfColumns[i].size(); k++)
            {
                if (_listsOfColumns[i][k] == colIndex)
                {
                    UpdateIndex(i, rowIndex);
                    break;
                }
            }
        }
    }

}


void ITTable::FillColumn(const vector<string>& col,
  const unsigned int colIndex)
{

    if (col.empty())
    {
        return;
    }

    FillColumn(colIndex, col.begin(), col.end());

}


void ITTable::FillColumn(const unsigned int colIndex,
  vector<string>::const_iterator colBeg, vector<string>::const_iterator colEnd)
{

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ITTable::FillColumn");
    }

    if (colBeg > colEnd)
    {
        // Wrong tuple range.
        throw out_of_range("Invalid column range in ITTable::FillColumn");
    }

    if (colBeg == colEnd)
    {
        return;    
    }

    vector<string> saveColumn;

    // VLAD - IMPROVE - if no indices present this is not needed !!!
    GetColumn(saveColumn, colIndex);

    unsigned int rowI = 0;
    for (vector<string>::const_iterator currIter = colBeg; currIter < colEnd;
      ++currIter, ++rowI)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable(colIndex, rowI) = *currIter;
        else
            _ttable(rowI, colIndex) = *currIter;
    }

    try
    {
        // VLAD - IMPROVE - NO NEED TO GO till GetNumRows(). It is sufficient
        // to go until colEnd - colBegin!!
        for (unsigned int rowI = 0; rowI < GetNumRows(); ++rowI)
            UpdateIndicesOnCellUpdate(rowI, colIndex);
    }
    catch (AlreadyExistsException& exc)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.FillTuple(colIndex, saveColumn);
        else
            _ttable.FillColumn(colIndex, saveColumn);

        throw;
    }

}


void ITTable::AppendToColumn(const unsigned int colIndex,
  const vector<string>& col)
{

    // VLAD: INDEX CORRUPTION
    // This is not safe after index/indices is/are created
    // indices will be corrupted

    if (colIndex >= _compare_opts.size())
    {
        _compare_opts.insert(_compare_opts.end(),
            colIndex - _compare_opts.size() + 1, DEFAULT_OPTIONS);
    }

    unsigned int numColumns;

    if (_orient == eCOLUMN_WISE)
        numColumns = _ttable.GetNumTuples();
    else
        numColumns = _ttable.GetNumColumns();

    if (colIndex < numColumns)
    {
        unsigned int oldNumRows = GetNumRows();

        for (unsigned int colI = 0; colI < col.size(); ++colI)
        {
            if (_orient == eCOLUMN_WISE)
                _ttable.AddColumn();
            else
                _ttable.AddTuple();
        }

        for (unsigned int colI = 0; colI < col.size(); ++colI)
        {
            if (_orient == eCOLUMN_WISE)
                _ttable(colIndex, oldNumRows - 1 + colI) = col[colI];
            else
                _ttable(oldNumRows - 1 + colI, colIndex) = col[colI];
        }
    }
    else
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.InsertTuple(colIndex, col);
        else
            _ttable.InsertColumn(colIndex, col);
    }

    for (unsigned int rowI = GetNumRows() - col.size(); rowI < GetNumRows();
      ++rowI)
    {
        try
        {
            InsertEntry(rowI);
        }
        catch (AlreadyExistsException& exc)
        {
            for (unsigned int delRowI = GetNumRows() - col.size(); delRowI <
              GetNumRows(); ++delRowI)
            {
                if (_orient == eCOLUMN_WISE)
                    _ttable.DeleteColumn(delRowI);
                else
                    _ttable.DeleteTuple(delRowI);
            }

            for (unsigned int delRowI = GetNumRows() - col.size(); delRowI <
              rowI; ++delRowI)
            {
                DeleteEntry(delRowI);
            }

            throw;
        }
    }

}


int ITTable::UpdateCell(const string& cell, const unsigned int colIndex,
  const unsigned int rowIndex)
{

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::UpdateCell");
    }

    string saveCell = operator()(rowIndex, colIndex);

    if (_orient == eCOLUMN_WISE)
        _ttable(colIndex, rowIndex) = cell;
    else
        _ttable(rowIndex, colIndex) = cell;

    try
    {
        UpdateIndicesOnCellUpdate(rowIndex, colIndex);
    }
    catch (AlreadyExistsException& exc)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable(colIndex, rowIndex) = saveCell;
        else
            _ttable(rowIndex, colIndex) = saveCell;

        throw;
    }

    return NO_TABLE_ERROR;

}


const string& ITTable::operator()(const unsigned int rowIndex,
  const unsigned int colIndex) const
{

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ITTable::operator()");
    }

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ITTable::operator()");
    }

    if (_orient == eCOLUMN_WISE)
        return(_ttable(colIndex, rowIndex));
    else
        return(_ttable(rowIndex, colIndex));

}


int ITTable::SetFlags(const unsigned char newOpts, const unsigned int colIndex)
{

    if (colIndex >= _compare_opts.size())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ITTable::SetFlags");
    }

    if (newOpts & DT_MASK)
    {
        _compare_opts[colIndex] &= 15;
        _compare_opts[colIndex] |= ((~0) & ( newOpts & DT_MASK ));
    }

    if (newOpts & CASE_INSENSE)
    {
        _compare_opts[colIndex] |= CASE_INSENSE;
    }
    else
        _compare_opts[colIndex] &= ~CASE_INSENSE;

    if (newOpts & W_SPACE_INSENSE)
        _compare_opts[colIndex] |= W_SPACE_INSENSE;
    else
        _compare_opts[colIndex] &= ~W_SPACE_INSENSE;

    ValidateOptions(colIndex);

    return NO_TABLE_ERROR;

}


void ITTable::CreateColumn(const unsigned int atColIndex,
  const vector<string>& col)
{

    CreateColumn(atColIndex, col.begin(), col.end());

}


void ITTable::CreateColumn(const unsigned int atColIndex,
  vector<string>::const_iterator colBeg, vector<string>::const_iterator colEnd)
{

    if (atColIndex > GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ITTable::CreateColumn");
    }

    if (colBeg > colEnd)
    {
        // Wrong tuple range.
        throw out_of_range("Invalid column range in TTable::CreateColumn");
    }

    unsigned int colSize = colEnd - colBeg;

    if ((GetNumRows() != 0) && (colSize > GetNumRows()))
    {
        // Wrong column index.
        throw out_of_range("Invalid column size in ITTable::FillColumn");
    }

    if (colSize != 0)
    {
        if (_orient == eCOLUMN_WISE)
            _ttable.InsertTuple(atColIndex, colBeg, colEnd);
        else
            _ttable.InsertColumn(atColIndex, colBeg, colEnd);

        if (atColIndex >= _compare_opts.size())
        {
            _compare_opts.insert(_compare_opts.end(),
                atColIndex - _compare_opts.size() + 1, DEFAULT_OPTIONS);
        }
        else
            _compare_opts.insert(_compare_opts.begin() + atColIndex,
              DEFAULT_OPTIONS);

        if (atColIndex != (GetNumColumns() - 1))
        {
            // If it is an insert, the column indices have to be updated.
            UpdateColListOnColInsert(atColIndex);
        }

    }

}


void ITTable::VerifyColumnsIndices(const vector<unsigned int>& colIndices)
{

    for (unsigned int index = 0; index < colIndices.size(); index++)
    {
        if (colIndices[index] >= GetNumColumns())
        {
            // Wrong column index.
            throw out_of_range("Invalid column index in "\
              " ITTable::VerifyColumnsIndices");
        }
    }

}


void ITTable::InsertColumn(const unsigned int colIndex,
  const vector<string>& col)
{

    InsertColumn(colIndex, col.begin(), col.end());

}


void ITTable::InsertColumn(const unsigned int colIndex,
  vector<string>::const_iterator colBeg, vector<string>::const_iterator colEnd)
{

    if (GetNumRows() != 0)
    {
        vector<string> newCol(GetNumRows(), string());
        CreateColumn(colIndex, newCol);
        FillColumn(colIndex, colBeg, colEnd);
    }
    else
    {
        CreateColumn(colIndex, colBeg, colEnd);
    }

#ifdef VLAD_IMPROVE
    // VLAD. It is better to introduce one data structure to
    // account for all column info. E.g. we could combine _compare_opts
    // and _precision to be two fields of a struct.
#endif

}


void ITTable::GetColumn(vector<string>& col, const unsigned int colIndex,
  const unsigned int indexIndex)
{

    col.clear();

    col.reserve(_indices[indexIndex].size());

    for (tIndex::iterator pos = _indices[indexIndex].begin();
      pos != _indices[indexIndex].end(); ++pos)
    {
        col.push_back(operator()(pos->second, colIndex));
    }

}
