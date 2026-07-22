//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file ISTable.C
**
** \brief Implementation file for ISTable class.
*/


#include <limits.h>

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "Exceptions.h"
#include "GenString.h"
#include "ISTable.h"


using std::find;
using std::sort;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::ofstream;
using std::ostream;
using std::cout;
using std::setw;
using std::endl;


const string ISTable::_version("V9");

const unsigned char ISTable::DEFAULT_OPTIONS = (ISTable::DT_STRING_VAL << 4);


ISTable::ISTable(const Char::eCompareType colCaseSense) :
  _orient(ITTable::eROW_WISE), _colCaseSense(colCaseSense),
  _colNames(StringLess(colCaseSense))
{

    Init();

}


ISTable::ISTable(ITTable::eOrientation orient,
  const Char::eCompareType colCaseSense) : _orient(orient),
  _colCaseSense(colCaseSense), _colNames(StringLess(colCaseSense))
{

    Init();

}


ISTable::ISTable(const string& name,
  const Char::eCompareType colCaseSense) : _name(name),
  _orient(ITTable::eROW_WISE), _colCaseSense(colCaseSense),
  _colNames(StringLess(colCaseSense))
{

    Init();

}


ISTable::ISTable(const string& name, ITTable::eOrientation orient,
  const Char::eCompareType colCaseSense) : _name(name),
  _orient(orient), _colCaseSense(colCaseSense),
  _colNames(StringLess(colCaseSense))
{

    Init();

}


ISTable::ISTable(const ISTable& inTable)
{

    _name = inTable._name;

    _ittables = inTable._ittables;

    _orient = inTable._orient;

    _colCaseSense = inTable._colCaseSense;

    _colNames = inTable._colNames;

    _precision = inTable._precision;
    _compare_opts = inTable._compare_opts;

    _indexNames = inTable._indexNames;
    _listsOfColumns = inTable._listsOfColumns;
    _unique = inTable._unique;

    _ser = inTable._ser;

    _modified = inTable._modified;

    _numRows = inTable._numRows;

    _rowIndexCache = inTable._rowIndexCache;
    _rowLocCache = inTable._rowLocCache;

}


ISTable::~ISTable()
{

    Clear();

}


ISTable& ISTable::operator=(const ISTable& inTable)
{

    if (this != &inTable)
    {
        Clear();

        _name = inTable._name;

        _ittables = inTable._ittables;

        _orient = inTable._orient;

        _colCaseSense = inTable._colCaseSense;

        _colNames = inTable._colNames;

        _precision = inTable._precision;
        _compare_opts = inTable._compare_opts;

        _indexNames = inTable._indexNames;
        _listsOfColumns = inTable._listsOfColumns;
        _unique = inTable._unique;

        _ser = inTable._ser;

        _modified = inTable._modified;

        _numRows = inTable._numRows;

        _rowIndexCache = inTable._rowIndexCache;
        _rowLocCache = inTable._rowLocCache;

    }

    return(*this);

}


ISTable::eTableDiff ISTable::operator==(ISTable& inTable)
{

    if (_colCaseSense != inTable._colCaseSense)
        return(ISTable::eCASE_SENSE);

    if (_colNames.size() > inTable._colNames.size())
        return(ISTable::eMORE_COLS);

    if (_colNames.size() < inTable._colNames.size())
        return(ISTable::eLESS_COLS);

    if (_colNames != inTable._colNames)
        return(ISTable::eCOL_NAMES);

    if (GetNumRows() > inTable.GetNumRows())
        return(ISTable::eMORE_ROWS);

    if (GetNumRows() < inTable.GetNumRows())
        return(ISTable::eLESS_ROWS);

    vector<string> firstTableCol, secondTableCol;

    for (unsigned int colI = 0; colI < _colNames.size(); ++colI)
    {
        GetColumn(firstTableCol, _colNames[colI]);
        inTable.GetColumn(secondTableCol, _colNames[colI]);
        if (firstTableCol != secondTableCol)
            return(ISTable::eCELLS);
    }

    return(ISTable::eNONE);

}


void ISTable::SetName(const string& name)
{

    _name = name;

}


const vector<string>& ISTable::GetColumnNames() const
{

    return(_colNames.get_vector());

}


bool ISTable::IsColumnPresent(const string& colName)
{

    try
    {
        GetColumnIndex(colName);

        return(true);
    }
    catch (NotFoundException& exc)
    {
        return(false);
    }

}


void ISTable::AddColumn(const string& colName, const vector<string>& col)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::AddColumn");
    }

    if ((GetNumRows() != 0) && (col.size() > GetNumRows()))
    {
        throw out_of_range("In table \"" + _name + "\", size of column \"" +
          colName + "\", " + String::IntToString(col.size()) + ", is greater "\
          "than the number of rows, " + String::IntToString(GetNumRows()) +
          ", generated at: ISTable::AddColumn");
    }

    InsertColumn(colName, _colNames.size(), col);

}


void ISTable::InsertColumn(const string& colName, const string& atColName,
  const vector<string>& col)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::InsertColumn");
    }

    if (atColName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty at column name",
          "ISTable::InsertColumn");
    }

    try
    {
        unsigned int atColIndex = GetColumnIndex(atColName);

        InsertColumn(colName, atColIndex, col);
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::InsertColumn");
        throw;
    }

}


void ISTable::FillColumn(const string& colName, const vector<string>& col)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::FillColumn");
    }

    if (col.empty())
    {
        return;
    }

    if ((GetNumRows() != 0) && (col.size() > GetNumRows()))
    {
        throw out_of_range("Column size greater than number of rows in "\
          " ISTable::FillColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        if (colIndex < _ittables[0].GetNumColumns())
        {
            // If column to be filled already exists in subtables,
            // just fill it.
            for (unsigned int tableI = 0, numRows = 0, prevNumRows = 0;
              tableI < _ittables.size(); ++tableI)
            {
                numRows += _ittables[tableI].GetNumRows();

                if (numRows > col.size())
                {
                    _ittables[tableI].FillColumn(colIndex, col.begin() +
                      prevNumRows, col.end());
                    break;
                }
                else
                {
                    _ittables[tableI].FillColumn(colIndex, col.begin() +
                      prevNumRows, col.begin() + numRows);
                }

                prevNumRows += _ittables[tableI].GetNumRows();
            }
        }
        else
        {
            // If column to be filled does not exist in subtables,
            // create it and fill it.
            CreateColumn(colIndex, col);
        }
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::FillColumn");
        throw;
    }

}


void ISTable::GetColumn(vector<string>& col, const string& colName)
{

    col.clear();

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::GetColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        vector<string> subCol;
        for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        {
            _ittables[tableI].GetColumn(subCol, colIndex);
            col.insert(col.end(), subCol.begin(), subCol.end());
        }
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::GetColumn");
        throw;
    }

}


void ISTable::GetColumn(vector<string>& col, const string& colName,
  const unsigned int fromRowIndex, unsigned int toRowIndex)
{

    col.clear();

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::GetColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        for (unsigned int rowI = fromRowIndex; rowI <= toRowIndex; ++rowI)
        {
            pair<unsigned int, unsigned int> rowLoc;
            GetRowLocation(rowLoc, rowI);

            col.push_back(_ittables[rowLoc.first](rowLoc.second, colIndex));
        }
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::GetColumn");
        throw;
    }
}


void ISTable::GetColumn(vector<string>& col, const string& colName,
  const vector<unsigned int>& rowIndex)
{

    col.clear();

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::GetColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        for (unsigned int rowI = 0; rowI < rowIndex.size(); ++rowI)
        {
            pair<unsigned int, unsigned int> rowLoc;
            GetRowLocation(rowLoc, rowIndex[rowI]);

            col.push_back(_ittables[rowLoc.first](rowLoc.second, colIndex));
        }
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::GetColumn");
        throw;
    }

}


void ISTable::RenameColumn(const string& oldColName, const string& newColName)
{

    if (oldColName.empty())
    {
        // Empty old column name.
        throw EmptyValueException("Empty old column name",
          "ISTable::RenameColumn");
    }

    if (newColName.empty())
    {
        // Empty new column name.
        throw EmptyValueException("Empty new column name",
          "ISTable::RenameColumn");
    }

    if (IsColumnPresent(newColName))
    {
        // Column already exists
        throw AlreadyExistsException("Duplicate column name",
          "ISTable::RenameColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(oldColName);

        _colNames.erase(oldColName);

        _colNames.insert(colIndex, newColName);
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Old column not found",
          "ISTable::RenameColumn");
        throw;
    }

}


void ISTable::ClearColumn(const string& colName)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::ClearColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
            _ittables[tableI].ClearColumn(colIndex);
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::ClearColumn");
        throw;
    }
}


void ISTable::DeleteColumn(const string& colName)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::DeleteColumn");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        if (GetNumRows() != 0)
        {
            for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
                _ittables[tableI].DeleteColumn(colIndex);
        }

        for (unsigned int i = 0; i < _indexNames.size(); ++i)
        {
            if (IsColumnInIndex(i, colIndex))
            {
                DeleteIndex(i);
            }
        }

        UpdateColListOnColDelete(colIndex);

        _precision.erase(_precision.begin() + colIndex);

        _compare_opts.erase(_compare_opts.begin() + colIndex);

        _colNames.erase(colName);

        if (_ittables[0].GetNumRows() == 0)
            _numRows = 0;
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::DeleteColumn");
        throw;
    }

}


unsigned int ISTable::AddRow(const vector<string>& row)
{

    return(InsertRow(GetNumRows(), row));

}


unsigned int ISTable::InsertRow(const unsigned int rowIndex,
  const vector<string>& row)
{

    if (_colNames.empty())
    {
        // Column not found.
        throw EmptyContainerException("No columns in table",
          "ISTable::InsertRow");
    }

    if ((!row.empty()) && (row.size() > _colNames.size()))
    {
        // Wrong row index.
        throw out_of_range("Invalid row size in ISTable::InsertRow");
    }

    if (GetNumRows() == 0)
    {
        for (unsigned int colI = 0; colI < _colNames.size(); ++colI)
        {
           vector<string> newCol;
           if (colI < row.size())
               newCol.push_back(row[colI]);
           else
               newCol.push_back(string());

           _ittables[0].InsertColumn(colI, newCol);
           _ittables[0].SetFlags(_compare_opts[colI], colI);
        }

        ++_numRows;
    }
    else
    {
        // If this is an row append to a non-empty table and the size of the
        // last table has reached its maximum, create a new table.
        if ((rowIndex == GetNumRows()) &&
          (_ittables[_ittables.size() - 1].GetNumRows() ==
          MAX_NUM_ITTABLE_ROWS))
        {
            // Begin create new table and set the row
            ITTable newTable(_orient);
            vector<string> newCol;
            newCol.push_back(string());

            for (unsigned int colI = 0; colI < _colNames.size(); ++colI)
            {
                newTable.InsertColumn(colI, newCol);
                newTable.SetFlags(_compare_opts[colI], colI);
            }

            for (unsigned int indI = 0; indI < _indexNames.size(); ++indI)
                newTable.CreateIndex(_listsOfColumns[indI], _unique[indI]);

            try
            {
                newTable.FillRow(0, row);

                _ittables.push_back(newTable);

                ++_numRows;
            }
            catch (AlreadyExistsException& exc)
            {
                throw;
            }
            // End of create new table and set the row
        }
        else
        {
            if (rowIndex == GetNumRows())
                // This is an add row
                _ittables[_ittables.size() - 1].AddRow(row);
            else
            {
                // This is an insert
                pair<unsigned int, unsigned int> rowLoc;
                GetRowLocation(rowLoc, rowIndex);

                _ittables[rowLoc.first].InsertRow(rowLoc.second, row);

                // VLAD - IT CAN BE SMARTER THAN THIS, BUT LEAVE IT FOR NOW
                // SINCE INSERTS ARE REALLY RARE
                CacheRowLocation(0);
            }

            ++_numRows;
        }
    }

    return(GetNumRows());

}


void ISTable::FillRow(const unsigned int rowIndex, const vector<string>& row)
{

    pair<unsigned int, unsigned int> rowLoc;
    GetRowLocation(rowLoc, rowIndex);

    _ittables[rowLoc.first].FillRow(rowLoc.second, row);

}


void ISTable::GetRow(vector<string>& row, const unsigned int rowIndex,
  const string& fromColName, const string& toColName)
{

    row.clear();

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ISTable::GetRow");
    }

    try
    {
        unsigned int fromColIndex = 0;
        if (!fromColName.empty())
        {
            fromColIndex = GetColumnIndex(fromColName);
        }

        unsigned int toColIndex = GetNumColumns() - 1;
        if (!toColName.empty())
        {
            toColIndex = GetColumnIndex(toColName);
        }

        pair<unsigned int, unsigned int> rowLoc;
        GetRowLocation(rowLoc, rowIndex);

        _ittables[rowLoc.first].GetRow(row, rowLoc.second, fromColIndex,
          toColIndex);
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::GetSubRow");
        throw;
    }

}


const vector<string>& ISTable::GetRow(const unsigned int rowIndex)
{

    if (_orient != eROW_WISE)
    {
        // Wrong table orientation.
        throw InvalidStateException("Cannot get row reference on "\
          "column-wise oriented table", "ISTable::GetRow");
    }

    if (rowIndex >= GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ISTable::GetRow");
    }

    pair<unsigned int, unsigned int> rowLoc;
    GetRowLocation(rowLoc, rowIndex);

    return(_ittables[rowLoc.first].GetRow(rowLoc.second));

}


void ISTable::ClearRow(const unsigned int rowIndex)
{

    pair<unsigned int, unsigned int> rowLoc;
    GetRowLocation(rowLoc, rowIndex);

    _ittables[rowLoc.first].ClearRow(rowLoc.second);

}


void ISTable::DeleteRow(const unsigned int rowIndex)
{

    pair<unsigned int, unsigned int> rowLoc;
    GetRowLocation(rowLoc, rowIndex);

    _ittables[rowLoc.first].DeleteRow(rowLoc.second);

    if (_ittables[rowLoc.first].GetNumRows() == 0)
        _ittables.erase(_ittables.begin() + rowLoc.first);

    --_numRows;

    if (_ittables.empty())
    {
        // Must have at least one internal table 
        ITTable firstTable(_orient);

        _ittables.push_back(firstTable);
    }

    // VLAD - IT CAN BE SMARTER THAN THIS, BUT LEAVE IT FOR NOW
    // SINCE INSERTS ARE REALLY RARE
    CacheRowLocation(0);

}


void ISTable::DeleteRows(const vector<unsigned int>& rows)
{

    // VLAD ERROR CHECKING ALL THE INDICES in rows
    // VLAD IMPROVE WITH REVERSE ITERATORS or more meaningfull index loops

    vector<unsigned int> sortedRows(rows);

    sort(sortedRows.begin(), sortedRows.end());

    for (unsigned int rowIndex = 0; rowIndex < sortedRows.size(); ++rowIndex)
    {
        DeleteRow(sortedRows[sortedRows.size() - 1 - rowIndex]);
    }

}


void ISTable::FindDuplicateRows(vector<pair<unsigned int, unsigned int> >& duplRows, const vector<string>& colNames, const bool keep,
  const eSearchDir searchDir)
{

    duplRows.clear();

    vector<unsigned int> colIndices;

    // VLAD: MAKE AND UTILIZE METHOD THAT CONVERTS FROM COLUMN NAMES
    // TO COLUMN INDICES
    // VLAD VALUE replace with GetColumnIndices()
    try
    {
        for (unsigned int i = 0; i < colNames.size(); ++i)
        {
            colIndices.push_back(GetColumnIndex(colNames[i]));
        }
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::FindDuplicateRows");
        throw;
    }

    for (unsigned int rowI = 0; rowI < GetNumRows(); rowI++)
    {
        unsigned int realRowIndex;

        if (searchDir == eFORWARD)
        {
            realRowIndex = rowI;
        }
        else
        {
            realRowIndex = GetNumRows() - rowI - 1;
        }

        bool cont = false;
        // See if this row has already been identified as duplicate and
        // ignore it
        for (unsigned int duplRowI = 0; duplRowI < duplRows.size(); ++duplRowI)
        {
            if (duplRows[duplRowI].second == realRowIndex)
            {
                cont = true;
                break;
            }
        }

        if (cont)
            continue;

        vector<string> row;

        for (unsigned int colI = 0; colI < colNames.size(); ++colI)
            row.push_back(operator()(realRowIndex, colNames[colI]));

        vector<unsigned int> res;

        if (searchDir == eFORWARD)
        {
            if (realRowIndex != (GetNumRows() - 1))
                Search(res, row, colNames, realRowIndex + 1);
        }
        else
        {
            if (realRowIndex != 0)
                Search(res, row, colNames, realRowIndex - 1, eBACKWARD);
        }

        for (unsigned int resI = 0; resI < res.size(); ++resI) 
        {
            duplRows.push_back(make_pair(realRowIndex, res[resI]));
        }
    }

    if ((!keep) && (!duplRows.empty()))
    {
        vector<unsigned int> duplRowIndices;

        for (unsigned int rowI = 0; rowI < duplRows.size(); rowI++)
        {
            duplRowIndices.push_back(duplRows[rowI].second);
        }

        DeleteRows(duplRowIndices);
    }

}


void ISTable::UpdateCell(const unsigned int rowIndex, const string& colName,
  const string& cell)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::UpdateCell");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        if (rowIndex >= GetNumRows())
        {
            // Wrong row index.
            throw out_of_range("Invalid row index in ISTable::UpdateCell");
        }

        pair<unsigned int, unsigned int> rowLoc;
        GetRowLocation(rowLoc, rowIndex);

        _ittables[rowLoc.first].UpdateCell(cell, colIndex, rowLoc.second);
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::UpdateCell");
        throw;
    }
}


const string& ISTable::operator()(const unsigned int rowIndex,
      const string& colName) const
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::operator()");
    }

    try
    {
        if (rowIndex >= GetNumRows())
        {
            // Wrong row index.
            throw out_of_range("Invalid row index in ISTable::operator()");
        }

        unsigned int colIndex = GetColumnIndex(colName);

        pair<unsigned int, unsigned int> rowLoc;
        GetRowLocation(rowLoc, rowIndex);

        return(_ittables[rowLoc.first](rowLoc.second, colIndex));
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::operator()");
        throw;
    }
}


void ISTable::CreateIndex(const string& indexName,
  const vector<string>& colNames, const unsigned int unique)
{

    //   if unique flag is set to 0, table rows with duplicate values
    //     are allowed
    //   if unique flag is set to 1, table must not have duplicate values
    //     and rows with repeated values are deleted.

    vector<unsigned int> colIndices;

    GetColumnsIndices(colIndices, colNames);

    CreateIndex(indexName, colIndices, unique);

}


void ISTable::CreateIndex(const string& indexName,
  const vector<unsigned int>& colIndices, const unsigned int unique)
{

    if (colIndices.empty())
    {
        throw EmptyValueException("Empty column indices",
          "ISTable::CreateIndex");
    }

    _ittables[0].VerifyColumnsIndices(colIndices);

    int indexIndex = FindIndex(indexName);
    if (indexIndex != -1)
    {
        throw AlreadyExistsException("Duplicate search index name",
          "ISTable::CreateIndex");
    }
    else
    {
        _indexNames.push_back(indexName);
        _listsOfColumns.push_back(colIndices);
        _unique.push_back(unique);

        for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
            _ittables[tableI].CreateIndex(colIndices, unique);
    }

}


#ifdef VLAD_REMOVE_FROM_API
void ISTable::UpdateIndex(const string& indexName, const unsigned int rowIndex)
{

    int indexIndex = FindIndex(indexName);
    if (indexIndex == -1)
    {
        // Search index not found.
        throw NotFoundException("Search index not found",
          "ISTable::UpdateIndex");
    }

    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        _ittables[tableI].UpdateIndex(indexIndex, rowIndex);

}

void ISTable::RebuildIndex(const string& indexName)
{

    int indexIndex = FindIndex(indexName);
    if (indexIndex == -1)
    {
        // Search index not found.
        throw NotFoundException("Search index not found",
          "ISTable::RebuildIndex");
    }

    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        _ittables[tableI].RebuildIndex(indexIndex);

}
#endif


void ISTable::DeleteIndex(const string& indexName)
{

    int indexIndex = FindIndex(indexName);
    if (indexIndex == -1)
    {
        // Search index not found.
        throw NotFoundException("Search index not found",
          "ISTable::DeleteIndex");
    }

    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        _ittables[tableI].DeleteIndex(indexIndex);

    DeleteIndex(indexIndex);
}


void ISTable::CreateKey(const vector<string>& colNames)
{

    vector<unsigned int> colIndices;

    GetColumnsIndices(colIndices, colNames);

    CreateKey(colIndices);

}


void ISTable::CreateKey(const vector<unsigned int>& colIndices)
{

    // VLAD HARDCODED CONST
    CreateIndex("__key__", colIndices, 1);

}


void ISTable::DeleteKey()
{

    DeleteIndex("__key__");

}


void ISTable::SetFlags(const string& colName, const unsigned char newOpts)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::SetFlags");
    }

    try {

    unsigned int colIndex = GetColumnIndex(colName);

    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        if (colIndex < _ittables[tableI].GetNumColumns())
            _ittables[tableI].SetFlags(newOpts, colIndex);

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
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::SetFlags");
        throw;
    }
}



void ISTable::ValidateOptions(unsigned int colIndex)
{

    if (colIndex >= _compare_opts.size())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ISTable::ValidateOptions");
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

unsigned char ISTable::GetDataType(const string& colName)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::GetDataType");
    }

    try {

    unsigned int colIndex = GetColumnIndex(colName);

    return((_compare_opts[colIndex] & DT_MASK) >> 4);

    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::GetDataType");
        throw;
    }
}


int ISTable::FindIndex(const string& indexName)
{

    // VLAD: FOR LOOP SEARCH: Use STL construct given below for this search
    // int index = distance(strings.begin(), find(strings.begin(),
    // strings.end(), string("xyz")));
    for (unsigned int i = 0; i < _indexNames.size(); i++)
    {
        if (_indexNames[i] == indexName)
        {
            return(i);
        }
    }

    return(-1);

}


int ISTable::FindKeyIndex()
{

    // VLAD HARDCODED CONST
    return(FindIndex("__key__"));

}



void ISTable::SetUnion(const vector<unsigned int>& a,
  const vector<unsigned int>& b, vector<unsigned int>& ret)
{

  ret.clear();

  unsigned int morea, moreb, ia, ib, lena, lenb, itema, itemb;

  ia = ib = morea = moreb = lena = lenb = 0;

  if (!a.empty())
  {
    lena  = a.size();
    morea = ia < lena;
  }

  if (morea)
    itema = a[ia];
  else
    itema = INT_MAX;

  if (!b.empty())
  {
    lenb  = b.size();
    moreb = ib < lenb;
  }

  if (moreb)
    itemb = b[ib];
  else
    itemb = INT_MAX;

  while (morea || moreb)
  {
    if (itema < itemb)
    {
      ret.push_back(itema);
      ia++; morea = ia < lena;
      if (morea)
        itema = a[ia];
      else
        itema = INT_MAX;
    }
    else if ( itema  == itemb)
    {
      ret.push_back(itema);

      ia++; morea = ia < lena;
      if (morea)
        itema = a[ia];
      else
        itema = INT_MAX;

      ib++; moreb = ib < lenb;
      if (moreb)
        itemb = b[ib];
      else
        itemb = INT_MAX;
    }
    else
    {
      ret.push_back(itemb);
      ib++; moreb = ib < lenb;
      if (moreb)
        itemb = b[ib];
      else
        itemb = INT_MAX;
    }
  }

}


void ISTable::SetIntersect(const vector<unsigned int>& a,
  const vector<unsigned int>& b, vector<unsigned int>& ret)
{

    ret.clear();

    int more, ia, ib, lena, lenb;


    //  cerr << "SetIntersect() Starting " << endl;

    if (a.empty() || b.empty())
        return;

    lena = a.size();
    lenb = b.size();

    ia = ib = 0;
    more = 1;

    while (more)
    {
        if (a[ia] < b[ib])
        {
            ia++;
            more = ia < lena;
        }
        else if (a[ia] == b[ib])
        {
            ret.push_back(a[ia]);
            ia++;
            ib++;
            more = (ib < lenb) && (ia < lena);
        }
        else
        {
            ib++;
            more = ib < lenb;
        }
    }

}


void ISTable::Init()
{

    _modified = false;

    _ser = NULL;

    _numRows = 0;

    _rowIndexCache = 0;
    _rowLocCache.first = 0;
    _rowLocCache.second = 0;

    ITTable firstTable(_orient);

    _ittables.push_back(firstTable);

}


void ISTable::Clear()
{

    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        _ittables[tableI].Clear();

    _colNames.clear();

    _modified = false;
    _numRows = 0;
    _rowIndexCache = 0;
    _rowLocCache.first = 0;
    _rowLocCache.second = 0;
    _indexNames.clear();
    _listsOfColumns.clear();
    _unique.clear();
    _precision.clear();
    _compare_opts.clear();

}

unsigned int ISTable::GetColumnIndex(const string& colName) const
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::GetColumnIndex");
    }

    unsigned int colIndex = _colNames.find(colName);

    if (colIndex == _colNames.size())
    {
        // Column not found.
        throw NotFoundException("Column \"" + colName + "\" not found in "\
          "the table \"" + _name + "\"", "ISTable::GetColumnIndex");
    }
    else
    {
        return(colIndex);
    }

}


void ISTable::GetColumnsIndices(vector<unsigned int>& colIndices,
  const vector<string>& colNames)
{

    colIndices.clear();

    try
    {
        for (unsigned int index = 0; index < colNames.size(); ++index)
        {
            colIndices.push_back(GetColumnIndex(colNames[index]));
        }
    }

    catch (RcsbException& rcsbException)
    {
        colIndices.clear();
        throw;
    }

}


string ISTable::CreateInternalIndexName(const unsigned int indexIndex)
{

    // VLAD HARCODED CONST
    string name("index_");

    name += String::IntToString(indexIndex);

    return(name);

}


void ISTable::Search(vector<unsigned int>& res, const string& target,
  const string& colName, const unsigned int fromRowIndex,
  const eSearchDir searchDir, const eSearchType searchType)
{

    res.clear();

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::Search");
    }

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        if (GetNumRows() == 0)
            return;

        vector<string> targets;
        targets.push_back(target);

        vector<unsigned int> colIds;

        colIds.push_back(colIndex);

        Search(res, targets, colIds, fromRowIndex, searchDir, searchType);
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::Search");
        throw;
    }

}


void ISTable::Search(vector<unsigned int>& res, const vector<string>& targets,
  const vector<string>& colNames, const unsigned int fromRowIndex,
  const eSearchDir searchDir, const eSearchType searchType,
  const string& indexName)
{

    res.clear();

    if (targets.size() != colNames.size())
    {
        throw out_of_range("colNames and targets have different size "\
          "in ISTable::Search");
    }

    if (GetNumRows() == 0)
        return;

    vector<unsigned int> colIds;

    GetColumnsIndices(colIds, colNames);

    Search(res, targets, colIds, fromRowIndex, searchDir, searchType,
      indexName);

}


void ISTable::Search(vector<unsigned int>& res, const vector<string>& targets,
  const vector<unsigned int>& colIds, const unsigned int fromRowIndex,
  const eSearchDir searchDir, const eSearchType searchType,
  const string& indexName)
{

    res.clear();

    int indexIndex;

    if (indexName.empty())
    {
        indexIndex = _ittables[0].FindIndex(colIds);

        if (indexIndex == -1)
        {
            // No search index found. Create a new one.
            indexIndex = _indexNames.size();
            string name = CreateInternalIndexName(indexIndex);
            CreateIndex(name, colIds);
        }
    }
    else
    {
        indexIndex = FindIndex(indexName);
        if (indexIndex == -1)
        {
            // Search index not found.
            throw NotFoundException("Search index not found",
              "ISTable::Search");
        }
    }

    vector<unsigned int> subRes;
 
    pair<unsigned int, unsigned int> rowLoc;
    GetRowLocation(rowLoc, fromRowIndex);
 
    if (searchDir == eFORWARD)
    {
        unsigned int prevNumRows = fromRowIndex - rowLoc.second;

        for (unsigned int tableI = rowLoc.first; tableI < _ittables.size();
          ++tableI)
        {
            _ittables[tableI].Search(subRes, targets, colIds,
              (unsigned int)indexIndex, searchType);

            for (unsigned int subResI = 0; subResI < subRes.size(); ++subResI)
            {
                if ((tableI == rowLoc.first) && ((prevNumRows +
                  subRes[subResI]) < fromRowIndex))
                {
                    continue;
                }

                res.push_back(prevNumRows + subRes[subResI]);
            }

            prevNumRows += _ittables[tableI].GetNumRows();
        }
    }
    else
    {
        unsigned int prevNumRows = fromRowIndex - rowLoc.second;

        for (unsigned int tableI = 0; tableI <= rowLoc.first; ++tableI)
        {
            _ittables[rowLoc.first - tableI].Search(subRes, targets, colIds,
              (unsigned int)indexIndex, searchType);

            for (unsigned int subResI = 0; subResI < subRes.size(); ++subResI)
            {
                if ((tableI == 0) && ((prevNumRows +
                  subRes[subRes.size() - 1 - subResI]) > fromRowIndex))
                {
                    continue;
                }

                res.push_back(prevNumRows + subRes[subRes.size() - 1 -
                  subResI]);
            }

            prevNumRows -= _ittables[rowLoc.first - tableI].GetNumRows();
        }
    }

}


unsigned int ISTable::FindFirst(const vector<string>& targets,
  const vector<string>& colNames, const string& indexName)
{

    if (targets.size() != colNames.size())
    {
        throw out_of_range("colNames and targets have different size "\
          "in ISTable::FindFirst");
    }

    if (GetNumRows() == 0)
        return(0);

    vector<unsigned int> colIds;

    GetColumnsIndices(colIds, colNames);

    return(FindFirst(targets, colIds, indexName));

}


unsigned int ISTable::FindFirst(const vector<string>& targets,
  const vector<unsigned int>& colIds, const string& indexName)
{

    int indexIndex;

    if (indexName.empty())
    {
        indexIndex = _ittables[0].FindIndex(colIds);
        if (indexIndex == -1)
        {
            indexIndex = _indexNames.size();
            string name = CreateInternalIndexName(indexIndex);
            CreateIndex(name, colIds);
        }
    }
    else
    {
        indexIndex = FindIndex(indexName);
        if (indexIndex == -1)
        {
            // Search index not found.
            throw NotFoundException("Search index not found",
              "ISTable::FindFirst");
        }
    }

    unsigned int subFind;
 
    for (unsigned int tableI = 0, prevNumRows = 0; tableI < _ittables.size();
      ++tableI)
    {
        subFind = _ittables[tableI].FindFirst(targets, colIds, indexIndex);

        if (subFind != _ittables[tableI].GetNumRows())
        {
            return(prevNumRows + subFind);
        }

        prevNumRows += _ittables[tableI].GetNumRows();
    }

    return(GetNumRows());

}


ISTable* ISTable::Merge(ISTable& firstTable, ISTable& secondTable,
  unsigned int typeOfMerge)
{
  unsigned int i,j;
  int index;
  string cell;


  ISTable* retTableP = new ISTable(firstTable);
  retTableP->DeleteKey();

  const vector<string>& firstTableColNames = firstTable.GetColumnNames();
  const vector<string>& secondTableColNames = secondTable.GetColumnNames();

  vector<string> newColumn;
  newColumn.insert(newColumn.end(), firstTable.GetNumRows(), string());

  for (i = 0; i < secondTableColNames.size(); ++i)
  {
    if (!retTableP->IsColumnPresent(secondTableColNames[i]))
    {
      retTableP->AddColumn(secondTableColNames[i], newColumn);
    }
  }

  for (i = 0; i < secondTable.GetNumRows(); i++)
  {
    vector<string> target;
    vector<string> newRow;

    const vector<string>& newColNames = retTableP->GetColumnNames();

    newRow.insert(newRow.end(), newColNames.size(), string());
    for (j = 0; j < newColNames.size(); j++)
    {
      if ((index = secondTable.GetColumnIndex(newColNames[j])) >= 0)
      {
	newRow[retTableP->GetColumnIndex(newColNames[j])] =
          secondTable(i, index);
      }
    }

    int firstKeyIndex = firstTable.FindKeyIndex();
    int secondKeyIndex = secondTable.FindKeyIndex();

    target.insert(target.end(),
      secondTable._listsOfColumns[secondKeyIndex].size(), string());

    // Find out what column ids are in key index
    for (unsigned int jj=0;
      jj < secondTable._listsOfColumns[secondKeyIndex].size(); jj++)
    {
      target[retTableP->GetColumnIndex(
        secondTableColNames[secondTable._listsOfColumns[secondKeyIndex][jj]])] =
        secondTable(i, secondTable._listsOfColumns[secondKeyIndex][jj]);
    }

    vector<unsigned int> result;
    result.clear();

    vector<string> searchCols;
    searchCols.insert(searchCols.end(),
      firstTable._listsOfColumns[firstKeyIndex].size(), string());

    for (unsigned int jj=0;
      jj < firstTable._listsOfColumns[firstKeyIndex].size(); jj++)
    {
      searchCols[retTableP->GetColumnIndex(
        firstTableColNames[firstTable._listsOfColumns[firstKeyIndex][jj]])] =
        firstTableColNames[firstTable._listsOfColumns[firstKeyIndex][jj]];
    }

    retTableP->Search(result, target, searchCols);
    if (result.empty())
    {
      retTableP->AddRow(newRow);
    }
    else
    {
      if (typeOfMerge == 1)
      {
	// overlap
	for (j=0; j<firstTable._colNames.size(); j++)
        {
	  cell.clear();
	  cell = firstTable.operator()(result[0], j);
	  if (cell.empty())
          {
	    if (secondTable.IsColumnPresent(firstTableColNames[j]))
            {
	      cell.clear();
	      cell = secondTable(i, firstTableColNames[j]);
	      retTableP->UpdateCell(cell,j,result[0]);
	    }
	  }
	}
      }
      else
      {
	// overwrite
        retTableP->FillRow(result[0], newRow);
      }
    }
  }

  return(retTableP);

}


ostream& operator<<(ostream& out, const ISTable& isTable)
{

    const vector<string>& colNames = isTable.GetColumnNames();

    unsigned int numRows = isTable.GetNumRows();

    out << endl << "Table " << isTable.GetName() <<  "  Columns: " <<
      colNames.size() << "  Rows: " << numRows << endl << endl;

    out << setw(5) << "RowNo";

    for (unsigned int i = 0; i < colNames.size(); i++)
    {
        // VLAD HARDCODED CONST
        out << setw(10) << colNames[i] << " ";
    }
    out << endl;

    for (unsigned int j = 0; j < numRows; j++)
    {
        out << setw(5) << j;

        for (unsigned int i = 0; i < colNames.size(); i++)
        {
            // VLAD HARDCODED CONST
            try
            {
                out << setw(10) << isTable(j, colNames[i]) << " ";
            }
            catch (out_of_range& exc)
            {
                break;
            }
        }

        out << endl;
    }

    return(out);

}


bool ISTable::PrintDiff(ISTable& inTable)
{
  bool ret = true;
  unsigned int i;
  unsigned int j;
  int index;
  vector<string> target;
  vector<unsigned int> sameCol1, sameCol2;
  vector<unsigned int> diffRow1, diffRow2;
  vector<unsigned int> sameRow1, sameRow2;
  string cell;
  string cell1, cell2;
  string Name1;
  string Name2;
  ofstream rpt;


  Name1=GetName();
  Name2=inTable.GetName();

  const vector<string>& ColNames1 = GetColumnNames();
  const vector<string>& ColNames2 = inTable.GetColumnNames();

  cout<<"** Compares table "<<Name1<<" and table "<<Name2<<" **"<<endl;

  cout<<"----------------- Columns report -----------------"<<endl;

  for (i=0; i<_colNames.size(); i++) {
    if (!inTable.IsColumnPresent(ColNames1[i])) {
      cout.width(5);
      cout<<i<<"  < ("<<ColNames1[i]<<")"<<endl;
      ret = false;
    }
    else {
      index = inTable.GetColumnIndex(ColNames1[i]);
      sameCol1.push_back(GetColumnIndex(ColNames1[i]));
      sameCol2.push_back(index);
      cout.width(5);
      cout<<i;
      cout.width(20);
      cout<<index<<"    ("<<ColNames1[i]<<")"<<endl;
    }
  }

  for (i=0; i<(inTable._colNames.size()); i++) {
    if (!IsColumnPresent(ColNames2[i])) {
      cout.width(25);
      cout<<i<<"  < ("<<ColNames2[i]<<")"<<endl;
      ret = false;
    }
  }

  cout<<endl;
  cout<<"------------------ Rows report ------------------"<<endl;
  for (i=0; i<GetNumRows(); i++) {
    // Find out what column ids are in key index
    int indexIndex = inTable.FindKeyIndex();
    for (unsigned int jj=0; jj< _listsOfColumns[indexIndex].size();
      jj++) {
      cell.clear();
      cell = operator()(i, _listsOfColumns[indexIndex][jj]);
      target.push_back(cell);
    }
    vector<unsigned int> result;
    
    vector<string> searchCol;
    for (unsigned int colI = 0; colI <
      inTable._listsOfColumns[indexIndex].size(); ++colI)
        searchCol.push_back(ColNames2[inTable._listsOfColumns[indexIndex][colI]]);

    inTable.Search(result, target, searchCol);
    if (result.empty()) {
      // value for key is different
      cout.width(5);
      cout<<i<<"  <"<<endl;
      ret = false;
    }
    else {
      // looking for other then key column for same key value
      j=0;
      cell1.clear();
      cell2.clear();
      cell1 = operator()(i, sameCol1[j]);
      cell2 = inTable(result[0], sameCol2[j]);
      while(cell1 == cell2 &&j<sameCol1.size()) {
	cell1.clear();
	cell2.clear();
	cell1 = operator()(i, sameCol1[j]);
	cell2 = inTable(result[0], sameCol2[j]);
	j++;
      }
      if (cell1 == cell2) {
	cout.width(5);
	cout<<i;
	cout.width(10);
	cout<<result[0]<<endl;
      }
      else {
	diffRow2.push_back(result[0]);
	cout.width(5);
	cout<<i<<"  *"<<endl;
	ret = false;
      }
    }
    target.clear();
  }
    for (i=0; i<diffRow2.size(); i++){
    cout.width(15);
    cout<<diffRow2[i]<<"  *"<<endl;
  }
  cout<<endl;

  return ret;
}


void ISTable::SetSerializer(Serializer* ser)
{

    _ser = ser;

}


int ISTable::WriteObject(Serializer* ser, int& size)
{

    return(WriteObjectV9(ser, size));

}


int ISTable::WriteObjectV9(Serializer* ser, int& size)
{

    unsigned int num;

    if (ser == NULL)
        return ERROR_NO_FILE_NAVIGATOR;

    UInt32 firstIndex = ser->WriteString(_version);

    UInt32 currIndex = ser->WriteString(_name);

    currIndex = ser->WriteUInt32(_ittables.size());

    int unused = 0;
    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        _ittables[tableI].Write(ser, unused);

    currIndex = ser->WriteUInt32(_colCaseSense);

    currIndex = ser->WriteUInt32(_colNames.size());

    if (_colNames.empty())
    {
        size = currIndex - firstIndex + 1;

        return (firstIndex);
    }

    currIndex = ser->WriteStrings(_colNames.get_vector());

    currIndex = ser->WriteUInt32s(_precision);

    string optsToWrite(_colNames.size(), ' ');
    for (unsigned int j = 0; j < _colNames.size(); ++j)
    {
        optsToWrite[j] = _compare_opts[j];
    }
    currIndex = ser->WriteString(optsToWrite);

    num = _indexNames.size();
    currIndex = ser->WriteUInt32(num);

    if (num != 0)
    {
        currIndex = ser->WriteStrings(_indexNames);

        currIndex = ser->WriteUInt32s(_unique);

        for (unsigned int l = 0; l < num; ++l)
        {
            currIndex = ser->WriteUInt32s(_listsOfColumns[l]);
        }
    }

    size = currIndex + firstIndex + 1;

    return (firstIndex);

}


void ISTable::Read(unsigned int indexInFile)
{
    GetObject(indexInFile, _ser);
}


int ISTable::Write()
{
    int size;

    return(WriteObject(_ser, size));
}


int ISTable::GetObject(UInt32 index, Serializer* ser) {

/*
  There are several GetObject methods to support reading tables
  from binary files saved with some of previous verson
*/
  string firstString;
  string::size_type verStringIndex;

  Clear();

  ser->ReadString(firstString, index);
  index++;

  // Reverse search of a table name for a version string

  verStringIndex = firstString.rfind(" $$$3");
  if (verStringIndex != string::npos)
  {
    // This version found. Remove the version string from the table name.
    firstString.erase(verStringIndex);
    SetName(firstString);
    return GetObjectV3(index,ser);
  }
  else
  {
    verStringIndex = firstString.rfind(" $$$2");
    if (verStringIndex != string::npos)
    {
      // This version found. Remove the version string from the table name.
      firstString.erase(verStringIndex);
      SetName(firstString);
      return GetObjectV2(index,ser);
    }
    else
    {
      // Reverse search of a table name for a version string
      verStringIndex = firstString.rfind(" $$$1");
      if (verStringIndex != string::npos)
      {
        // This version found. Remove the version string from the table name
        firstString.erase(verStringIndex);
        SetName(firstString);
        return GetObjectV1_1(index,ser);
      }
      else
      {
        if (firstString == _version)
        {
            return GetObjectV9(index,ser);
        }
        else if (firstString == "V8")
        {
            return GetObjectV8(index,ser);
        }
        else if (firstString == "V7")
        {
            return GetObjectV7(index,ser);
        }
        else if (firstString == "V6")
        {
            return GetObjectV6(index,ser);
        }
        else
        {
            SetName(firstString);
            return GetObjectV1(index,ser);
        }
      }
    }
  }
}


int ISTable::GetObjectV9(UInt32 index, Serializer* ser)
{

    if (ser == NULL)
        return ERROR_NO_FILE_NAVIGATOR;

    SInt32 ret = 0;
    UInt32 numI = 0;
    UInt32 num;

    ser->ReadString(_name, index);
    index++;

    unsigned int numTables = ser->ReadUInt32(index);
    index++;

    // Insert one less, since one has already been created at construction
    _ittables.insert(_ittables.end(), numTables - 1, ITTable());

    for (unsigned int tableI = 0; tableI < numTables; ++tableI)
    {
        index = _ittables[tableI].Read(index, ser);
    }

    // Calculate the number of rows
    _numRows = 0;
    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
    {
        _numRows += _ittables[tableI].GetNumRows();
    }

    _orient = _ittables[0].GetOrientation();

    _colCaseSense = (Char::eCompareType)ser->ReadUInt32(index);
    index++;

    unsigned int numColumns = ser->ReadUInt32(index);
    index++;

    if (numColumns == 0)
    {
        return ret;
    }

    ser->ReadStrings(_colNames.get_vector(), index);
    index++;
    _colNames.index_it();

    if (numColumns != _colNames.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    ser->ReadUInt32s(_precision, index);
    index++;

    if (numColumns != _precision.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    string oToGet;
    ser->ReadString(oToGet, index);
    index++;

    if (numColumns != oToGet.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    for (unsigned int j = 0; j < oToGet.size(); j++)
    {
        _compare_opts.push_back(oToGet[j]);
        for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
            _ittables[tableI].SetFlags(oToGet[j], j);
    }

    // Number of search indices
    num = ser->ReadUInt32(index);
    index++;

    if (num != 0)
    {
        vector<string> idxNamesToGet;
        ser->ReadStrings(idxNamesToGet, index);
        index++;
        numI = idxNamesToGet.size();

        if (numI != num)
            return INTERNAL_INCONSISTENCY_ERROR;

        vector<unsigned int> unique;
        ser->ReadUInt32s(unique, index);
        index++;

        // This version of odb file ( $$$3) and ISTables (_version == 5)
        // does not take serialized indices from the *.sdb files, but
        // re-generates them.
        vector<UInt32> listToGet;
        for (unsigned int j = 0; j < num; j++)
        {
            ser->ReadUInt32s(listToGet, index);
            index++;
            CreateIndex(idxNamesToGet[j], listToGet, unique[j]);
            listToGet.clear();
        }

    }

    return ret;
}


int ISTable::GetObjectV8(UInt32 index, Serializer* ser)
{

    if (ser == NULL)
        return ERROR_NO_FILE_NAVIGATOR;

    SInt32 ret = 0;
    UInt32 numI = 0;
    UInt32 num;

    ser->ReadString(_name, index);
    index++;

    index = _ittables[0].Read(index, ser);

    // Calculate the number of rows
    _numRows = 0;
    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
    {
        _numRows += _ittables[tableI].GetNumRows();
    }

    _orient = _ittables[0].GetOrientation();

    _colCaseSense = (Char::eCompareType)ser->ReadUInt32(index);
    index++;

    unsigned int numColumns = ser->ReadUInt32(index);
    index++;

    if (numColumns == 0)
    {
        return ret;
    }

    ser->ReadStrings(_colNames.get_vector(), index);
    index++;
    _colNames.index_it();

    if (numColumns != _colNames.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    ser->ReadUInt32s(_precision, index);
    index++;

    if (numColumns != _precision.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    string oToGet;
    ser->ReadString(oToGet, index);
    index++;

    if (numColumns != oToGet.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    for (unsigned int j = 0; j < oToGet.size(); j++)
    {
        _compare_opts.push_back(oToGet[j]);
        _ittables[0].SetFlags(oToGet[j], j);
    }

    // Number of search indices
    num = ser->ReadUInt32(index);
    index++;

    if (num != 0)
    {
        vector<string> idxNamesToGet;
        ser->ReadStrings(idxNamesToGet, index);
        index++;
        numI = idxNamesToGet.size();

        if (numI != num)
            return INTERNAL_INCONSISTENCY_ERROR;

        vector<unsigned int> unique;
        ser->ReadUInt32s(unique, index);
        index++;

        // This version of odb file ( $$$3) and ISTables (_version == 5)
        // does not take serialized indices from the *.sdb files, but
        // re-generates them.
        vector<UInt32> listToGet;
        for (unsigned int j = 0; j < num; j++)
        {
            ser->ReadUInt32s(listToGet, index);
            index++;
            CreateIndex(idxNamesToGet[j], listToGet, unique[j]);
            listToGet.clear();
        }

    }

    return ret;
}


int ISTable::GetObjectV7(UInt32 index, Serializer* ser)
{

    if (ser == NULL)
        return ERROR_NO_FILE_NAVIGATOR;

    SInt32 ret = 0;
    UInt32 numI = 0;
    UInt32 num;

    ser->ReadString(_name, index);
    index++;

    index = _ittables[0].Read(index, ser);

    // Calculate the number of rows
    _numRows = 0;
    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
    {
        _numRows += _ittables[tableI].GetNumRows();
    }

    _orient = _ittables[0].GetOrientation();

    _colCaseSense = (Char::eCompareType)ser->ReadUInt32(index);
    index++;

    unsigned int numColumns = ser->ReadUInt32(index);
    index++;

    if (numColumns == 0)
    {
        return ret;
    }

    ser->ReadStrings(_colNames.get_vector(), index);
    index++;
    _colNames.index_it();

    if (numColumns != _colNames.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    ser->ReadUInt32s(_precision, index);
    index++;

    if (numColumns != _precision.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    string oToGet;
    ser->ReadString(oToGet, index);
    index++;

    if (numColumns != oToGet.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    for (unsigned int j = 0; j < oToGet.size(); j++)
    {
        _compare_opts.push_back(oToGet[j]);
        _ittables[0].SetFlags(oToGet[j], j);
    }

    // Number of search indices
    num = ser->ReadUInt32(index);
    index++;

    if (num != 0)
    {
        vector<string> idxNamesToGet;
        ser->ReadStrings(idxNamesToGet, index);
        index++;
        numI = idxNamesToGet.size();

        if (numI != num)
            return INTERNAL_INCONSISTENCY_ERROR;

        // This version of odb file ( $$$3) and ISTables (_version == 5)
        // does not take serialized indices from the *.sdb files, but
        // re-generates them.
        vector<vector<UInt32> > listToGet;
        vector<UInt32> oneList;
        for (unsigned int j = 0; j < num; j++)
        {
            listToGet.push_back(oneList);
            ser->ReadUInt32s(listToGet[j], index);
            index++;
        }

        vector<unsigned int> unique;
        ser->ReadUInt32s(unique, index);
        index++;

        for (unsigned int j = 0; j < num; j++)
        {
            CreateIndex(idxNamesToGet[j], listToGet[j], unique[j]);
        }
    }

    return ret;
}


int ISTable::GetObjectV6(UInt32 index, Serializer* ser)
{

    if (ser == NULL)
        return ERROR_NO_FILE_NAVIGATOR;

    SInt32 ret = 0;
    UInt32 numI = 0;
    UInt32 num;

    ser->ReadString(_name, index);
    index++;

    unsigned int numColumns = ser->ReadUInt32(index);
    index++;

    if (numColumns == 0)
    {
        return ret;
    }

    // unsigned int numRows = ser->ReadUInt32(index);
    // This is not used and is not read, only skipped
    index++;

    ser->ReadStrings(_colNames.get_vector(), index);
    index++;
    _colNames.index_it();

    if (numColumns != _colNames.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    ser->ReadUInt32s(_precision, index);
    index++;

    if (numColumns != _precision.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    string oToGet;
    ser->ReadString(oToGet, index);
    index++;

    if (numColumns != oToGet.size())
        return INTERNAL_INCONSISTENCY_ERROR;

    for (unsigned int j = 0; j < oToGet.size(); j++)
    {
        _compare_opts.push_back(oToGet[j]);
    }

    for (unsigned int k = 0; k < _colNames.size(); ++k)
    {
        vector<string> tmpCol;
        ser->ReadStrings(tmpCol, index);
        index++;
        _ittables[0].InsertColumn(k, tmpCol);
    }

    // Calculate the number of rows
    _numRows = 0;
    for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
    {
        _numRows += _ittables[tableI].GetNumRows();
    }

    for (unsigned int j = 0; j < oToGet.size(); j++)
    {
        _ittables[0].SetFlags(oToGet[j], j);
    }

    // Number of search indices
    num = ser->ReadUInt32(index);
    index++;

    if (num != 0)
    {
        vector<string> idxNamesToGet;
        ser->ReadStrings(idxNamesToGet, index);
        index++;
        numI = idxNamesToGet.size();

        if (numI != num)
            return INTERNAL_INCONSISTENCY_ERROR;

        // This version of odb file ( $$$3) and ISTables (_version == 5)
        // does not take serialized indices from the *.sdb files, but
        // re-generates them.
        vector<vector<UInt32> > listToGet;
        vector<UInt32> oneList;
        for (unsigned int j = 0; j < num; j++)
        {
            listToGet.push_back(oneList);
            ser->ReadUInt32s(listToGet[j], index);
            index++;
        }

        vector<unsigned int> unique;
        ser->ReadUInt32s(unique, index);
        index++;

        for (unsigned int j = 0; j < num; j++)
        {
            CreateIndex(idxNamesToGet[j], listToGet[j], unique[j]);
        }
    }

    return ret;
}


int ISTable::GetObjectV3(UInt32 index, Serializer* ser)
{

  if (ser == NULL)
      return ERROR_NO_FILE_NAVIGATOR;

  SInt32 ret = 0;
  UInt32 numI = 0;
  UInt32 num;

  ser->ReadUInt32(index);
  index++;

  unsigned int numColumns = ser->ReadUInt32(index);
  index++;

  // unsigned int numRows = ser->ReadUInt32(index);
  // This is not used and is not read, only skipped
  index++;

  ser->ReadStrings(_colNames.get_vector(), index);
  index++;
  _colNames.index_it();

  if (numColumns != _colNames.size())
      return INTERNAL_INCONSISTENCY_ERROR;

  ser->ReadUInt32s(_precision, index);
  index++;

  string oToGet;
  ser->ReadString(oToGet, index);
  index++;
  if (!_colNames.empty())
  {
    for (unsigned int j = 0; j < _colNames.size(); j++)
    {
      _compare_opts.push_back(oToGet[j]);
    }
  }

  for (unsigned int k = 0; k < _colNames.size(); ++k)
  {
    vector<string> tmpCol;
    ser->ReadStrings(tmpCol, index);
    index++;
    _ittables[0].InsertColumn(k, tmpCol);
  }

  // Calculate the number of rows
  _numRows = 0;
  for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
  {
      _numRows += _ittables[tableI].GetNumRows();
  }

  for (unsigned int j = 0; j < _colNames.size(); j++)
  {
      _ittables[0].SetFlags(oToGet[j], j);
  }
  // EnlargeRowMap(numRows);

  // Number of search indices
  num = ser->ReadUInt32(index);
  index++;

  vector<string> idxNamesToGet;
  ser->ReadStrings(idxNamesToGet, index);
  index++;
  numI = idxNamesToGet.size();

  if (numI != num)
      return INTERNAL_INCONSISTENCY_ERROR;

  vector<unsigned int> unique;
  ser->ReadUInt32s(unique, index);
  index++;

  // This version of odb file ( $$$3) and ISTables (_version == 5)
  // does not take serialized indices from the *.sdb files, but
  // re-generates them.
  vector<UInt32> listToGet;
  for (unsigned int j = 0; j < num; j++)
  {
    ser->ReadUInt32s(listToGet, index);
    index++;
    CreateIndex(idxNamesToGet[j], listToGet, unique[j]);
    listToGet.clear();
  }

  return ret;
}


int ISTable::GetObjectV2(UInt32 index, Serializer* ser) {

  if (!ser) return ERROR_NO_FILE_NAVIGATOR;
  unsigned int j;
  SInt32 ret = 0;
  UInt32 numC = 0;
  UInt32 numI = 0;
  UInt32 num;

  unsigned int version = ser->ReadUInt32(index); index++;

  if ((version==3) || (version==4))
  {
    num = ser->ReadUInt32(index);index++;
    vector<UInt32> kToGet;
    ser->ReadUInt32s(kToGet, index);index++;
    kToGet.clear();
  }

  unsigned int numColumns = ser->ReadUInt32(index); index++;
  unsigned int numRows    = ser->ReadUInt32(index); index++;
  //EnlargeRowMap(numRows);
  unsigned int numDels = ser->ReadUInt32(index); index++;

  // _colAlloc field is next. It is ignored since it is not used with
  // this version of ISTable
  index++;

  vector<string> colNamesToGet;
  ser->ReadStrings(colNamesToGet, index);
  numC = colNamesToGet.size();
  index++;

  if (numC != numColumns) return INTERNAL_INCONSISTENCY_ERROR;
  if (numC > 0) {
    string  temp;
    for (j = 0; j < numColumns; j++) {
      if (colNamesToGet[j].empty()) {
	  temp = "";
      } else {
	  temp = colNamesToGet[j];
      }
      _colNames.push_back(temp);
    }
  }


  vector<UInt32> pToGet;
  ser->ReadUInt32s(pToGet, index);
  index++;
  if (!pToGet.empty()) {
    for (j = 0; j < _colNames.size(); j++) {
      _precision.push_back(pToGet[j]);
    }
  }
  if (!pToGet.empty()) pToGet.clear();

  string oToGet;
  ser->ReadString(oToGet, index); index++;
  if (_colNames.size() > 0) {
    for (j = 0; j < _colNames.size(); j++) {
      _compare_opts.push_back(oToGet[j]);
    }
  }

  // VLAD: THIS CODE HAS NOT BEEN TESTED
  // NOTE: The data read here is a former ISTable member _deleted. Values
  // for _deleted used to be 0 or 1, where 0 indicated that row is not
  // deleted and 1 indicated that row was deleted.

  // Out of order read of _deleted, since it is needed for proper
  // construction of _data. Index is not advanced.
  vector<unsigned int> deleted;
  vector<UInt32> dToGet;
  ser->ReadUInt32s(dToGet, index + _colNames.size());
  if (!dToGet.empty()) {
    for (j = 0; j < numRows; j++) {
      deleted.push_back(dToGet[j]);
    }
  }
  if (!dToGet.empty()) dToGet.clear();

  unsigned int k, l;
  vector<string> tempColumns;
  for (k = 0; k < _colNames.size(); k++) {
    vector<string> tmpCol;
    ser->ReadStrings(tempColumns, index);
    numC = tempColumns.size();
    index++;
    if (numC > 0) {
      string * tempString = new string[numC];
      for (l = 0; l < numC; l++) {
        if (tempColumns[l].empty()) {
          tempString[l] = "";
        }
        else {
          tempString[l] = tempColumns[l];
        }
          if (!deleted[l])
              tmpCol.push_back(tempString[l]);
      }

      _ittables[0].InsertColumn(k, tmpCol);

      if (!tempColumns.empty()) {
        for (unsigned int l3 = 0; l3 < numC; l3++) {
        }
        tempColumns.clear();
      }
      delete[] tempString;
    }
  }

  index++;

  // Calculate the number of rows
  _numRows = 0;
  for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
  {
      _numRows += _ittables[tableI].GetNumRows();
  }

  for (j = 0; j < _colNames.size(); j++)
  {
        _ittables[0].SetFlags(oToGet[j], j);
  }

  num = ser->ReadUInt32(index); index++;

  vector<string> idxNamesToGet;
  ser->ReadStrings(idxNamesToGet, index);
  numI = idxNamesToGet.size();
  index++;

  if (numI != num) return INTERNAL_INCONSISTENCY_ERROR;
  if (numI > 0) {
    string  temp;
    for (j = 0; j < numI; j++) {
      if (idxNamesToGet[j].empty()) {
	  temp = "";
      } else {
	  temp = idxNamesToGet[j];
      }
        //_indexNames.push_back(temp);
    }
  }

  vector<UInt32> uniqueToGet;
  ser->ReadUInt32s(uniqueToGet, index);
  index++;

  // This version of odb file ( $$$3) and ISTables (_version == 5)
  // does not take serialized indices from the *.sdb files, but
  // re-generates them.
  vector<UInt32> listToGet;
  for (j=0; j<num; j++)
  {
    ser->ReadUInt32s(listToGet, index);
    index++;
    vector<unsigned int> list;
    for (unsigned int l2=0; l2<listToGet.size(); l2++){
      list.push_back(listToGet[l2]);
    }

    if (j < uniqueToGet.size())
        CreateIndex(idxNamesToGet[j], list, uniqueToGet[j]);
    else
        CreateIndex(idxNamesToGet[j], list);

    //_listsOfColumns.push_back(list);
    list.clear();
    if (!listToGet.empty()) listToGet.clear();
  }

   // Correct the members
  numRows -= numDels;

  // EnlargeRowMap(numRows);

  return ret;
}


int ISTable::GetObjectV1_1(UInt32 index, Serializer* ser) {

  if (!ser) return ERROR_NO_FILE_NAVIGATOR;
  unsigned int j;
  SInt32 ret = 0;
  UInt32 numC = 0;

  ser->ReadUInt32(index);
  index++;

  vector<UInt32> kToGet;
  ser->ReadUInt32s(kToGet, index);
  index++;
  kToGet.clear();

  unsigned int numColumns = ser->ReadUInt32(index);
  index++;
  // unsigned int numRows    = ser->ReadUInt32(index);
  // This is not used and is not read, only skipped
  index++;

  // _colAlloc field is next. It is ignored since it is not used with
  // this version of ISTable
  // This is not used and is not read, only skipped
  index++;

  // unsigned int treeAlloc  = (int) ser->ReadUInt32(index);
  // This is not used and is not read, only skipped
  index++;

  unsigned int numTrees = ser->ReadUInt32(index); index++;

  vector<UInt32> cMap;
  ser->ReadUInt32s(cMap, index);index++;
  cMap.clear();


  vector<string> colNamesToGet;
  ser->ReadStrings(colNamesToGet, index);
  numC = colNamesToGet.size();
  index++;

  if (numC != numColumns) return INTERNAL_INCONSISTENCY_ERROR;
  if (numC > 0) {
    string  temp;
    for (j = 0; j < numColumns; j++) {
      if (colNamesToGet[j].empty()) {
	  temp = "";
      } else {
	  temp = colNamesToGet[j];
      }
        _colNames.push_back(temp);
    }
  }
  if (!colNamesToGet.empty()) {
    if (_colNames.size() > 0) colNamesToGet.clear();
  }

  vector<UInt32> pToGet;
  ser->ReadUInt32s(pToGet, index);
  index++;
  if (!pToGet.empty()) {
    for (j = 0; j < _colNames.size(); j++) {
      _precision.push_back(pToGet[j]);
    }
  }
  if (!pToGet.empty()) pToGet.clear();

  string oToGet;
  ser->ReadString(oToGet, index); index++;
  if (_colNames.size() > 0) {
    for (j = 0; j < _colNames.size(); j++) {
      _compare_opts.push_back(oToGet[j]);
    }
  }

  unsigned int k, l;
  vector<string> tempColumns;
  for (k = 0; k < _colNames.size(); k++) {
    vector<string> tmpCol;
    ser->ReadStrings(tempColumns, index);
    numC = tempColumns.size();
    index++;
    if (numC > 0) {
      string * tempString = new string[numC];
      for (l = 0; l < numC; l++) {
        if (tempColumns[l].empty()) {
          tempString[l] = "";
        }
        else {
          tempString[l] = tempColumns[l];
        }
          tmpCol.push_back(tempString[l]);
      }

      _ittables[0].InsertColumn(k, tmpCol);

      if (!tempColumns.empty()) {
        tempColumns.clear();
      }
      delete[] tempString;
    }
  }

  // Calculate the number of rows
  _numRows = 0;
  for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
  {
      _numRows += _ittables[tableI].GetNumRows();
  }

  for (j = 0; j < _colNames.size(); j++)
  {
      _ittables[0].SetFlags(oToGet[j], j);
  }
/********************************************/
  int tableIndex = index ;

  vector<UInt32> tInts;

  for (k = 0; k < numTrees; k++) {
    ser->ReadUInt32s(tInts, tableIndex + k); index++;
    tInts.clear();
  }

  return ret;

}



int ISTable::GetObjectV1(UInt32 index, Serializer* ser) {

  if (!ser) return ERROR_NO_FILE_NAVIGATOR;
  unsigned int j;
  SInt32 ret = 0;
  UInt32 numC = 0;

  unsigned int numColumns = ser->ReadUInt32(index);
  index++;

  // unsigned int numRows    = ser->ReadUInt32(index);
  // This is not used and is not read, only skipped
  index++;

  // _colAlloc field is next. It is ignored since it is not used with
  // this version of ISTable
  index++;

  // treeAlloc is not used, is not read, only skipped
  // unsigned int treeAlloc  = (int) ser->ReadUInt32(index);
  index++;

  unsigned int numTrees = ser->ReadUInt32(index);
  index++;

  vector<UInt32> cMap;
  ser->ReadUInt32s(cMap, index);index++;
  cMap.clear();


  vector<string> colNamesToGet;
  ser->ReadStrings(colNamesToGet, index);
  numC = colNamesToGet.size();
  index++;

  if (numC != numColumns) return INTERNAL_INCONSISTENCY_ERROR;
  if (numC > 0) {
    string  temp;
    for (j = 0; j < numColumns; j++) {
      if (colNamesToGet[j].empty()) {
	  temp = "";
      } else {
	  temp = colNamesToGet[j];
      }
        _colNames.push_back(temp);
    }
  }
  if (!colNamesToGet.empty()) {
    if (_colNames.size() > 0) colNamesToGet.clear();
  }

  vector<UInt32> pToGet;
  ser->ReadUInt32s(pToGet, index);
  index++;
  if (!pToGet.empty()) {
    for (j = 0; j < _colNames.size(); j++) {
      _precision.push_back(pToGet[j]);
    }
  }
  if (!pToGet.empty()) pToGet.clear();

  string oToGet;
  ser->ReadString(oToGet, index); index++;
  if (_colNames.size() > 0) {
    for (j = 0; j < _colNames.size(); j++) {
      _compare_opts.push_back(oToGet[j]);
    }
  }

  unsigned int k, l;
  vector<string> tempColumns;
  for (k = 0; k < _colNames.size(); k++) {
    vector<string> tmpCol;
    ser->ReadStrings(tempColumns, index);
    numC = tempColumns.size();
    index++;
    if (numC > 0) {
      string * tempString = new string[numC];
      for (l = 0; l < numC; l++) {
        if (tempColumns[l].empty()) {
          tempString[l] = "";
        }
        else {
          tempString[l] = tempColumns[l];
        }
          tmpCol.push_back(tempString[l]);
      }

      _ittables[0].InsertColumn(k, tmpCol);

      if (!tempColumns.empty()) {
        tempColumns.clear();
      }
      delete[] tempString;
    }
  }

  // Calculate the number of rows
  _numRows = 0;
  for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
  {
      _numRows += _ittables[tableI].GetNumRows();
  }

  for (j = 0; j < _colNames.size(); j++)
  {
      _ittables[0].SetFlags(oToGet[j], j);
  }
/********************************************/
  int tableIndex = index ;

  vector<UInt32> tInts;

  for (k = 0; k < numTrees; k++) {
    ser->ReadUInt32s(tInts, tableIndex + k); index++;
    tInts.clear();
  }

  return ret;
}


int ISTable::UpdateCell(const string& cell, const unsigned int colIndex,
  const unsigned int rowIndex)
{

  pair<unsigned int, unsigned int> rowLoc;
  GetRowLocation(rowLoc, rowIndex);

  return(_ittables[rowLoc.first].UpdateCell(cell, colIndex, rowLoc.second));

}


const string& ISTable::operator()(const unsigned int rowIndex,
  const unsigned int colIndex) const
{

  pair<unsigned int, unsigned int> rowLoc;
  GetRowLocation(rowLoc, rowIndex);

  return(_ittables[rowLoc.first](rowLoc.second, colIndex));

}


void ISTable::CreateColumn(const string& colName, const unsigned int atColIndex,
  const vector<string>& col)
{

    if (atColIndex > _colNames.size())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in ISTable::CreateColumn");
    }

    if ((GetNumRows() != 0) && (col.size() > GetNumRows()))
    {
        // Wrong column index.
        throw out_of_range("Invalid column size in ISTable::FillColumn");
    }

    if (IsColumnPresent(colName))
    {
        // Column already exists
        throw AlreadyExistsException("Duplicate column name \"" + colName +
          "\" in table \"" + _name + "\"", "ISTable::CreateColumn");
    }

    _colNames.insert(atColIndex, colName);

    // VLAD. It is better to introduce one data structure to
    // account for all column info. E.g. we could combine _compare_opts
    // and _precision to be two fields of a struct.
    _compare_opts.insert(_compare_opts.begin() + atColIndex,
      DEFAULT_OPTIONS);

    // From ANSI C++ standard point of view, the re-casting of
    // DEFAULT_PRECISION to (unsigned int) is not necessary, since its type
    // is unsigned int. However, the linker issues unresolved symbol error
    // if this re-casting is omitted.
    _precision.insert(_precision.begin() + atColIndex,
      (unsigned int)DEFAULT_PRECISION);

    if (!col.empty())
    {
        // If column is not empty, create it in subtables and fill it.
        CreateColumn(atColIndex, col);
    }

}


void ISTable::InsertColumn(const string& colName,
  const unsigned int atColIndex, const vector<string>& col)
{

    if (GetNumRows() != 0)
    {
        if (col.size() > GetNumRows())
        {
            // Wrong column index.
            throw out_of_range("Invalid column size in ISTable::InsertColumn");
        }

        vector<string> newCol(GetNumRows(), string());
        CreateColumn(colName, atColIndex, newCol);
        FillColumn(colName, col);
    }
    else
    {
        CreateColumn(colName, atColIndex, col);
    }

}


void ISTable::GetColumn(vector<string>& col, const string& colName,
  const string& indexName)
{

    if (colName.empty())
    {
        // Empty column name.
        throw EmptyValueException("Empty column name",
          "ISTable::GetColumn");
    }

    col.clear();

    try
    {
        unsigned int colIndex = GetColumnIndex(colName);

        unsigned int indexIndex = FindIndex(indexName);

        vector<string> subCol;
        for (unsigned int tableI = 0; tableI < _ittables.size(); ++tableI)
        {
            _ittables[tableI].GetColumn(subCol, colIndex, indexIndex);
            col.insert(col.end(), subCol.begin(), subCol.end());
        }
    }
    catch (NotFoundException& notFound)
    {
        notFound.AppendMessage("Column not found",
          "ISTable::GetColumn");
        throw;
    }

}


void ISTable::GetRowLocation(pair<unsigned int, unsigned int>& rowLoc,
  const unsigned int rowIndex) const
{

    if (rowIndex != _rowIndexCache)
        CacheRowLocation(rowIndex);

    rowLoc = _rowLocCache;

}


void ISTable::CacheRowLocation(const unsigned int rowIndex) const
{

    if (rowIndex > GetNumRows())
    {
        // Wrong row index.
        throw out_of_range("Invalid row index in ISTable::CacheRowLocation");
    }

    if (rowIndex == 0)
    {
        _rowIndexCache = 0;
        _rowLocCache.first = 0;
        _rowLocCache.second = 0;

        return;
    }

    if (rowIndex > _rowIndexCache)
    {
        unsigned int currNumRows = _rowIndexCache - _rowLocCache.second;

        for (unsigned int tableI = _rowLocCache.first;
          tableI < _ittables.size(); ++tableI)
        {
            currNumRows += _ittables[tableI].GetNumRows();

            if (rowIndex < currNumRows)
            {
                _rowLocCache.first = tableI;
                _rowLocCache.second = _ittables[tableI].GetNumRows() -
                  (currNumRows - rowIndex);

                _rowIndexCache = rowIndex;

                return;
            }
        }
    }

    if (rowIndex < _rowIndexCache)
    {
        //unsigned int currNumRows = _rowIndexCache - _rowLocCache.second;
        unsigned int currNumRows = _rowIndexCache +
          (_ittables[_rowLocCache.first].GetNumRows() - _rowLocCache.second);

        for (unsigned int tableI = 0; tableI <= _rowLocCache.first; ++tableI)
        {
            currNumRows -= _ittables[_rowLocCache.first - tableI].GetNumRows();

            if (rowIndex >= currNumRows)
            {
                _rowLocCache.first = _rowLocCache.first - tableI;
                _rowLocCache.second = rowIndex - currNumRows;

                _rowIndexCache = rowIndex;

                return;
            }

        }
    }

}


void ISTable::CreateSubtables(const unsigned int numRows)
{

    // This method creates the needed number of ITTables to hold
    // given number of rows.

    // One ITTable was created at construction time.
    unsigned int currNumRows = MAX_NUM_ITTABLE_ROWS;

    // Create the rest of the ITTables, if needed.
    while (currNumRows < numRows)
    {
        _ittables.push_back(ITTable(_orient));
        currNumRows += MAX_NUM_ITTABLE_ROWS;
    }

}


void ISTable::CreateSubtableColumns(const unsigned int atColIndex, 
  const vector<string>& col)
{

    // This method creates new columns in ITTables and appropriatelly
    // fills them

    for (unsigned int tableI = 0, numRows = 0, prevNumRows = 0;
      tableI < _ittables.size(); ++tableI)
    {
        if (_ittables[tableI].GetNumRows() == 0)
            numRows += MAX_NUM_ITTABLE_ROWS;
        else
            numRows += _ittables[tableI].GetNumRows();

        if (numRows > col.size())
        {
            _ittables[tableI].InsertColumn(atColIndex, col.begin() +
              prevNumRows, col.end());
            break;
        }
        else
        {
            _ittables[tableI].InsertColumn(atColIndex, col.begin() +
              prevNumRows, col.begin() + numRows);
        }

        prevNumRows += _ittables[tableI].GetNumRows();

        _ittables[tableI].SetFlags(_compare_opts[atColIndex],
          atColIndex);
    }

}


void ISTable::CreateColumn(const unsigned int atColIndex, 
  const vector<string>& col)
{

    if (GetNumRows() == 0)
    {
        // If no rows exist, create the subtables first so that whole
        // column can fit.
        CreateSubtables(col.size());
    }

    CreateSubtableColumns(atColIndex, col);

    if (GetNumRows() == 0)
    {
        // If no rows existed prior to this, set the number of rows.
        _numRows = col.size();
    }

}


void ISTable::DeleteIndex(const unsigned int indexIndex)
{
    _unique.erase(_unique.begin() + indexIndex);
    _listsOfColumns.erase(_listsOfColumns.begin() + indexIndex);
    _indexNames.erase(_indexNames.begin() + indexIndex);
}


void ISTable::UpdateColListOnColDelete(const unsigned int colIndex)
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


bool ISTable::IsColumnInIndex(const unsigned int indexIndex,
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

